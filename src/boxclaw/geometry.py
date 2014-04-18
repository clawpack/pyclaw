#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing boxclaw.geometry.
"""

from clawpack import pyclaw
from clawpack.pyclaw import geometry as pyclaw_geometry

import logging
logger = logging.getLogger('root')

class Patch(pyclaw_geometry.Patch):
    """Parallel Patch class."""

    __doc__ += pyclaw.util.add_parent_doc(pyclaw_geometry.Patch)

    def __deepcopy__(self,memo):
        # don't recreate boxarrays
        return self

    def __init__(self,dimensions):

        import boxlib
        import numpy as np

        super(Patch,self).__init__(dimensions)

        self._ba, self._gbox = self._create_boxarray()

        is_per = np.asarray(self.num_dim * [ 1 ], np.int32)
        rb     = boxlib.RealBox(boxlib.lo(self._gbox), boxlib.hi(self._gbox))
        self._geom = boxlib.bl[self.num_dim].Geometry(self._gbox, rb, 0, is_per)

        # XXX: create a multifab from the boxarray to get geometry information
        tmp = boxlib.MultiFab(self._ba)

        fab = None
        for i in range(tmp.size()):
            fab = tmp[i]
            if fab is not None:
                break

        assert(fab is not None)

        lo = boxlib.lo(fab.box())
        hi = boxlib.hi(fab.box())
        grid_dimensions = []
        for i in range(self.num_dim):
            lower = (lo[i]+0) * self.delta[i]
            upper = (hi[i]+1) * self.delta[i]
            num_cells = hi[i]-lo[i]+1

            grid_dimensions.append(pyclaw_geometry.Dimension(lower,upper,
                                        num_cells,name=dimensions[i].name))

            if lower == self.lower_global[i]:
                grid_dimensions[-1].on_lower_boundary = True
            else:
                grid_dimensions[-1].on_lower_boundary = False

            if upper == self.upper_global[i]:
                grid_dimensions[-1].on_upper_boundary = True
            else:
                grid_dimensions[-1].on_upper_boundary = False

        self.grid = pyclaw_geometry.Grid(grid_dimensions)

        del tmp


    def _create_boxarray(self):
        """Returns a BoxLib BoxArray."""

        # note that boxlib supports more than one box per processor

        import boxlib
        import math
        import numpy as np

        dx = self.delta
        lg = self.lower_global
        nc = self.num_cells_global

        lo = [ int(lg[d]/dx[d]) for d in range(self.num_dim) ]
        hi = [ lo[d]+nc[d]-1    for d in range(self.num_dim) ]

        box = boxlib.Box(lo, hi)
        ba  = boxlib.BoxArray([box])

        max_sizes = self.num_dim * [ 1 ]

        nprocs = boxlib.size()
        if nprocs % 2**self.num_dim == 0:
            # divide domain into cubes
            nproc_per_dim = nprocs / 2**self.num_dim + 1
            print nproc_per_dim
            for d in range(self.num_dim):
                max_sizes[d] = int(math.ceil(float(hi[d] - lo[d] + 1) / nproc_per_dim))
        else:
            # divide domain into slabs
            max_sizes[0] = int(math.ceil(float(hi[0] - lo[0] + 1) / nprocs))
            for d in range(1, self.num_dim):
                max_sizes[d] = hi[d] - lo[d] + 1

        max_sizes = boxlib.bl[self.num_dim].IntVect(*max_sizes)
        ba.maxSize(max_sizes)

        if boxlib.rank() == 0:
            logger.info("max_sizes: " + str(max_sizes))
            logger.info("boxarray:  " + str(ba))

        return ba, box


class Domain(pyclaw_geometry.Domain):
    """Parallel Domain Class"""

    __doc__ += pyclaw.util.add_parent_doc(pyclaw.ClawSolver2D)

    def __init__(self,geom):
        if not isinstance(geom,list):
            geom = [geom]
        if isinstance(geom[0],Patch):
            self.patches = geom
        elif isinstance(geom[0],pyclaw_geometry.Dimension):
            self.patches = [Patch(geom)]
