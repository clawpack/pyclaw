#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing boxclaw.geometry.
"""

from clawpack import pyclaw
from clawpack.pyclaw import geometry as pyclaw_geometry

class Patch(pyclaw_geometry.Patch):
    """Parallel Patch class."""

    __doc__ += pyclaw.util.add_parent_doc(pyclaw_geometry.Patch)

    def __init__(self,dimensions):

        import boxlib
        import numpy as np

        super(Patch,self).__init__(dimensions)

        self._ba, box   = self._create_boxarray()

        is_per = np.asarray(self.num_dim * [ 1 ], np.int32)
        rb = boxlib.RealBox(boxlib.lo(box), boxlib.hi(box))
        self._geom = boxlib.bl[self.num_dim].Geometry(box, rb, 0, is_per)

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
            lower = self.lower_global[i] + lo[i]*self.delta[i]
            upper = self.lower_global[i] + (hi[i]+1)*self.delta[i]
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

        import boxlib
        import numpy as np

        dx = self.delta
        lg = self.lower_global
        nc = self.num_cells_global

        lo = [ int(lg[d]/dx[d]) for d in range(self.num_dim) ]
        hi = [ lo[d]+nc[d]-1    for d in range(self.num_dim) ]

        box = boxlib.Box(lo, hi)
        ba  = boxlib.BoxArray([box])

        # XXX: breaking up the box array to distribute it should be
        # done more intelligently

        # boxlib supports more than one box per processor; however,
        # for pyclaw we will enforce one box per processor
        nproc_per_dim = int(float(boxlib.size())**(1./self.num_dim))

        max_sizes = [ (hi[d] - lo[d]) / nproc_per_dim + 1 for d in range(self.num_dim) ]
        max_sizes = boxlib.bl[self.num_dim].IntVect(*max_sizes)
        ba.maxSize(max_sizes)

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


