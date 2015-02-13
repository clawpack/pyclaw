#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing BoxClaw geometry.
"""

import fboxlib
import math

from clawpack import pyclaw
from clawpack.pyclaw import geometry as pyclaw_geometry

import logging
logger = logging.getLogger('root')

class Patch(pyclaw_geometry.Patch):
    """Parallel Patch class."""

    __doc__ += pyclaw.util.add_parent_doc(pyclaw_geometry.Patch)

    def __deepcopy__(self,memo):
        return self

    def __init__(self,dimensions):

        super(Patch,self).__init__(dimensions)

        self._la = self._create_layout()

        bxs = self._la.local_boxes
        lo, hi = self._la.get_box(bxs[0])

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


    def _create_layout(self):
        """Returns a FBoxLib layout."""

        dx = self.delta
        lg = self.lower_global
        nc = self.num_cells_global

        lo = [ int(lg[d]/dx[d]) for d in range(self.num_dim) ]
        hi = [ lo[d]+nc[d]-1    for d in range(self.num_dim) ]

        max_sizes = self.num_dim * [ 1 ]
        nprocs = fboxlib.mpi_size()
        if (self.num_dim > 1) and (nprocs % 2**self.num_dim == 0):
            # divide domain into cubes
            nproc_per_dim = nprocs / 2**self.num_dim + 1
            for d in range(self.num_dim):
                max_sizes[d] = int(math.ceil(float(hi[d] - lo[d] + 1) / nproc_per_dim))
        else:
            # divide domain into slabs
            max_sizes[0] = int(math.ceil(float(hi[0] - lo[0] + 1) / nprocs))
            for d in range(1, self.num_dim):
                max_sizes[d] = hi[d] - lo[d] + 1

        ba = fboxlib.boxarray(boxes=[[lo, hi]])
        ba.maxsize(max_sizes)
        la = fboxlib.layout(ba)
        return la


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
