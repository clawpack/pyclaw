#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing forestclaw.geometry.
All we do here is override the pyclaw.Patch class to add a block_number and mpi_rank.
"""
from clawpack import pyclaw
from clawpack.pyclaw import geometry as pyclaw_geometry
from clawpack.pyclaw.geometry import Dimension
from clawpack.pyclaw.geometry import Domain


class Patch(pyclaw_geometry.Patch):
    """Patch class with specific ForestClaw attributes.
    """

    __doc__ += pyclaw.util.add_parent_doc(pyclaw_geometry.Patch)

    def __init__(self, dimensions):

        super(Patch, self).__init__(dimensions)

        self.block_number = 0
        r"""(int) - Block number of current patch, ``default = 0``"""
        self.mpi_rank = 0
        r"""(int) - MPI rank this patch belongs to, ``default = 0``"""
