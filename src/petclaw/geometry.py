#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing petclaw.geometry.
"""

import pyclaw.geometry

class Patch(pyclaw.geometry.Patch):
    def __init__(self,dimensions):
        
        super(Patch,self).__init__(dimensions)

        self._da = self._create_DA()
        ranges = self._da.getRanges()
        grid_dimensions = []
        for i,nrange in enumerate(ranges):
            lower = self.lower[i] + nrange[0]*self.delta[i]
            upper = self.lower[i] + nrange[1]*self.delta[i]
            num_cells   = nrange[1]-nrange[0]
            grid_dimensions.append(pyclaw.geometry.Dimension(lower,upper,
                                        num_cells,name=dimensions[i].name))

        self.grid = pyclaw.geometry.Grid(grid_dimensions)


    def _create_DA(self):
        r"""Returns a PETSc DA and associated global Vec.
        Note that no local vector is returned.
        """
        from petsc4py import PETSc

        if hasattr(PETSc.DA, 'PeriodicType'):
            if self.num_dim == 1:
                periodic_type = PETSc.DA.PeriodicType.X
            elif self.num_dim == 2:
                periodic_type = PETSc.DA.PeriodicType.XY
            elif self.num_dim == 3:
                periodic_type = PETSc.DA.PeriodicType.XYZ
            else:
                raise Exception("Invalid number of dimensions")

            DA = PETSc.DA().create(dim=self.num_dim,
                                          dof=1,
                                          sizes=self.num_cells,
                                          periodic_type = periodic_type,
                                          stencil_width=0,
                                          comm=PETSc.COMM_WORLD)
        else:
            DA = PETSc.DA().create(dim=self.num_dim,
                                          dof=1,
                                          sizes=self.num_cells,
                                          boundary_type = PETSc.DA.BoundaryType.PERIODIC,
                                          stencil_width=0,
                                          comm=PETSc.COMM_WORLD)

        return DA

# ============================================================================
#  PetClaw Domain object definition
# ============================================================================
class Domain(pyclaw.geometry.Domain):
    r"""
    A Domain is a list of Patches.
    
    Need to add functionality to accept a list of patches as input.
    """
    def __init__(self,geom):
        if not isinstance(geom,list):
            geom = [geom]
        if isinstance(geom[0],Patch):
            self.patches = geom
        elif isinstance(geom[0],pyclaw.geometry.Dimension):
            self.patches = [Patch(geom)]


