#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing SharpClaw solvers for PyClaw/PetClaw

#  File:        sharpclaw.py
#  Created:     2010-03-20
#  Author:      David Ketcheson
"""
# ============================================================================
#      Copyright (C) 2010 David I. Ketcheson <david.ketcheson@kaust.edu.sa>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import numpy as np

# Solver superclass
from pyclaw.evolve.clawpack import start_step, src
from solver import PetClawSolver
from pyclaw.evolve.sharpclaw import SharpClawSolver1D
from petsc4py import PETSc

class RKStage(object):
    """
    A single Runge-Kutta stage.
    """
    def __init__(self,grid):
        periodic = False
        for dimension in grid.dimensions:
            if dimension.mthbc_lower == 2 or dimension.mthbc_upper == 2:
                periodic = True
                break
                
        if hasattr(PETSc.DA,'PeriodicType'):
            if grid.ndim == 1:
                periodic_type = PETSc.DA.PeriodicType.X
            elif grid.ndim == 2:
                periodic_type = PETSc.DA.PeriodicType.XY
            elif grid.ndim == 3:
                periodic_type = PETSc.DA.PeriodicType.XYZ
            else:
                raise Exception("Invalid number of dimensions")

        self.da = PETSc.DA().create(dim=grid.ndim,
                                    dof=grid.meqn,
                                    sizes=grid.n, 
                                    periodic_type = periodic_type,
                                    stencil_width=grid.mbc,
                                    comm=PETSc.COMM_WORLD)
        self.gVec = self.da.createGlobalVector()
        self.lVec = self.da.createLocalVector()


class SharpPetClawSolver1D(SharpClawSolver1D,PetClawSolver):
    """SharpClaw evolution routine in 1D
    
    This class represents the 1d SharpClaw solver.  Note that there are 
    routines here for interfacing with the fortran time stepping routines and
    the python time stepping routines.  The ones used are dependent on the 
    argument given to the initialization of the solver (defaults to fortran).
    
    """
    
    def qbc(self,grid,q,t):
        return PetClawSolver.qbc(self,grid,q,t)

