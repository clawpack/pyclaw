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

from pyclaw.evolve.clawpack import start_step, src
from clawpack import PetClawSolver
from pyclaw.evolve.sharpclaw import SharpClawSolver1D
from petsc4py import PETSc

class RKStageDA(object):
    """
    A single Runge-Kutta stage.
    """
    def __init__(self,grid):
        self.t = grid.t
        periodic = False
        for dimension in grid.dimensions:
            if dimension.mthbc_lower == 2 or dimension.mthbc_upper == 2:
                periodic = True
                break
                
        if hasattr(PETSc.DA,'PeriodicType'):  # PETSc-3.1
            if grid.ndim == 1:
                periodic_type = PETSc.DA.PeriodicType.X
            elif grid.ndim == 2:
                periodic_type = PETSc.DA.PeriodicType.XY
            elif grid.ndim == 3:
                periodic_type = PETSc.DA.PeriodicType.XYZ
            else:
                raise Exception("Invalid number of dimensions")

            self.q_da = PETSc.DA().create(dim=grid.ndim,
                                          dof=grid.meqn,
                                          sizes=grid.n,
                                          periodic_type = periodic_type,
                                          stencil_width=grid.mbc,
                                          comm=PETSc.COMM_WORLD)
        else:
            self.q_da = PETSc.DA().create(dim=grid.ndim,
                                          dof=grid.meqn,
                                          sizes=grid.n,
                                          boundary_type = PETSc.DA.BoundaryType.PERIODIC,
                                          stencil_width=grid.mbc,
                                          comm=PETSc.COMM_WORLD)

        if grid.ndim == 1:
            self.q_da.setUniformCoordinates(xmin=grid.x.lower,xmax=grid.x.upper)
        elif grid.ndim == 2:
            self.q_da.setUniformCoordinates(xmin=grid.x.lower,xmax=grid.x.upper,ymin=grid.y.lower,ymax=grid.y.upper)
        elif grid.ndim == 3:
            self.q_da.setUniformCoordinates(xmin=grid.x.lower,xmax=grid.x.upper,ymin=grid.y.lower,ymax=grid.y.upper,zmin=grid.z.lower,zmax=grid.z.upper)
        self.gqVec = self.q_da.createGlobalVector()
        self.lqVec = self.q_da.createLocalVector()
        self.mbc  = grid.mbc
        self.meqn = grid.meqn
        self.ndim = grid.ndim

    def local_n():
        def fget(self):
            shape = [i[1]-i[0] for i in self.q_da.getRanges()]
            return shape
        return locals()
    def q():
        def fget(self):
            q_dim = self.local_n
            q_dim.insert(0,self.meqn)
            q=self.gqVec.getArray().reshape(q_dim, order = 'F')
            return q
        def fset(self,q):
            self.gqVec.setArray(q.reshape([-1], order = 'F'))
        return locals()

    def ghosted_q():
        def fget(self):
            self.q_da.globalToLocal(self.gqVec, self.lqVec)
            q_dim = [self.local_n[i] + 2*self.mbc for i in xrange(self.ndim)]
            q_dim.insert(0,self.meqn)
            ghosted_q=self.lqVec.getArray().reshape(q_dim, order = 'F')
            return ghosted_q
        def fset(self,ghosted_q):
            self.lqVec.setArray(ghosted_q.reshape([-1], order = 'F'))
        return locals()

    local_n     = property(**local_n())
    q           = property(**q())
    ghosted_q   = property(**ghosted_q())
 
class SharpPetClawSolver1D(SharpClawSolver1D,PetClawSolver):
    """SharpClaw evolution routine in 1D
    
    This class represents the 1d SharpClaw solver.  Note that there are 
    routines here for interfacing with the fortran time stepping routines and
    the python time stepping routines.  The ones used are dependent on the 
    argument given to the initialization of the solver (defaults to fortran).
    
    """
    
    def qbc(self,grid,q,t):
        return PetClawSolver.qbc(self,grid,q,t)

    def setup(self,solutions):

        if self.time_integrator == 'Euler': nregisters=1
        elif self.time_integrator == 'SSP33': nregisters=2
 
        grid = solutions['n'].grids[0]
        self.rk_stages = []
        for i in range(nregisters-1):
            self.rk_stages.append(RKStageDA(grid))



    def dqdt(self,grid,rk_stage):
        """
        Evaluate dq/dt
        """

        q = self.qbc(grid,rk_stage.q,rk_stage.t)

        self.dt = 1
        deltaq = self.dq_homogeneous(grid,q,rk_stage.t)

        # Godunov Splitting -- really the source term should be called inside rkstep
        if self.src_term == 1:
            deltaq+=self.src(grid,q,rk_stage.t)

        return deltaq.flatten('f')
            

