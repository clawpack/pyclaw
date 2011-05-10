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
from clawpack import PetSolver
from pyclaw.evolve.sharpclaw import SharpClawSolver1D
from petsc4py import PETSc

class RKStageState(object):
    """
    A single Runge-Kutta stage.
    """
    def __init__(self,grid):
        #Probably we should pass a solver instead of a grid to this function.
        self.t = grid.t
        self.q_da=grid.q_da
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
 
class PetSharpClawSolver1D(PetSolver,SharpClawSolver1D):
    """SharpClaw evolution routine in 1D
    
    This class represents the 1d SharpClaw solver.  Note that there are 
    routines here for interfacing with the fortran time stepping routines and
    the python time stepping routines.  The ones used are dependent on the 
    argument given to the initialization of the solver (defaults to fortran).
    
    """
    
    def setup(self,solutions):
        """
        Allocate RK stage arrays.
        """

        if self.time_integrator == 'Euler': nregisters=1
        elif self.time_integrator == 'SSP33': nregisters=2
 
        grid = solutions['n'].grids[0]
        self.rk_stages = []
        for i in range(nregisters-1):
            self.rk_stages.append(RKStageState(grid))


    def dqdt(self,grid,rk_stage):
        """
        Evaluate dq/dt
        """

        q = self.qbc(grid,rk_stage)

        self.dt = 1
        deltaq = self.dq_homogeneous(grid,q,rk_stage.t)

        # Godunov Splitting -- really the source term should be called inside rkstep
        if self.src_term == 1:
            deltaq+=self.src(grid,q,rk_stage.t)

        return deltaq.flatten('f')
            

