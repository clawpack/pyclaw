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

from clawpack import PetSolver
from pyclaw.evolve.sharpclaw import SharpClawSolver1D, SharpClawSolver2D

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

    local_n     = property(**local_n())
    q           = property(**q())
 

class PetSharpClawSolver1D(PetSolver,SharpClawSolver1D):
    """
    
    1D parallel SharpClaw solver.

    Note that there are routines here for interfacing with the fortran time
    stepping routines and the python time stepping routines.  The ones used are
    dependent on the kernel_language argument given to the initialization of
    the solver (defaults to fortran).
    
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

        if self.kernel_language=='Fortran':
            from sharpclaw1 import clawparams, workspace, reconstruct
            self.set_fortran_parameters(grid,clawparams,workspace,reconstruct)

           
class PetSharpClawSolver2D(PetSolver,SharpClawSolver2D):
    """
    
    2D parallel SharpClaw solver.  
    
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

        if self.kernel_language=='Fortran':
            from sharpclaw2 import clawparams, workspace, reconstruct
            self.set_fortran_parameters(grid,clawparams,workspace,reconstruct)
