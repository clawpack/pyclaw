#!/usr/bin/env python
# encoding: utf-8
r"""
This module contains the most abstract parallel solver class, PetSolver.
All parallel solvers inherit from this class.

:Authors:
    Amal Alghamdi
    David Ketcheson
"""
# ============================================================================
#      Copyright (C) 2010 David I. Ketcheson <david.ketcheson@kaust.edu.sa>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

from pyclaw.solver import Solver

# ============================================================================
#  Generic PetClaw solver class
# ============================================================================
class PetSolver(Solver):
    r"""
    Generic PetClaw solver
    
    All PetClaw solvers inherit from this base class.
    See superclass Solver for documentation of attributes.
    
    :Initialization:
    
    Input:
     - *data* - (:class:`~petclaw.data.Data`) Data object, the solver will look 
       for the named variables to instantiate itself.    
    Output:
     - (:class:`PetClawSolver`) - Initialized petclaw solver
    """
    
    # ========== Boundary Conditions ==================================
    def append_ghost_cells(self,grid,state,q):
        """
        Returns q with ghost cells attached.  For PetSolver,
        this means returning the local vector.  
        The arguments 'q' and 'grid' are
        not used here; they are passed only in order to have a common
        interface for the petclaw and pyclaw versions of this function.
        """
        state.q_da.globalToLocal(state.gqVec, state.lqVec)
        q_dim = [state.local_n[i] + 2*self.mbc for i in xrange(state.ndim)]
        q_dim.insert(0,state.meqn)
        ghosted_q=state.lqVec.getArray().reshape(q_dim, order = 'F')
        return ghosted_q
 
    def communicateCFL(self):
        from petsc4py import PETSc

        if self.dt_variable:
            cflVec = PETSc.Vec().createWithArray([self.cfl])
            self.cfl = cflVec.max()[1]
