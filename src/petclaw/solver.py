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

import pyclaw.solver

# ============================================================================
#  Generic PetClaw solver class
# ============================================================================
class PetSolver(pyclaw.solver.Solver):
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
    def append_ghost_cells(self,state):
        """
        Returns q with ghost cells attached.  For PetSolver,
        this means returning the local vector.  
        """
        state.q_da.globalToLocal(state.gqVec, state.lqVec)
        q_dim = [n + 2*self.mbc for n in state.grid.ng]
        q_dim.insert(0,state.meqn)
        ghosted_q=state.lqVec.getArray().reshape(q_dim, order = 'F')
        return ghosted_q
    
    def update_global_q(self,state,ghosted_q):
        """
        Update the value of q. for PetSolver, this involves setting ghosted_q
        as the local vector array then perform a local to global communication. 
        """
        state.lqVec.placeArray(ghosted_q)
        state.q_da.localToGlobal(state.lqVec,state.gqVec)
        state.lqVec.resetArray() # This call is required because placeArray is
                                 # intended to be temporarly placement

    def append_ghost_cells_to_aux(self,state):
        """
        Returns aux with ghost cells attached.  For PetSolver,
        this means returning the local vector.  

        The globalToLocal call here could be wasteful if the aux variables
        don't change in time.  We should add a flag for this and just
        do it once.
        """
        state.aux_da.globalToLocal(state.gauxVec, state.lauxVec)
        aux_dim = [n + 2*self.mbc for n in state.grid.ng]
        aux_dim.insert(0,state.maux)
        ghosted_aux=state.lauxVec.getArray().reshape(aux_dim, order = 'F')
        return ghosted_aux
 
    def communicateCFL(self):
        from petsc4py import PETSc

        if self.dt_variable:
            cflVec = PETSc.Vec().createWithArray([self.cfl])
            self.cfl = cflVec.max()[1]

    def allocate_rk_stages(self,solutions):
        r"""We could eliminate this function and just use
        the version in pyclaw.solver.Solver, if we were willing to
        check there whether the solver is a PetSolver.  But this
        would mean putting parallel-aware code in PyClaw, so for
        now we duplicate the function here.
        """
        from state import State

        if self.time_integrator   == 'Euler':  nregisters=1
        elif self.time_integrator == 'SSP33':  nregisters=2
        elif self.time_integrator == 'SSP104': nregisters=3
 
        state = solutions['n'].states[0]
        self.rk_stages = []
        for i in range(nregisters-1):
            self.rk_stages.append(State(state.grid))
            self.rk_stages[-1].meqn = state.meqn
            self.rk_stages[-1].maux = state.maux
            self.rk_stages[-1].set_stencil_width(self.mbc)
            self.rk_stages[-1].aux_global       = state.aux_global
            self.rk_stages[-1].t                = state.t
            if state.maux > 0:
                self.rk_stages[-1].aux              = state.aux



