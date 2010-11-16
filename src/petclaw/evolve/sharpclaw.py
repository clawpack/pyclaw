#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing SharpClaw solvers for PyClaw/PetClaw

#  File:        sharpclaw.py
#  Created:     2009-04-07
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
from pyclaw.evolve.clawpack import ClawSolver, ClawSolver1D, start_step, src
from petsc4py import PETSc

# Limiters
import recon

class SharpClawSolver1D(ClawSolver1D):
    """SharpClaw evolution routine in 1D
    
    This class represents the 1d SharpClaw solver.  Note that there are 
    routines here for interfacing with the fortran time stepping routines and
    the python time stepping routines.  The ones used are dependent on the 
    argument given to the initialization of the solver (defaults to fortran).
    
    """
    
    attributes = ['dt_variable','max_steps','cfl_max','cfl_desired',
            'dt_max','dt_initial','mwaves','lim_type','mthlim',
            'time_integrator','char_decomp','t0','src_term']
    
    # ========================================================================
    #   Initialization routines
    # ========================================================================
    def __init__(self, kernelsType, data=None):
        r"""
        Here we just set the flag for using Python or Fortran kernels.
        """
        
        self.kernelsType=kernelsType
        
        # Call general initialization function
        super(SharpClawSolver1D,self).__init__(data)
 

    # ========== Time stepping routines ======================================
    def step(self,solutions):
        """Evolve q over one time step.

        Takes the appropriate times steps to reach tend.  This is done based 
        on the current solution attributes set in the object.

        Arguments:
          tend - Target time step to reach
        """


        self.start_step(self,solutions)

        # Take a step on the homogeneous problem
        self.rkstep(solutions)

        # Check here if we violated the CFL condition, if we did, return 
        # immediately to evolve_to_time and let it deal with picking a new
        # dt
        if self.cfl >= self.cfl_max:
            return False

        # Godunov Splitting -- really the source term should be called inside rkstep
        if self.src_term == 1:
            self.src(solutions,solutions['n'].t,self.dt)

        
    def rkstep(self,solutions):
        """
        Take one Runge-Kutta time step.
        Right now this is just Euler for debugging.
        Each RK stage should be a Solution object.
        """
        # Grid we will be working on
        grid = solution.grids[0]
        # Number of equations
        meqn = grid.meqn
        # Q with appended boundary conditions
        q = grid.qbc()

        maux = grid.maux
        aux=grid.aux

        capa = grid.capa
        d = grid.d
        mbc = grid.mbc
        aux_global = grid.aux_global
        local_n = q.shape[0]
 
        qold=solutions['n'].q
        told=solutions['n'].t
        dq=self.dq1(solutions['n'],qold,told)
        solutions['n'].q = qold+dq
        #q2 = qold+dq
        #t2 = told+self.dt
        #dq=self.dq1(solutions['n'],q2,t2)
        #solutions['n'].q = 0.5*(qold+q2+dq)


    def dq1(self,q, aux, capa, d, meqn, maux, mbc, aux_global):
        """Compute dq/dt * (delta t) for the homogeneous hyperbolic system

        Note that the capa array, if present, should be located in the aux
        variable.

        Indexing works like this:  here mbc=2 as an example
         0     1     2     3     4     mx+mbc-2     mx+mbc      mx+mbc+2
                     |                        mx+mbc-1 |  mx+mbc+1
         |     |     |     |     |   ...   |     |     |     |     |
            0     1  |  2     3            mx+mbc-2    |mx+mbc       
                                                  mx+mbc-1   mx+mbc+1

        The top indices represent the values that are located on the grid
        cell boundaries such as waves, s and other Riemann problem values, 
        the bottom for the cell centered values such as q.  In particular
        the ith grid cell boundary has the following related information:
                          i-1         i         i+1
                           |          |          |
                           |   i-1    |     i    |
                           |          |          |
        Again, grid cell boundary quantities are at the top, cell centered
        values are in the cell.

        """
    
        # Limiter to use in the pth family
        limiter = np.array(self.mthlim,ndmin=1)  
        lim_type=self.lim_type
        char_decomp=self.char_decomp
    
        local_n = q.shape[0]

        # Flux vector
        f = np.empty( (local_n, meqn) )
        dtdx = np.zeros( (local_n) )

        # Find local value for dt/dx
        if capa is not None:
            dtdx = self.dt / (d[0] * capa)
        else:
            dtdx += self.dt/d[0]

        dtdx = np.zeros( (2*grid.mbc+grid.mx) )
        dq = np.empty(q.shape)

        if aux is not None:
            aux_l=aux[:-1,:]
            aux_r=aux[1: ,:]
        else:
            aux_l = None
            aux_r = None
   
        #Reconstruct (wave reconstruction uses a Riemann solve)
        if lim_type==0: #Unlimited reconstruction; need to write
            pass
        elif lim_type==1: #TVD Reconstruction; need to fill in
            pass
        elif lim_type==2: #WENO Reconstruction
            if char_decomp==0: #No characteristic decomposition
                ql,qr=recon.weno5(q)
            elif char_decomp==1: #Wave-based reconstruction
                q_l=q[:-1,:]
                q_r=q[1: ,:]
                wave,s,amdq,apdq = self.rp(q_l,q_r,aux_l,aux_r,aux_global)
                ql,qr=recon.weno5_wave(q,wave,s)
            elif char_decomp==2: #Characteristic-wise reconstruction
                pass

        # Solve Riemann problem at each interface
        q_l=qr[:-1,:]
        q_r=ql[1: ,:]
        wave,s,amdq,apdq = self.rp(q_l,q_r,aux_l,aux_r,aux_global)

        # Loop limits for local potion of grid
        LL = mbc - 1
        UL = local_n - mbc + 1

        # Compute maximum wave speed
        self.cfl = 0.0
        for mw in xrange(self.mwaves):
            smax1 = max(dtdx[LL:UL]*s[LL-1:UL-1,mw])
            smax2 = max(-dtdx[LL-1:UL-1]*s[LL-1:UL-1,mw])
            self.cfl = max(self.cfl,smax1,smax2)

       
        #Find total fluctuation within each cell
        wave,s,amdq2,apdq2 = self.rp(ql,qr,aux,aux,aux_global)

        # Compute dq
        for m in xrange(meqn):
            dq[LL:UL,m] = -dtdx[LL:UL]*(amdq[LL:UL,m] + apdq[LL-1:UL-1,m] \
                            + apdq2[LL:UL,m] + amdq2[LL:UL,m])
    
        return dq[LL+1:UL-1]
