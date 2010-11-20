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

class qrk(object):
    """
    A single Runge-Kutta stage
    """
    def __init__(self,grid):
        periodic = False
        for dimension in grid.dimensions:
            if dimension.mthbc_lower == 2 or dimension.mthbc_upper == 2:
                periodic = True
                break
                
        if grid.ndim == 1:
            if periodic: periodic_type = PETSc.DA.PeriodicType.X
            else: periodic_type = PETSc.DA.PeriodicType.GHOSTED_XYZ
        elif grid.ndim == 2:
            if periodic: periodic_type = PETSc.DA.PeriodicType.XY
            else: periodic_type = PETSc.DA.PeriodicType.GHOSTED_XYZ
        elif grid.ndim == 3:
            if periodic: periodic_type = PETSc.DA.PeriodicType.XYZ
            else: periodic_type = PETSc.DA.PeriodicType.GHOSTED_XYZ
        else:
            raise Exception("Invalid number of dimensions")

        self.da = PETSc.DA().create(dim=grid.ndim,
                                    dof=grid.meqn,
                                    sizes=grid.n, 
                                    periodic_type = periodic_type,
                                    #stencil_type=grid.STENCIL,
                                    stencil_width=grid.mbc,
                                    comm=PETSc.COMM_WORLD)
        self.gVec = self.da.createGlobalVector()
        self.lVec = self.da.createLocalVector()

        if self.time_integrator=='SSP33':
            self.qrk = Solution(grid)


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
        self.src_term=0
        
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
        from pyclaw.solution import Solution
        # Grid we will be working on
        grid = solutions['n'].grids[0]

        if self.time_integrator=='Euler':
            deltaq=self.dq(solutions)
            grid.q+=deltaq
        elif self.time_integrator=='SSP33':
            #solutions['rk1']=Solution(grid.__deepcopy__())
            qold = grid.q.copy()
            told = solutions['n'].t
            deltaq=self.dq(solutions)
            grid.q+=deltaq
            solutions['n'].t=told+self.dt
            deltaq=self.dq(solutions)
            grid.q = 0.75*qold + 0.25*(grid.q+deltaq)
            solutions['n'].t=told+0.5*self.dt
            deltaq=self.dq(solutions)
            grid.q = 1./3.*qold + 2./3.*(grid.q+deltaq)

        else:
            raise Exception('Unrecognized time integrator')

        
    def dq(self,solutions):
        """
        Take one Runge-Kutta time step.
        Right now this is just Euler for debugging.
        Each RK stage should be a Solution object.
        """

        self.start_step(self,solutions)

        # Take a step on the homogeneous problem

        # Grid we will be working on
        grid = solutions['n'].grids[0]
        # Number of equations
        meqn = grid.meqn
        # Q with appended boundary conditions

        maux = grid.maux
        aux=grid.aux

        capa = grid.capa
        d = grid.d
        mbc = grid.mbc
        aux_global = grid.aux_global
 
        q = self.qbc(grid.lqVec,grid)
        local_n = q.shape[0]

        t=solutions['n'].t
        deltaq=self.dqhyp1(q,t,aux,capa,d,meqn,maux,mbc,aux_global)

        # Check here if we violated the CFL condition, if we did, return 
        # immediately to evolve_to_time and let it deal with picking a new
        # dt
        if self.cfl >= self.cfl_max:
            return False

        # Godunov Splitting -- really the source term should be called inside rkstep
        if self.src_term == 1:
            deltaq+=self.src(solutions['n'],solutions['n'].t,self.dt)


        return deltaq


    def dqhyp1(self,q, t, aux, capa, d, meqn, maux, mbc, aux_global):
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

        dq = np.empty(q.shape)
        mwaves = meqn # amal: need to be modified

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
                #HACK
                ql=q; qr=q;

        # Solve Riemann problem at each interface
        q_l=qr[:-1,:]
        q_r=ql[1: ,:]
        wave,s,amdq,apdq = self.rp(q_l,q_r,aux_l,aux_r,aux_global)

        # Loop limits for local potion of grid
        LL = mbc - 1
        UL = local_n - mbc + 1

        # Compute maximum wave speed
        self.cfl = 0.0
        for mw in xrange(mwaves):
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

    # ========== Boundary Conditions ==================================
    def qbc(self,lqVec,grid):
        """
        Accepts an array qbc that includes ghost cells.
        Returns an array with the ghost cells filled.
        It would be nice to do the ghost cell array fetch in here, but
        we need to think about how to associate q_da and gqVec, lqVec.

        For now, grid and dim are passed in for backward compatibility.
        We should think about what makes the most sense.
        """
        #THIS ONLY WORKS IN 1D:
        qbc=lqVec.getArray().reshape([-1,grid.meqn])
        for i in xrange(len(grid._dimensions)):
            dim = getattr(grid,grid._dimensions[i])
            #If a user defined boundary condition is being used, send it on,
            #otherwise roll the axis to front position and operate on it
            if dim.mthbc_lower == 0:
                self.qbc_lower(qbc,grid,dim)
            else:
                self.qbc_lower(np.rollaxis(qbc,i),grid,dim)
            if dim.mthbc_upper == 0:
                self.qbc_upper(qbc,grid,dim)
            else:
                self.qbc_upper(np.rollaxis(qbc,i),grid,dim)
        return qbc

    def qbc_lower(self,qbc,grid,dim):
        r"""
        
        """
        # User defined functions
        if dim.mthbc_lower == 0:
            self.user_bc_lower(grid,dim,qbc)
        # Zero-order extrapolation
        elif dim.mthbc_lower == 1:
            ##specify rank
            rank = PETSc.Comm.getRank(PETSc.COMM_WORLD) # Amal: hardcoded communicator
            print rank
            if rank == 0:
                qbc[:grid.mbc,...] = qbc[grid.mbc,...]
        # Periodic
        elif dim.mthbc_lower == 2:
            pass # Amal: this is implemented automatically by petsc4py
            
        # Solid wall bc
        elif dim.mthbc_lower == 3:
            raise NotImplementedError("Solid wall upper boundary condition not implemented.")
        else:
            raise NotImplementedError("Boundary condition %s not implemented" % x.mthbc_lower)

    def qbc_upper(self,qbc,grid,dim):
        r"""
        
        """
        # User defined functions
        if dim.mthbc_upper == 0:
            self.user_bc_upper(grid,dim,qbc)
        # Zero-order extrapolation
        elif dim.mthbc_upper == 1:
            rank = PETSc.Comm.getRank(PETSc.COMM_WORLD) # Amal: hardcoded communicator
            size = PETSc.Comm.getSize(PETSc.COMM_WORLD)
            
            if rank == size-1:
                local_n = grid.q.shape[0]
                list_from =[local_n - grid.mbc -1]*grid.mbc    
                list_to = range(local_n - grid.mbc, local_n )
                grid.q[list_to,:]=grid.q[list_from,:]
 	    
        elif dim.mthbc_upper == 2:
            # Periodic
            pass # Amal: this is implemented automatically by petsc4py

        # Solid wall bc
        elif dim.mthbc_upper == 3:
            raise NotImplementedError("Solid wall upper boundary condition not implemented.")

        else:
            raise NotImplementedError("Boundary condition %s not implemented" % x.mthbc_lower)


