#!/usr/bin/env python
# encoding: utf-8
r"""
Module containg the PetClaw solvers

This module contains the pure and wrapped PetClaw solvers.  All 
PetClaw solvers inherit from the :class:`ClawSolver` superclass which in turn 
inherits from the :class:`~petclaw.evolve.solver.Solver` superclass.  As such, 
the only solver classes that should be directly used should be the 
dimensionally dependent ones such as :class:`PetClawSolver1D`.

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

import numpy as np

from pyclaw.evolve.clawpack import ClawSolver, ClawSolver1D, start_step, src
from pyclaw.evolve import limiters

from petsc4py import PETSc


# ============================================================================
#  Generic PetClaw solver class
# ============================================================================
class PetClawSolver(ClawSolver):
    r"""
    Generic PetClaw solver
    
    All PetClaw solvers inherit from this base class.
    
    .. attribute:: mthlim 
    
        Limiter to be used on each wave.  ``Default = [1]``
    
    .. attribute:: order
    
        Order of the solver, either 1 for first order or 2 for second order 
        corrections.  ``Default = 2``
    
    .. attribute:: src_split
    
        Whether to use a source splitting method, 0 for none, 1 for first 
        order Godunov splitting and 2 for second order Strang splitting.
        ``Default = 0``
        
    .. attribute:: fwave
    
        Whether to split the flux into waves, requires that the Riemann solver
        performs the splitting.  ``Default = False``
        
    .. attribute:: src
    
        Source term function.  Default is the stub function.
    
    .. attribute:: start_step
    
        Function called before each time step is taken.  Default is the stub
        function
        
    
    :Initialization:
    
    Input:
     - *data* - (:class:`~petclaw.data.Data`) Data object, the solver will look 
       for the named variables to instantiate itself.    
    Output:
     - (:class:`PetClawSolver`) - Initialized petclaw solver
    """
    
    # ========== Generic Init Routine ========================================
    def __init__(self, kernelsType, data=None):
        r"""
        See :class:`ClawSolver` for full documentation.
        """
        
        self.kernelsType=kernelsType
        
        # Call general initialization function
        super(PetClawSolver,self).__init__(data)
    
         
    # ========== Time stepping routines ======================================
    def step(self,solutions):
        r"""
        Evolve solutions one time step

        This routine encodes the generic order in a full time step in this
        order:
        
        1. The :meth:`start_step` function is called
        
        2. A half step on the source term :func:`src` if Strang splitting is 
           being used (:attr:`src_split` = 2)
        
        3. A step on the homogeneous problem :math:`q_t + f(q)_x = 0` is taken
        
        4. A second half step or a full step is taken on the source term
           :func:`src` depending on whether Strang splitting was used 
           (:attr:`src_split` = 2) or Godunov splitting 
           (:attr:`src_split` = 1)

        This routine is called from the method evolve_to_time defined in the
        pyclaw.evolve.solver.Solver superclass.

        :Input:
         - *solutions* - (:class:`~petclaw.solution.Solution`) Dictionary of 
           solutions to be evolved
         
        :Output: 
         - (bool) - True if full step succeeded, False otherwise
        """
        # Grid we will be working on
        grid = solutions['n'].grids[0]
        # Number of equations
        meqn = solutions['n'].meqn
        maux = grid.maux
          
        q = self.qbc(grid.lqVec,grid)
        aux=grid.aux

        capa = grid.capa
        d = grid.d
        mbc = grid.mbc
        aux_global = grid.aux_global
        local_n = q.shape[0]
        

        # Call b4step, petclaw should be subclassed if this is needed
        self.start_step(self,solutions)

        # Source term splitting, petclaw should be subclassed if this 
        # is needed
        if self.src_split == 2:
            self.src(self,solutions,solutions['n'].t, self.dt/2.0)
    
        # Take a step on the homogeneous problem

        
        if(self.kernelsType == 'F'):
            from step1 import step1
            
            
            dt = self.dt
            dx = d[0]
            dtdx = np.zeros( (local_n) ) + dt/dx
            
            maux = grid.maux
            maxmx = local_n -mbc*2
            mx = maxmx
            
            if(aux == None):
                aux = np.empty( (local_n , maux) )
        
            method =np.ones(7, dtype=int) # hardcoded 7
            method[0] = self.dt_variable  # fixed or adjustable timestep
            method[1] = self.order  # order of the method
            method[2] = 0  # hardcoded 0, case of 2d or 3d
            method[3] = 0  # hardcoded 0 design issue: contorller.verbosity
            method[4] = self.src_split  # src term
            if (capa == None):
                method[5] = 0  #capa
            else:
                method[5] = 1  #capa. amal: mcapa no longer points to the capa componenets of the aux array as in fortran. capa now is a separate arry.
            method[6] = grid.maux  # aux
        
            mthlim = self.mthlim
        
            cfl = self.cfl
            f =  np.zeros( (local_n , meqn) )
            mwaves = meqn # amal: need to be modified

            wave = np.empty( (local_n,meqn,mwaves) )
            s = np.empty( (local_n,mwaves) )
            amdq = np.zeros( (local_n, meqn) )
            apdq = np.zeros( (local_n, meqn) )
        
            q,self.cfl = step1(maxmx,mbc,mx,q,aux,dx,dt,method,mthlim,f,wave,s,amdq,apdq,dtdx, -1)


        elif(self.kernelsType == 'P'):
            
            q = self.homogeneous_step( q, aux, capa, d, meqn,maux, mbc, aux_global)
        
        grid.q=q[mbc:-mbc,:]
        
        # Check here if we violated the CFL condition, if we did, return 
        # immediately to evolve_to_time and let it deal with picking a new
        # dt
        if self.cfl >= self.cfl_max:
            return False

        # Strang splitting
        if self.src_split == 2:
            self.src(self,solutions,solutions['n'].t + self.dt/2.0, self.dt/2.0)

        # Godunov Splitting
        if self.src_split == 1:
            self.src(self,solutions,solutions['n'].t,self.dt)
            
        return True
            
    def homogeneous_step(self,q, aux, capa, d, meqn, mbc, aux_global):
        r"""
        Take one homogeneous step on the solutions
        
        This is a dummy routine and must be overridden.
        """
        raise Exception("Dummy routine, please override!")
            

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


# ============================================================================
#  ClawPack 1d Solver Class
# ============================================================================
class PetClawSolver1D(PetClawSolver,ClawSolver1D):
    r"""
    PetClaw evolution routine in 1D
    
    This class represents the 1d clawpack solver on a single grid.  Note that 
    there are routines here for interfacing with the fortran time stepping 
    routines and the python time stepping routines.  The ones used are 
    dependent on the argument given to the initialization of the solver 
    (defaults to python).
    
    .. attribute:: rp
    
        Riemann solver function.
        
    :Initialization:
    
    Input:
     - *data* - (:class:`~petclaw.data.Data`) An instance of a Data object whose
       parameters can be used to initialize this solver
    Output:
     - (:class:`ClawSolver1D`) - Initialized 1d clawpack solver
        
    Need to check if we can simplify using multiple inheritance.

    :Authors:
        Amal Alghamdi
        David Ketcheson
    """

    def __init__(self,kernelsType,data=None):
        r"""
        Create 1d PetClaw solver
        
        See :class:`PetClawSolver1D` for more info.
        """   
        
        super(PetClawSolver1D,self).__init__(kernelsType,data)

    # ========== Python Homogeneous Step =====================================
    def homogeneous_step(self,q, aux, capa, d, meqn, maux, mbc, aux_global):
        r"""
        Take one time step on the homogeneous hyperbolic system

        Takes one time step of size dt on the hyperbolic system defined in the
        appropriate Riemann solver rp.
        """
    
        # Limiter to use in the pth family
        limiter = np.array(self.mthlim,ndmin=1)
         
        local_n = q.shape[0]
        # Flux vector
        f = np.empty( (local_n, meqn) )
    
        dtdx = np.zeros( (local_n) )

        # Find local value for dt/dx
        if capa is not None:
            dtdx = self.dt / (d[0] * capa)
        else:
            dtdx += self.dt/d[0]

        # Solve Riemann problem at each interface
        q_l=q[:-1,:]
        q_r=q[1:,:]
        if aux is not None:
            aux_l=aux[:-1,:]
            aux_r=aux[1:,:]
        else:
            aux_l = None
            aux_r = None
        wave,s,amdq,apdq = self.rp(q_l,q_r,aux_l,aux_r,aux_global)
        
        
        # Update loop limits, these are the limits for the Riemann solver
        # locations, which then update a grid cell value
        # We include the Riemann problem just outside of the grid so we can
        # do proper limiting at the grid edges
        #        LL    |                               |     UL
        #  |  LL |     |     |     |  ...  |     |     |  UL  |     |
        #              |                               |

       
        LL = mbc - 1
        UL = local_n - mbc + 1

        # Update q for Godunov update
        for m in xrange(meqn):
            q[LL:UL,m] -= dtdx[LL:UL]*apdq[LL-1:UL-1,m]
            q[LL-1:UL-1,m] -= dtdx[LL-1:UL-1]*amdq[LL-1:UL-1,m]
    
        # Compute maximum wave speed
        self.cfl = 0.0
        for mw in xrange(wave.shape[2]):
            smax1 = max(dtdx[LL:UL]*s[LL-1:UL-1,mw])
            smax2 = max(-dtdx[LL-1:UL-1]*s[LL-1:UL-1,mw])
            self.cfl = max(self.cfl,smax1,smax2)

        # comunicate max cfl
        if self.dt_variable:
            try:
                from mpi4py import MPI  
                comm = MPI.COMM_WORLD #Amal:should be consistent with petsc commworld
                size = comm.Get_size()
                rank = comm.Get_rank()
                max_cfl = 0
                max_cfl =comm.reduce( sendobj=self.cfl, op=MPI.MAX,  root=0)
                #max_cfl =comm.Reduce( self.cfl, max_cfl, op=MPI.MAX,  root=0)
                self.cfl = comm.bcast(max_cfl, root=0)
            except:
                raise Exception("Unable to communicate cfl")

        # If we are doing slope limiting we have more work to do
        if self.order == 2:
            # Initialize flux corrections
            f = np.zeros( (local_n + 2*mbc, meqn) )
        
            # Apply Limiters to waves
            if (limiter > 0).any():
                wave = limiters.limit(meqn,wave,s,limiter,dtdx)

            # Compute correction fluxes for second order q_{xx} terms
            dtdxave = 0.5 * (dtdx[LL-1:UL-1] + dtdx[LL:UL])
            if self.fwave:
                for mw in xrange(wave.shape[2]):
                    sabs = np.abs(s[LL-1:UL-1,mw])
                    om = 1.0 - sabs*dtdxave[:UL-LL]
                    ssign = np.sign(s[LL-1:UL-1,mw])
                    for m in xrange(meqn):
                        f[LL:UL,m] += 0.5 * ssign * om * wave[LL-1:UL-1,m,mw]
            else:
                for mw in xrange(wave.shape[2]):
                    sabs = np.abs(s[LL-1:UL-1,mw])
                    om = 1.0 - sabs*dtdxave[:UL-LL]
                    
                    for m in xrange(meqn):
                        f[LL:UL,m] += 0.5 * sabs * om * wave[LL-1:UL-1,m,mw]

            # Update q by differencing correction fluxes
            for m in xrange(meqn):
                q[LL:UL-1,m] -= dtdx[LL:UL-1] * (f[LL+1:UL,m] - f[LL:UL-1,m]) 
            
        return q
