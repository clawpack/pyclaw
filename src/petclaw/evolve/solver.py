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

from pyclaw.evolve.clawpack import ClawSolver, ClawSolver1D, ClawSolver2D, start_step, src
from pyclaw.evolve import limiters

from petsc4py import PETSc

#This should be modified so we don't depend on mpi4py:
try:
  from mpi4py import MPI
except:
  raise Exception("Unable to communicate cfl")

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
    def __init__(self, kernelsType='F', data=None):
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

        # Call b4step, pyclaw should be subclassed if this is needed
        self.start_step(self,solutions)

        # Source term splitting, pyclaw should be subclassed if this 
        # is needed
        if self.src_split == 2:
            self.src(self,solutions,solutions['n'].t, self.dt/2.0)
    
        # Take a step on the homogeneous problem
        self.homogeneous_step(solutions)

        # comunicate max cfl
        if self.dt_variable:
          comm = MPI.COMM_WORLD #Amal:should be consistent with petsc commworld
          max_cfl = np.array([0.])
          cfl1 = np.array([self.cfl])
          comm.Allreduce(cfl1, max_cfl, MPI.MAX)
          self.cfl = max_cfl[0]
          

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
      
            
    def homogeneous_step(self,solutions):
        r"""
        Take one homogeneous step on the solutions
        
        This is a dummy routine and must be overridden.
        """
        raise Exception("Dummy routine, please override!")
            

    # ========== Boundary Conditions ==================================
    def qbc(self,grid):
        """
        Accepts an array qbc that includes ghost cells.
        Returns an array with the ghost cells filled.
        It would be nice to do the ghost cell array fetch in here, but
        we need to think about how to associate q_da and gqVec, lqVec.

        For now, grid and dim are passed in for backward compatibility.
        We should think about what makes the most sense.
        """
        
        qbc=  grid.ghosted_q 
        for i in xrange(len(grid._dimensions)):
            dim = getattr(grid,grid._dimensions[i])
            #If a user defined boundary condition is being used, send it on,
            #otherwise roll the axis to front position and operate on it
            if dim.mthbc_lower == 0:
                self.qbc_lower(qbc,grid,dim)
            else:
                self.qbc_lower(np.rollaxis(qbc,i+1,1),grid,dim)
            if dim.mthbc_upper == 0:
                self.qbc_upper(qbc,grid,dim)
            else:
                self.qbc_upper(np.rollaxis(qbc,i+1,1),grid,dim)
        return qbc

    def qbc_lower(self,qbc,grid,dim):
        r"""
        
        """
        # User defined functions
        if dim.mthbc_lower == 0:
            self.user_bc_lower(grid,dim,qbc)
        # Zero-order extrapolation
        elif dim.mthbc_lower == 1:
            if dim.nstart == 0:
                for i in xrange(grid.mbc):
                    qbc[:,i,...] = qbc[:,grid.mbc,...]
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
            if dim.nend == dim.n :
                for i in xrange(grid.mbc):
                    qbc[:,-i-1,...] = qbc[:,-grid.mbc-1,...] 
 	    
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
    def homogeneous_step(self,solutions):
        r"""
        Take one time step on the homogeneous hyperbolic system

        Takes one time step of size dt on the hyperbolic system defined in the
        appropriate Riemann solver rp.
        """
        # Grid we will be working on
        grid = solutions['n'].grids[0]
        # Number of equations
        meqn,maux,mwaves,mbc,aux = grid.meqn,grid.maux,self.mwaves,grid.mbc,grid.aux
          
        q = self.qbc(grid)

        local_n = q.shape[1]

        if(self.kernelsType == 'F'):
            from step1 import step1
            
            dx,dt = grid.d[0],self.dt
            dtdx = np.zeros( (local_n) ) + dt/dx
            maxmx = local_n -mbc*2
            mx = maxmx
            
            if(aux == None): aux = np.empty( (maux,local_n) )
        
            method =np.ones(7, dtype=int) # hardcoded 7
            method[0] = self.dt_variable  # fixed or adjustable timestep
            method[1] = self.order  # order of the method
            method[2] = 0  # hardcoded 0, case of 2d or 3d
            method[3] = 0  # hardcoded 0 design issue: contorller.verbosity
            method[4] = self.src_split  # src term
            if (grid.capa == None):
                method[5] = 0  #capa
            else:
                method[5] = 1  #capa. amal: mcapa no longer points to the capa componenets of the aux array as in fortran. capa now is a separate arry.
            method[6] = maux  # aux
        
            f    = np.empty( (meqn,local_n) )
            wave = np.empty( (meqn,mwaves,local_n) )
            s    = np.empty( (mwaves,local_n) )
            amdq = np.empty( (meqn,local_n) )
            apdq = np.empty( (meqn,local_n) )
        
            q,self.cfl = step1(maxmx,mbc,mx,q,aux,dx,dt,method,self.mthlim,f,wave,s,amdq,apdq,dtdx)

        elif(self.kernelsType == 'P'):
           
            # Limiter to use in the pth family
            limiter = np.array(self.mthlim,ndmin=1)  
            # Q with appended boundary conditions
            q = self.qbc(grid)
            # Flux vector
            f = np.empty( (2*grid.mbc + grid.n[0], meqn) )
        
            dtdx = np.zeros( (2*grid.mbc+grid.n[0]) )

            # Find local value for dt/dx
            if grid.capa is not None:
                dtdx = self.dt / (grid.d[0] * grid.capa)
            else:
                dtdx += self.dt/grid.d[0]
        
            # Solve Riemann problem at each interface
            q_l=q[:-1,:]
            q_r=q[1:,:]
            if grid.aux is not None:
                aux_l=grid.aux[:-1,:]
                aux_r=grid.aux[1:,:]
            else:
                aux_l = None
                aux_r = None
            wave,s,amdq,apdq = self.rp(q_l,q_r,aux_l,aux_r,grid.aux_global)
            
            # Update loop limits, these are the limits for the Riemann solver
            # locations, which then update a grid cell value
            # We include the Riemann problem just outside of the grid so we can
            # do proper limiting at the grid edges
            #        LL    |                               |     UL
            #  |  LL |     |     |     |  ...  |     |     |  UL  |     |
            #              |                               |
            LL = grid.mbc - 1
            UL = grid.mbc + grid.n[0] + 1 

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

            # If we are doing slope limiting we have more work to do
            if self.order == 2:
                # Initialize flux corrections
                f = np.zeros( (grid.n[0] + 2*grid.mbc, meqn) )
            
                # Apply Limiters to waves
                if (limiter > 0).any():
                    wave = limiters.limit(grid.meqn,wave,s,limiter,dtdx)

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
                
        
        grid.q=q[:,mbc:-mbc]



# ============================================================================
#  PetClaw 2d Solver Class
# ============================================================================
class PetClawSolver2D(PetClawSolver,ClawSolver2D):
    r"""
    PetClaw evolution routine in 2D
    
    This class represents the 2d clawpack solver on a single grid.  Note that 
    only the fortran routines are supported for now in 2D.
    
    :Initialization:
    
    Input:
     - *data* - (:class:`~petclaw.data.Data`) An instance of a Data object whose
       parameters can be used to initialize this solver
    Output:
     - (:class:`PetClawSolver1D`) - Initialized 1d clawpack solver
        
    Need to check if we can simplify using multiple inheritance.

    :Authors:
        Amal Alghamdi
        David Ketcheson
    """

    def __init__(self,kernelsType='F',data=None):
        r"""
        Create 2D PetClaw solver.
        See :class:`PetClawSolver2D` for more info.
        """   
        
        super(PetClawSolver2D,self).__init__(kernelsType,data)

    # ========== Python Homogeneous Step =====================================
    def homogeneous_step(self,solutions):
        r"""
        Take one time step on the homogeneous hyperbolic system.
        Only the dimensionally split algorithm is supported for now.
        """
        
        # Grid we will be working on
        grid = solutions['n'].grids[0]
        # Number of equations
        meqn,maux,mbc,mwaves = grid.meqn,grid.maux,grid.mbc,self.mwaves

        if(self.kernelsType == 'F'):
            from dimsp2 import dimsp2
            maxmx,maxmy = grid.local_n[0],grid.local_n[1]
            maxm = max(maxmx, maxmy)
            mx,my = maxmx,maxmy
            aux = grid.aux
            
            #Old workaround
            #if(aux == None): aux=np.empty([0]*(grid.ndim+1))

            #New workaround
            #The following is an awful hack to work around an issue
            #with f2py.  It involves wastefully allocating a large array.
            #It could be avoided by using the other workaround below, but
            #that makes the build process different for every example
            #because we would need to generate and modify the .pyf file.
            if(aux == None): 
                maux=1
                aux=np.empty((maux,maxmx+2*mbc,maxmy+2*mbc))
                
            dx,dy,dt = grid.d[0],grid.d[1],self.dt

            method =np.ones(7, dtype=int)
            method[0] = self.dt_variable
            method[1] = self.order
            method[2] = -1  # hardcoded 0, case of 2d or 3d
            method[3] = 0  # hardcoded 0 design issue: controller.verbosity
            method[4] = self.src_split  # src term

            # mcapa no longer points to the capa components of the aux 
            # array as in fortran. capa now is a separate array.
            if (grid.capa == None): method[5] = 0
            else: method[5] = 1  
            method[6] = maux
            
            cflv = np.zeros(4)
            cflv[0:2] = [self.cfl_max,self.cfl_desired]
            #cflv[2] and cflv[3] are output values.

            if method[4] < 2: narray = 1
            else: narray = 2

            mwork = (maxm+2*mbc) * (5*meqn + mwaves + meqn*mwaves) \
                       + (narray-1) * (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn 
            work = np.empty((mwork))
            
            qold = self.qbc(grid)
            qnew = qold #(input/output)

            #Old Workaround for f2py bug (?)
            #f2py Doesn't like fortran arrays with first dimension zero,
            # so if maux=0 we just pass empty 1x1 arrays for aux1,aux2,aux3.
            #if maux==0:
            #    q, cfl = dimsp2(maxm,maxmx,maxmy,mbc,mx,my, \
            #              qold,qnew,aux,dx,dy,dt,method,self.mthlim,self.cfl,cflv, \
            #              np.empty((meqn,maxm+2*mbc)), \
            #              np.empty((meqn,maxm+2*mbc)), \
            #              np.empty((meqn,2,maxm+2*mbc)), \
            #              np.empty((meqn,maxm+2*mbc)), \
            #              np.empty((maxmx+2*mbc)), np.empty((maxmy+2*mbc)), \
            #              np.empty((1,1)), np.empty((1,1)), np.empty((1,1)), \
            #              work)
            #else:
            q, cfl = dimsp2(maxm,maxmx,maxmy,mbc,mx,my, \
                      qold,qnew,aux,dx,dy,dt,method,self.mthlim,self.cfl,cflv, \
                      np.empty((meqn,maxm+2*mbc)), \
                      np.empty((meqn,maxm+2*mbc)), \
                      np.empty((meqn,2,maxm+2*mbc)), \
                      np.empty((meqn,maxm+2*mbc)), \
                      np.empty((maxm+2*mbc)), np.empty((maxm+2*mbc)), \
                      np.empty((maux,maxm+2*mbc)), \
                      np.empty((maux,maxm+2*mbc)), \
                      np.empty((maux,maxm+2*mbc)), \
                      work)

            self.cfl = cfl
            grid.q=q[:,mbc:grid.local_n[0]+mbc,mbc:grid.local_n[1]+mbc]

        elif(self.kernelsType == 'P'):
            raise NotImplementedError("No python implementation for homogeneous_step in case of 2D.")

        
    
