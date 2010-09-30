#!/usr/bin/env python
# encoding: utf-8
r"""
Module containg the classic Clawpack solvers

This module contains the pure and wrapped classic clawpack solvers.  All 
clawpack solvers inherit from the :class:`ClawSolver` superclass which in turn 
inherits from the :class:`~pyclaw.evolve.solver.Solver` superclass.  As such, 
the only solver classes that should be directly used should be the 
dimensionally dependent ones such as :class:`ClawSolver1D`.

:Authors:
    Kyle T. Mandli (2008-09-11) Initial version
"""
# ============================================================================
#      Copyright (C) 2008 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import numpy as np

from pyclaw.evolve.solver import Solver
from petsc4py import PETSc

import limiters

# ========================================================================
#  User-defined routines
# ========================================================================
def start_step(solver,solutions):
    r"""
    Dummy routine called before each step
    
    Replace this routine if you want to do something before each time step.
    """
    pass

def src(solver,solutions,t,dt):
    r"""
    Dummy routine called to calculate a source term
    
    Replace this routine if you want to include a source term.
    """
    pass

# ============================================================================
#  Generic Clawpack solver class
# ============================================================================
class ClawSolver(Solver):
    r"""
    Generic classic Clawpack solver
    
    All Clawpack solvers inherit from this base class.
    
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
     - *data* - (:class:`~pyclaw.data.Data`) Data object, the solver will look 
       for the named variables to instantiate itself.    
    Output:
     - (:class:`ClawSolver`) - Initialized clawpack solver
    
    :Version: 1.0 (2009-06-01)
    """
    
    # ========== Generic Init Routine ========================================
    def __init__(self,data=None):
        r"""
        See :class:`ClawSolver` for full documentation.
        """
        
        # Required attributes for this solver
        for attr in ['mthlim','order','src_split','fwave','src','start_step']:
            self._required_attrs.append(attr)
        
        # Default required attributes
        self._default_attr_values['mthlim'] = [1]
        self._default_attr_values['order'] = 2
        self._default_attr_values['src_split'] = 0
        self._default_attr_values['fwave'] = False
        self._default_attr_values['src'] = src
        self._default_attr_values['start_step'] = start_step

        # Call general initialization function
        super(ClawSolver,self).__init__(data)
    
    # ========== Setup Routine ===============================================
    def setup(self):
        r"""
        Called before any set of time steps.
        
        This routine will be called once before the solver is used via the
        :class:`~pyclaw.controller.Controller`.  In the case of 
        :class:`ClawSolver` we make sure that the :attr:`mthlim` is a list.
        """
    
        # Change mthlim to be an array regardless of how long it is
        if not isinstance(self.mthlim,list) and self.mthlim is not None:
            self.mthlim = [self.mthlim]
    
    # ========== Riemann solver library routines =============================   
    def list_riemann_solvers(self):
        r"""
        List available Riemann solvers 
        
        This routine returns a list of available Riemann solvers which is
        constructed in the Riemann solver package (:ref:`pyclaw_rp`).  In this 
        case it lists all Riemann solvers.
        
        :Output:
         - (list) - List of Riemann solver names valid to be used with
           :meth:`set_riemann_solver`
        
        .. note::
            These Riemann solvers are currently only accessible to the python 
            time stepping routines.
        """
        rp_solver_list = []
        
        # Construct list from each dimension list
        for rp_solver in rp_solver_list_1d:
            rp_solver_list.append('%s_1d' % rp_solver)
        for rp_solver in rp_solver_list_2d:
            rp_solver_list.append('%s_2d' % rp_solver)
        for rp_solver in rp_solver_list_3d:
            rp_solver_list.append('%s_3d' % rp_solver)
        
        return rp_solver_list
    
    def set_riemann_solver(self,solver_name):
        r"""
        Assigns the library solver solver_name as the Riemann solver.
        
        :Input:
         - *solver_name* - (string) Name of the solver to be used, raises a 
           NameError if the solver does not exist.
        """
        raise Exception("Cannot set a Riemann solver with this class," +
                                        " use one of the derived classes.")
         
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
         - *solutions* - (:class:`~pyclaw.solution.Solution`) Dictionary of 
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
            

# ============================================================================
#  ClawPack 1d Solver Class
# ============================================================================
class ClawSolver1D(ClawSolver):
    r"""
    Clawpack evolution routine in 1D
    
    This class represents the 1d clawpack solver on a single grid.  Note that 
    there are routines here for interfacing with the fortran time stepping 
    routines and the python time stepping routines.  The ones used are 
    dependent on the argument given to the initialization of the solver 
    (defaults to python).
    
    .. attribute:: rp
    
        Riemann solver function.
        
    :Initialization:
    
    Input:
     - *data* - (:class:`~pyclaw.data.Data`) An instance of a Data object whose
       parameters can be used to initialize this solver
    Output:
     - (:class:`ClawSolver1D`) - Initialized 1d clawpack solver
        
    :Authors:
        Kyle T. Mandli (2008-09-11) Initial version
    """

    def __init__(self,data=None):
        r"""
        Create 1d Clawpack solver
        
        See :class:`ClawSolver1D` for more info.
        """   
        
        # Add the functions as required attributes
        self._required_attrs.append('rp')
        self._default_attr_values['rp'] = None
        
        # Import Riemann solvers
        exec('import pyclaw.evolve.rp as rp',globals())
            
        super(ClawSolver1D,self).__init__(data)

    # ========== Riemann solver library routines =============================   
    def list_riemann_solvers(self):
        r"""
        List available Riemann solvers 
        
        This routine returns a list of available Riemann solvers which is
        constructed in the Riemann solver package (_pyclaw_rp).  In this case
        it lists only the 1D Riemann solvers.
        
        :Output:
         - (list) - List of Riemann solver names valid to be used with
           :meth:`set_riemann_solver`
        
        .. note::
            These Riemann solvers are currently only accessible to the python 
            time stepping routines.
        """
        return rp.rp_solver_list_1d
    
    def set_riemann_solver(self,solver_name):
        r"""
        Assigns the library solver solver_name as the Riemann solver.
        
        :Input:
         - *solver_name* - (string) Name of the solver to be used, raises a 
           ``NameError`` if the solver does not exist.
        """
        if solver_name in rp.rp_solver_list_1d:
            exec("self.rp = rp.rp_%s_1d" % solver_name)
        else:
            error_msg = 'Could not find Riemann solver with name %s' % solver_name
            logger.warning(error_msg)
            raise NameError(error_msg)


    # ========== Setting Boundary Conditions ==================================
    def set_bc_lower(self,grid):
        r"""
        
        """

        # User defined functions
        x = grid.dimensions[0]
        if x.mthbc_lower == 0:
            #self.user_bc_lower(grid,self,qbc)
            temp =1 #this is meaningless
        # Zero-order extrapolation
        elif x.mthbc_lower == 1:
            #qbc[:grid.mbc,...] = qbc[grid.mbc,...]
        # Periodic
            temp =1 #this is meaningless
        elif x.mthbc_lower == 2:
            #qbc[:grid.mbc,...] = qbc[-2*grid.mbc:-grid.mbc,...]
            # no of elmts = grid.mbc [0, .., mbc-1] [grid.n, ...,grid.n+grid.mbc-1]
            list_from = [i for i in xrange(grid.x.n, grid.x.n + grid.mbc)]
            list_to = [i for i in xrange(grid.mbc)]
            is_from = PETSc.IS().createGeneral(list_from)
            is_to = PETSc.IS().createGeneral(list_to )
        
            q_scatter = PETSc.Scatter().create(grid.gqVec, is_from, grid.gqVec, is_to)
            q_scatter.scatterBegin(grid.gqVec, grid.gqVec, False, PETSc.Scatter.Mode.FORWARD)	
 	    q_scatter.scatterEnd( grid.gqVec, grid.gqVec, False, PETSc.Scatter.Mode.FORWARD)
 	    q_scatter.destroy()
            
        # Solid wall bc
        elif x.mthbc_lower == 3:
            raise NotImplementedError("Solid wall upper boundary condition not implemented.")
        else:
            raise NotImplementedError("Boundary condition %s not implemented" % x.mthbc_lower)


        
        


    # ========== Python Homogeneous Step =====================================
    def homogeneous_step(self,solutions):
        r"""
        Take one time step on the homogeneous hyperbolic system

        Takes one time step of size dt on the hyperbolic system defined in the
        appropriate Riemann solver rp.

        :Input:
         - *solutions* - (:class:`~pyclaw.solution.Solution`) Solution that 
           will be evolved

        :Version: 1.0 (2009-07-01)
        """
    
        # Grid we will be working on
        grid = solutions['n'].grids[0]
        # Number of equations
        meqn = solutions['n'].meqn
        # Limiter to use in the pth family
        limiter = np.array(self.mthlim,ndmin=1)  
        # Q with appended boundary conditions
        # q = grid.qbc()

        
        grid.da.globalToLocal(grid.gqVec, grid.lqVec)

        
        q = grid.lqVec.getArray()
         
        local_n = q.size
        # Flux vector
        f = np.empty( (local_n, meqn) )
    
        dtdx = np.zeros( (local_n) )

        # Find local value for dt/dx
        if grid.capa is not None:
            dtdx = self.dt / (grid.d[0] * grid.capa)
        else:
            dtdx += self.dt/grid.d[0]



        q= np.reshape(q, (local_n,meqn)) #remove value
        #why still need reshaping while I've set the dof?
        
        
        
        


        
    
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

       
        #LL = grid.mbc - 1
        #UL = grid.mbc + grid.n[0] + 1

        # Is this should be anny different?
        if PETSc.Comm.getRank(PETSc.COMM_WORLD) == 0:
            LL =  1
            UL =  local_n - 1
        elif PETSc.Comm.getRank(PETSc.COMM_WORLD) == (PETSc.Comm.getSize(PETSc.COMM_WORLD) -1):
            LL = 1
            UL = local_n - 1
        else:
            LL = 1
            UL = local_n - 1
        

        

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
            
        # Reset q update
        #grid.q = q[grid.mbc:-grid.mbc][:]

        



        
        
        grid.lqVec.setArray(q)
        grid.da.localToGlobal(grid.lqVec, grid.gqVec)

        # this part to deal with bc, needs to be separated in a method and deals with different cases of bc, not periodic only
        #if PETSc.Comm.getRank(PETSc.COMM_WORLD) == (PETSc.Comm.getSize(PETSc.COMM_WORLD) -1 or PETSc.Comm.getRank(PETSc.COMM_WORLD) == 0):
            #communicate bc points from the end to bc points at the begining of process 0
            
            #grid.gqVec.setValues([0,1], [q[-2,0], q[-1,0]] )
            #grid.gqVec.assemblyBegin()
            #grid.gqVec.assemblyEnd()
        
        self.set_bc_lower(grid)
        grid.q = grid.gqVec.getArray()
        grid.q= np.reshape(grid.q, (grid.q.size,grid.meqn))
        
        
        
