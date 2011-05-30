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

from pyclaw.evolve.solver import Solver

import limiters
import riemann

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
        self._default_attr_values['mbc'] = 2
        self._default_attr_values['mthlim'] = [1]
        self._default_attr_values['order'] = 2
        self._default_attr_values['src_split'] = 0
        self._default_attr_values['fwave'] = False
        self._default_attr_values['src'] = src
        self._default_attr_values['start_step'] = start_step
        self._default_attr_values['kernel_language'] = 'Fortran'
        self._default_attr_values['verbosity'] = 0

        # Call general initialization function
        super(ClawSolver,self).__init__(data)
    
    # ========== Setup Routine ===============================================
    def setup(self,solutions):
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

        # Putting this here now for PetClaw.  We should think about the best way
        # to handle CFL communication.
        self.communicateCFL()

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
        exec('import riemann',globals())
            
        super(ClawSolver1D,self).__init__(data)


    # ========== Setup routine =============================   
    def setup(self,solutions):
        r"""
        See setup doc string in the super class.
        We are initializing (allocating) the working arrays needed by fortran kernels 
        in this routine. These arrays are passed in each call to the fortran kernel classic.
        """
        import numpy as np

        # Grid we will be working on
        grid = solutions['n'].grids[0]

        # Number of equations
        meqn,maux,mwaves,mbc = grid.meqn,grid.maux,self.mwaves,self.mbc
        mx = grid.q.shape[1]

        #Set up mthlim array
        if not isinstance(self.mthlim,list): self.mthlim=[self.mthlim]
        if len(self.mthlim)==1: self.mthlim = self.mthlim * self.mwaves
        if len(self.mthlim)!=self.mwaves:
            raise Exception('Length of solver.mthlim is not 1 nor is it equal to solver.mwaves')
 
        self.method =np.ones(7, dtype=int) # hardcoded 7
        self.method[0] = self.dt_variable  # fixed or adjustable timestep
        self.method[1] = self.order  # order of the method
        self.method[2] = 0  # Not used in 1D
        self.method[3] = self.verbosity
        self.method[4] = self.src_split  # src term
        if (grid.capa == None):
            self.method[5] = 0  
        else:
            self.method[5] = 1  
        self.method[6] = maux  # aux
 
        if(self.kernel_language == 'Fortran'):
            import classic1
            grid.set_cparam(classic1)
            self.f    = np.empty( (meqn,mx+2*mbc) )
            self.wave = np.empty( (meqn,mwaves,mx+2*mbc) )
            self.s    = np.empty( (mwaves,mx+2*mbc) )
            self.amdq = np.empty( (meqn,mx+2*mbc) )
            self.apdq = np.empty( (meqn,mx+2*mbc) )

    def teardown(self):
        if(self.kernel_language == 'Fortran'):
            import classic1
            del classic1


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
        return riemann.rp_solver_list_1d
    
    def set_riemann_solver(self,solver_name):
        r"""
        Assigns the library solver solver_name as the Riemann solver.
        
        :Input:
         - *solver_name* - (string) Name of the solver to be used, raises a 
           ``NameError`` if the solver does not exist.
        """
        import logging
        if solver_name in riemann.rp_solver_list_1d:
            exec("self.rp = riemann.rp_%s_1d" % solver_name)
        else:
            logger = logging.getLogger('solver')
            error_msg = 'Could not find Riemann solver with name %s -- perhaps you need to add the solver to riemann.__init__.py.' % solver_name
            logger.warning(error_msg)
            raise NameError(error_msg)

    # ========== Homogeneous Step =====================================
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
        import numpy as np

        grid = solutions['n'].grids[0]
        q = self.qbc(grid,grid)

        meqn,maux,mwaves,mbc,aux = grid.meqn,grid.maux,self.mwaves,self.mbc,grid.aux
          
        if(self.kernel_language == 'Fortran'):
            from classic1 import step1
            
            mx = grid.q.shape[1]
            dx,dt = grid.d[0],self.dt
            dtdx = np.zeros( (mx+2*mbc) ) + dt/dx
            
            if(aux == None): aux = np.empty( (maux,mx+2*mbc) )
        
       
            q,self.cfl = step1(mx,mbc,mx,q,aux,dx,dt,self.method,self.mthlim,self.f,self.wave,self.s,self.amdq,self.apdq,dtdx)

        elif(self.kernel_language == 'Python'):
 
            # Limiter to use in the pth family
            limiter = np.array(self.mthlim,ndmin=1)  
        
            dtdx = np.zeros( (2*self.mbc+grid.n[0]) )

            # Find local value for dt/dx
            if grid.capa is not None:
                dtdx = self.dt / (grid.d[0] * grid.capa)
            else:
                dtdx += self.dt/grid.d[0]
        
            # Solve Riemann problem at each interface
            q_l=q[:,:-1]
            q_r=q[:,1:]
            if grid.aux is not None:
                aux_l=grid.aux[:,:-1]
                aux_r=grid.aux[:,1:]
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
            LL = self.mbc - 1
            UL = self.mbc + grid.n[0] + 1 

            # Update q for Godunov update
            for m in xrange(meqn):
                q[m,LL:UL] -= dtdx[LL:UL]*apdq[m,LL-1:UL-1]
                q[m,LL-1:UL-1] -= dtdx[LL-1:UL-1]*amdq[m,LL-1:UL-1]
        
            # Compute maximum wave speed
            self.cfl = 0.0
            for mw in xrange(wave.shape[1]):
                smax1 = np.max(dtdx[LL:UL]*s[mw,LL-1:UL-1])
                smax2 = np.max(-dtdx[LL-1:UL-1]*s[mw,LL-1:UL-1])
                self.cfl = max(self.cfl,smax1,smax2)

            # If we are doing slope limiting we have more work to do
            if self.order == 2:
                # Initialize flux corrections
                f = np.zeros( (meqn,grid.n[0] + 2*self.mbc) )
            
                # Apply Limiters to waves
                if (limiter > 0).any():
                    wave = limiters.limit(grid.meqn,wave,s,limiter,dtdx)

                # Compute correction fluxes for second order q_{xx} terms
                dtdxave = 0.5 * (dtdx[LL-1:UL-1] + dtdx[LL:UL])
                if self.fwave:
                    for mw in xrange(wave.shape[1]):
                        sabs = np.abs(s[mw,LL-1:UL-1])
                        om = 1.0 - sabs*dtdxave[:UL-LL]
                        ssign = np.sign(s[mw,LL-1:UL-1])
                        for m in xrange(meqn):
                            f[m,LL:UL] += 0.5 * ssign * om * wave[m,mw,LL-1:UL-1]
                else:
                    for mw in xrange(wave.shape[1]):
                        sabs = np.abs(s[mw,LL-1:UL-1])
                        om = 1.0 - sabs*dtdxave[:UL-LL]
                        for m in xrange(meqn):
                            f[m,LL:UL] += 0.5 * sabs * om * wave[m,mw,LL-1:UL-1]

                # Update q by differencing correction fluxes
                for m in xrange(meqn):
                    q[m,LL:UL-1] -= dtdx[LL:UL-1] * (f[m,LL+1:UL] - f[m,LL:UL-1]) 

        else: raise Exception("Unrecognized kernel_language; choose 'Fortran' or 'Python'")
            
        grid.q = q[:,self.mbc:-self.mbc]

   

# ============================================================================
#  ClawPack 2d Solver Class
# ============================================================================
class ClawSolver2D(ClawSolver):
    r"""
    Clawpack evolution routine in 2D
    
    Note that only the fortran routines are supported for now in 2D.
    """

    def __init__(self,data=None):
        r"""
        Create 2d Clawpack solver
        
        See :class:`ClawSolver2D` for more info.
        """   
        
        # Add the functions as required attributes
        self._required_attrs.append('rp')
        self._default_attr_values['rp'] = None
        
        # Import Riemann solvers
        exec('import riemann',globals())
            
        self._default_attr_values['dim_split'] = True

        super(ClawSolver2D,self).__init__(data)

    # ========== Setup routine =============================   
    def setup(self,solutions):
        r"""
        See setup doc string in the super class.
        We are initializing (allocating) the working arrays needed by fortran kernels 
        in this routine. These arrays are passed in each call to the fortran kernel dimsp2.
        """

        #Set up mthlim array
        if not isinstance(self.mthlim,list): self.mthlim=[self.mthlim]
        if len(self.mthlim)==1: self.mthlim = self.mthlim * self.mwaves
        if len(self.mthlim)!=self.mwaves:
            raise Exception('Length of solver.mthlim is not 1 nor is it equal to solver.mwaves')
 

        if(self.kernel_language == 'Fortran'):
            import numpy as np

            # Grid we will be working on
            grid = solutions['n'].grids[0]

            # The reload here is necessary because otherwise the common block
            # cparam in the Riemann solver doesn't get flushed between running
            # different tests in a single Python session.
            import classic2
            reload(classic2)
            grid.set_cparam(classic2)

            # Number of equations
            meqn,maux,mwaves,mbc,aux = grid.meqn,grid.maux,self.mwaves,self.mbc,grid.aux
            maxmx,maxmy = grid.q.shape[1],grid.q.shape[2]
            maxm = max(maxmx, maxmy)

            #We ought to put method and cflv and many other things in a Fortran
            #module and set the fortran variables directly here.
            self.method =np.ones(7, dtype=int)
            self.method[0] = self.dt_variable
            self.method[1] = self.order
            if self.dim_split:
                self.method[2] = -1  # Godunov dimensional splitting
            else:
                self.method[2] = 1   # 
            self.method[3] = self.verbosity
            self.method[4] = self.src_split  # src term

            if (grid.capa == None): 
                self.method[5] = 0
            else: 
                self.method[5] = 1  
            self.method[6] = maux
                
            self.cflv = np.zeros(4)
            self.cflv[0:2] = [self.cfl_max,self.cfl_desired]
            #cflv[2] and cflv[3] are output values.

            #The following is a hack to work around an issue
            #with f2py.  It involves wastefully allocating a three arrays.
            #f2py seems not able to handle multiple zero-size arrays being passed.
            # it appears the bug is related to f2py/src/fortranobject.c line 841.
            if(aux == None): maux=1

            if self.src_split < 2: narray = 1
            else: narray = 2

            # These work arrays really ought to live inside a fortran module
            # as is done for sharpclaw
            self.qadd = np.empty((meqn,maxm+2*mbc))
            self.fadd = np.empty((meqn,maxm+2*mbc))
            self.gadd = np.empty((meqn,2,maxm+2*mbc))
            self.q1d  = np.empty((meqn,maxm+2*mbc))
            self.dtdx1d = np.empty((maxm+2*mbc))
            self.dtdy1d = np.empty((maxm+2*mbc))
            self.aux1 = np.empty((maux,maxm+2*mbc))
            self.aux2 = np.empty((maux,maxm+2*mbc))
            self.aux3 = np.empty((maux,maxm+2*mbc))
            mwork = (maxm+2*mbc) * (5*meqn + mwaves + meqn*mwaves) \
                  + (narray-1) * (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn
            # Amal: I think no need for the term
            # (narray-1) * (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn
            # this extra q array should be created and handled in function
            # step in case we have src term with strange splitting (Do not
            # think the fortran code will complain, but not sure)
            self.work = np.empty((mwork))

        else: raise Exception('Only Fortran kernels are supported in 2D.')


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
        return riemann.rp_solver_list_1d
    
    def set_riemann_solver(self,solver_name):
        r"""
        Assigns the library solver solver_name as the Riemann solver.
        
        :Input:
         - *solver_name* - (string) Name of the solver to be used, raises a 
           ``NameError`` if the solver does not exist.
        """
        import logging
        if solver_name in riemann.rp_solver_list_1d:
            exec("self.rp = riemann.rp_%s_1d" % solver_name)
        else:
            logger = logging.getLogger('solver')
            error_msg = 'Could not find Riemann solver with name %s' % solver_name
            logger.warning(error_msg)
            raise NameError(error_msg)

    # ========== Homogeneous Step =====================================
    def homogeneous_step(self,solutions):
        r"""
        
        """
        import numpy as np

        # Grid we will be working on
        grid = solutions['n'].grids[0]
        # Number of equations
        meqn,maux,mwaves,mbc,aux = grid.meqn,grid.maux,self.mwaves,self.mbc,grid.aux


        if(self.kernel_language == 'Fortran'):
            from classic2 import dimsp2, step2
            mx,my = grid.q.shape[1],grid.q.shape[2]
            maxm = max(mx,my)
            aux = grid.aux
            
            #The following is a hack to work around an issue
            #with f2py.  It involves wastefully allocating a three arrays.
            #f2py seems not able to handle multiple zero-size arrays being passed.
            # it appears the bug is related to f2py/src/fortranobject.c line 841.
            if(aux == None): 
                maux=1
                aux=np.empty((0,mx+2*mbc,my+2*mbc))
                
            dx,dy,dt = grid.d[0],grid.d[1],self.dt

            qold = self.qbc(grid,grid)
            qnew = qold.copy('F') #(input/output)

            if self.dim_split:
                q, cfl = dimsp2(maxm,mx,my,mbc,mx,my, \
                      qold,qnew,aux,dx,dy,dt,self.method,self.mthlim,self.cfl,self.cflv, \
                      self.qadd,self.fadd,self.gadd,self.q1d,self.dtdx1d,\
                      self.dtdy1d,self.aux1,self.aux2,self.aux3,self.work)
            else:
                q, cfl = step2(maxm,mx,my,mbc,mx,my, \
                      qold,qnew,aux,dx,dy,dt,self.method,self.mthlim,self.cfl, \
                      self.qadd,self.fadd,self.gadd,self.q1d,self.dtdx1d,\
                      self.dtdy1d,self.aux1,self.aux2,self.aux3,self.work)


            self.cfl = cfl
            grid.q=q[:,mbc:mx+mbc,mbc:my+mbc]

        else:
            raise NotImplementedError("No python implementation for homogeneous_step in case of 2D.")
