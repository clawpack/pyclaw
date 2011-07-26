r"""
Module containing the classic Clawpack solvers.

This module contains the pure and wrapped classic clawpack solvers.  All 
clawpack solvers inherit from the :class:`ClawSolver` superclass which in turn 
inherits from the :class:`~pyclaw.solver.Solver` superclass.  These
are both pure virtual classes; the only solver classes that should be instantiated
are the dimension-specific ones, :class:`ClawSolver1D` and :class:`ClawSolver2D`.

:Authors:
    Kyle T. Mandli (2008-09-11) Initial version
    Amal Alghamdi (2010-2011)   Wrapped Fortran routines
    David I. Ketcheson (2011)   Various refinements, including removing BCs from grid
"""

from pyclaw.solver import Solver

import limiters.tvd
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
    
        Limiter to be used on each wave.  ``Default = limiters.tvd.minmod``
    
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
        
    .. attribute:: kernel_language

        Specifies whether to use wrapped Fortran routines ('Fortran')
        or pure Python ('Python').  ``Default = 'Fortran'``.
    
    .. attribute:: verbosity

        The level of detail of logged messages from the Fortran solver.
        ``Default = 0``.

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
        for attr in ['limiters','order','src_split','fwave','src','start_step']:
            self._required_attrs.append(attr)
        
        # Default required attributes
        self._default_attr_values['mbc'] = 2
        self._default_attr_values['limiters'] = limiters.tvd.minmod
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
        :class:`~pyclaw.controller.Controller`.
        """
        pass 
    
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
        pyclaw.solver.Solver superclass.

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

    def set_mthlim(self):
        self.mthlim = self.limiters
        if not isinstance(self.limiters,list): self.mthlim=[self.mthlim]
        if len(self.mthlim)==1: self.mthlim = self.mthlim * self.mwaves
        if len(self.mthlim)!=self.mwaves:
            raise Exception('Length of solver.limiters is not equal to 1 or to solver.mwaves')
 
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
            
        self.ndim = 1

        super(ClawSolver1D,self).__init__(data)


    # ========== Setup routine =============================   
    def setup(self,solutions):
        r"""
        Perform essential solver setup.  This routine must be called before
        solver.step() may be called.
        """
        self.set_mthlim()
        if(self.kernel_language == 'Fortran'):
            self.set_fortran_parameters(solutions)


    def set_fortran_parameters(self,solutions):
        r"""
        Pack parameters into format recognized by Clawpack (Fortran) code.

        Sets the method array and the cparam common block for the Riemann solver.
        """
        import numpy as np

        state = solutions['n'].state

        self.method =np.ones(7, dtype=int)
        self.method[0] = self.dt_variable
        self.method[1] = self.order 
        self.method[2] = 0  # Not used in 1D
        self.method[3] = self.verbosity
        self.method[4] = self.src_split
        if (state.capa == None):
            self.method[5] = 0  
        else:
            self.method[5] = 1  
        self.method[6] = state.maux  # aux
 
        if self.fwave:
            import classic1fw as classic1
        else:
            import classic1
        state.set_cparam(classic1)


    def teardown(self):
        r"""
        Delete Fortran objects, which otherwise tend to persist in Python sessions.
        """
        if(self.kernel_language == 'Fortran'):
            if self.fwave:
                import classic1fw as classic1
            else:
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

        state = solutions['n'].states[0]
        grid = solutions['n'].grids[0]

        q = state.qbc
            
        meqn,maux,mwaves,mbc = state.meqn,state.maux,self.mwaves,self.mbc
          
        if maux>0:
            aux = self.auxbc(grid,state)

        if(self.kernel_language == 'Fortran'):
            if self.fwave:
                import classic1fw as classic1
            else:
                import classic1

            mx = grid.n[0]
            dx,dt = grid.d[0],self.dt
            dtdx = np.zeros( (mx+2*mbc) ) + dt/dx
            
            if(maux == 0): aux = np.empty( (maux,mx+2*mbc) )
        
       
            q,self.cfl = classic1.step1(mx,mbc,mx,q,aux,dx,dt,self.method,self.mthlim)

        elif(self.kernel_language == 'Python'):
 
            # Limiter to use in the pth family
            limiter = np.array(self.mthlim,ndmin=1)  
        
            dtdx = np.zeros( (2*self.mbc+grid.n[0]) )

            # Find local value for dt/dx
            if state.capa is not None:
                dtdx = self.dt / (grid.d[0] * state.capa)
            else:
                dtdx += self.dt/grid.d[0]
        
            # Solve Riemann problem at each interface
            q_l=q[:,:-1]
            q_r=q[:,1:]
            if state.aux is not None:
                aux_l=aux[:,:-1]
                aux_r=aux[:,1:]
            else:
                aux_l = None
                aux_r = None
            wave,s,amdq,apdq = self.rp(q_l,q_r,aux_l,aux_r,state.aux_global)
            
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
                    wave = limiters.tvd.limit(state.meqn,wave,s,limiter,dtdx)

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
        # Amal: this line need to be replcaed by set_global_q    
        state.q = q[:,self.mbc:-self.mbc]

   

# ============================================================================
#  ClawPack 2d Solver Class
# ============================================================================
class ClawSolver2D(ClawSolver):
    r"""
    2D Classic (Clawpack) solver.

    Solve using the wave propagation algorithms of Randy LeVeque's
    Clawpack code (www.clawpack.org).

    See also the documentation for ClawSolver1D.
    In addition to the attributes of ClawSolver1D, ClawSolver2D
    also has the following options:
    
    .. attribute:: dim_split
    
        If True, use dimensional splitting (Godunov splitting).
        Dimensional splitting with Strang splitting is not supported
        at present but could easily be enabled if necessary.
        If False, use unsplit Clawpack algorithms, possibly including
        transverse Riemann solves.

    .. attribute:: order_trans
    
        If dim_split is True, this option has no effect.  If
        dim_plit is False, then order_trans should be one of
        the following values:

        ClawSolver2D.no_trans: Transverse Riemann solver
        not used.  The stable CFL for this algorithm is 0.5.  Not recommended.
        
        ClawSolver2D.trans_inc: Transverse increment waves are computed
        and propagated.

        ClawSolver2D.trans_cor: Transverse increment waves and transverse
        correction waves are computed and propagated.

    Note that only the fortran routines are supported for now in 2D.
    """

    no_trans  = 0
    trans_inc = 1
    trans_cor = 2

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
        self._default_attr_values['order_trans'] = self.trans_inc

        self.ndim = 2

        super(ClawSolver2D,self).__init__(data)

    # ========== Setup routine =============================   
    def setup(self,solutions):
        r"""
        Perform essential solver setup.  This routine must be called before
        solver.step() may be called.
        """
        self.set_mthlim()

        if (not self.dim_split) and (self.order_trans==0):
            cfl_recommended = 0.5
        else:
            cfl_recommended = 1.0

        if self.cfl_max > cfl_recommended:
            import warnings
            warnings.warn('cfl_max is set higher than the recommended value of %s' % cfl_recommended)

        if(self.kernel_language == 'Fortran'):
            self.set_fortran_parameters(solutions)
        else: raise Exception('Only Fortran kernels are supported in 2D.')


    def set_fortran_parameters(self,solutions):
        r"""
        Pack parameters into format recognized by Clawpack (Fortran) code.

        Sets the method array and the cparam common block for the Riemann solver.
        """
        import numpy as np

        # Grid we will be working on
        state = solutions['n'].states[0]
        grid  = solutions['n'].grids[0]

        # The reload here is necessary because otherwise the common block
        # cparam in the Riemann solver doesn't get flushed between running
        # different tests in a single Python session.
        if self.fwave:
            import classic2fw as classic2
        else:
            import classic2
        reload(classic2)
        state.set_cparam(classic2)

        # Number of equations
        meqn,maux,mwaves,mbc,aux = state.meqn,state.maux,self.mwaves,self.mbc,state.aux
        maxmx,maxmy = grid.n[0],grid.n[1]
        maxm = max(maxmx, maxmy)

        #We ought to put method and cflv and many other things in a Fortran
        #module and set the fortran variables directly here.
        self.method =np.ones(7, dtype=int)
        self.method[0] = self.dt_variable
        self.method[1] = self.order
        if self.dim_split:
            self.method[2] = -1  # Godunov dimensional splitting
        else:
            self.method[2] = self.order_trans
        self.method[3] = self.verbosity
        self.method[4] = self.src_split  # src term

        if (state.capa == None): 
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
        Take a step on the homogeneous hyperbolic system using the Clawpack
        algorithm.

        Clawpack is based on the Lax-Wendroff method, combined with Riemann
        solvers and TVD limiters applied to waves.
        """
        import numpy as np


        if(self.kernel_language == 'Fortran'):
            state = solutions['n'].states[0]
            grid  = solutions['n'].grids[0]
            meqn,maux,mwaves,mbc = state.meqn,state.maux,self.mwaves,self.mbc
            mx,my = grid.n[0],grid.n[1]
            maxm = max(mx,my)
            
            #The following is a hack to work around an issue
            #with f2py.  It involves wastefully allocating a three arrays.
            #f2py seems not able to handle multiple zero-size arrays being passed.
            # it appears the bug is related to f2py/src/fortranobject.c line 841.
            if maux == 0: 
                aux=np.empty((0,mx+2*mbc,my+2*mbc))
            else:
                aux = self.auxbc(grid,state)

            
            dx,dy,dt = grid.d[0],grid.d[1],self.dt
            
            qnew = state.qbc #(input/output)
            if self.dt_variable:
                qold = self.qbc_backup # Solver should quarantee that 
                                        # qbc_backup will not be
                                        # changed so that it can be used in
                                        # case of step rejection.
            else:
                qold = qnew.copy('F')
            
            if self.fwave:
                import classic2fw as classic2
            else:
                import classic2

            if self.dim_split:
                q, cfl = classic2.dimsp2(maxm,mx,my,mbc,mx,my, \
                      qold,qnew,aux,dx,dy,dt,self.method,self.mthlim,self.cfl,self.cflv, \
                      self.aux1,self.aux2,self.aux3,self.work)
            else:
                q, cfl = classic2.step2(maxm,mx,my,mbc,mx,my, \
                      qold,qnew,aux,dx,dy,dt,self.method,self.mthlim,self.cfl, \
                      self.aux1,self.aux2,self.aux3,self.work)

            self.cfl = cfl
            self.set_global_q(state, q)

        else:
            raise NotImplementedError("No python implementation for homogeneous_step in case of 2D.")
