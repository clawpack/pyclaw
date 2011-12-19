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

# ============================================================================
#  Generic Clawpack solver class
# ============================================================================
class ClawSolver(Solver):
    r"""
    Generic classic Clawpack solver
    
    All Clawpack solvers inherit from this base class.
    
    .. attribute:: mthlim 
    
        Limiter(s) to be used.  Specified either as one value or a list.
        If one value, the specified limiter is used for all wave families.
        If a list, the specified values indicate which limiter to apply to
        each wave family.  Take a look at pyclaw.limiters.tvd for an enumeration.
        ``Default = limiters.tvd.minmod``
    
    .. attribute:: order
    
        Order of the solver, either 1 for first order (i.e., Godunov's method)
        or 2 for second order (Lax-Wendroff-LeVeque).
        ``Default = 2``
    
    .. attribute:: src_split
    
        Which source splitting method to use: 1 for first 
        order Godunov splitting and 2 for second order Strang splitting.
        ``Default = 1``
        
    .. attribute:: fwave
    
        Whether to split the flux jump (rather than the jump in Q) into waves; 
        requires that the Riemann solver performs the splitting.  
        ``Default = False``
        
    .. attribute:: step_src
    
        Handle for function that evaluates the source term.  
        The required signature for this function is:

        def step_src(solver,state,dt)
    
    .. attribute:: start_step
    
        Function called before each time step is taken.
        The required signature for this function is:
        
        def start_step(solver,solution)

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
    """
    # ========== Generic Init Routine ========================================
    def __init__(self,data=None):
        r"""
        See :class:`ClawSolver` for full documentation.
        """
        self.mbc = 2
        self.limiters = limiters.tvd.minmod
        self.order = 2
        self.src_split = 1
        self.fwave = False
        self.step_src = None
        self.start_step = None
        self.kernel_language = 'Fortran'
        self.verbosity = 0
        self.cfl_max = 1.0
        self.cfl_desired = 0.9

        # Call general initialization function
        super(ClawSolver,self).__init__(data)
    
    # ========== Time stepping routines ======================================
    def step(self,solution):
        r"""
        Evolve solution one time step

        The elements of the algorithm for taking one step are:
        
        1. The :meth:`start_step` function is called
        
        2. A half step on the source term :func:`step_src` if Strang splitting is 
           being used (:attr:`src_split` = 2)
        
        3. A step on the homogeneous problem :math:`q_t + f(q)_x = 0` is taken
        
        4. A second half step or a full step is taken on the source term
           :func:`step_src` depending on whether Strang splitting was used 
           (:attr:`src_split` = 2) or Godunov splitting 
           (:attr:`src_split` = 1)

        This routine is called from the method evolve_to_time defined in the
        pyclaw.solver.Solver superclass.

        :Input:
         - *solution* - (:class:`~pyclaw.solution.Solution`) solution to be evolved
         
        :Output: 
         - (bool) - True if full step succeeded, False otherwise
        """

        if self.start_step is not None:
            self.start_step(self,solution)

        if self.src_split == 2 and self.step_src is not None:
            self.step_src(self,solution.states[0],self.dt/2.0)
    
        self.step_hyperbolic(solution)

        # Check here if the CFL condition is satisfied. 
        # If not, return # immediately to evolve_to_time and let it deal with
        # picking a new step size (dt).
        if self.cfl.get_cached_max() >= self.cfl_max:
            return False

        if self.step_src is not None:
            # Strang splitting
            if self.src_split == 2:
                self.step_src(self,solution.states[0],self.dt/2.0)

            # Godunov Splitting
            if self.src_split == 1:
                self.step_src(self,solution.states[0],self.dt)
                
        return True
            
    def check_cfl_settings(self):
        pass

    def allocate_workspace(self,solution):
        pass

    def step_hyperbolic(self,solution):
        r"""
        Take one homogeneous step on the solution.
        
        This is a dummy routine and must be overridden.
        """
        raise Exception("Dummy routine, please override!")

    def set_mthlim(self):
        r"""
        Convenience routine to convert user's limiter specification to 
        the format understood by the Fortran code (i.e., a list of length mwaves).
        """
        self.mthlim = self.limiters
        if not isinstance(self.limiters,list): self.mthlim=[self.mthlim]
        if len(self.mthlim)==1: self.mthlim = self.mthlim * self.mwaves
        if len(self.mthlim)!=self.mwaves:
            raise Exception('Length of solver.limiters is not equal to 1 or to solver.mwaves')
 
    def set_method(self,state):
        r"""
        Set values of the method() array required by the Fortran code.
        These are algorithmic parameters.
        """
        import numpy as np
        #We ought to put method and many other things in a Fortran
        #module and set the fortran variables directly here.
        self.method =np.empty(7, dtype=int,order='F')
        self.method[0] = self.dt_variable
        self.method[1] = self.order
        if self.ndim==1:
            self.method[2] = 0  # Not used in 1D
        elif self.dim_split:
            self.method[2] = -1  # First-order dimensional splitting
        else:
            self.method[2] = self.order_trans
        self.method[3] = self.verbosity
        self.method[4] = 0  # Not used for PyClaw (would be self.src_split)
        self.method[5] = state.mcapa + 1
        self.method[6] = state.maux

    def setup(self,solution):
        r"""
        Perform essential solver setup.  This routine must be called before
        solver.step() may be called.
        """
        if(self.kernel_language == 'Fortran'):
            # Set name of expected .so file with compile Fortran routines
            self.so_name = 'classic'+str(self.ndim)
            if self.fwave: self.so_name+='fw'

        # This is a hack to deal with the fact that petsc4py
        # doesn't allow us to change the stencil_width (mbc)
        solution.state.set_mbc(self.mbc)
        # End hack

        self.check_cfl_settings()

        self.set_mthlim()
        if(self.kernel_language == 'Fortran'):
            self.set_fortran_parameters(solution)
            self.allocate_workspace(solution)
        elif self.ndim>1:
            raise Exception('Only Fortran kernels are supported in multi-D.')

        self.allocate_bc_arrays(solution.states[0])


    def set_fortran_parameters(self,solution):
        r"""
        Pack parameters into format recognized by Clawpack (Fortran) code.

        Sets the method array and the cparam common block for the Riemann solver.
        """
        self.set_method(solution.state)
        classic = __import__(self.so_name)
        # The reload here is necessary because otherwise the common block
        # cparam in the Riemann solver doesn't get flushed between running
        # different tests in a single Python session.
        reload(classic)
        solution.state.set_cparam(classic)

    def teardown(self):
        r"""
        Delete Fortran objects, which otherwise tend to persist in Python sessions.
        """
        if(self.kernel_language == 'Fortran'):
            classic = __import__(self.so_name)
            del classic



# ============================================================================
#  ClawPack 1d Solver Class
# ============================================================================
class ClawSolver1D(ClawSolver):
    r"""
    Clawpack evolution routine in 1D
    
    This class represents the 1d clawpack solver on a single grid.  Note that 
    there are routines here for interfacing with the fortran time stepping 
    routines and the Python time stepping routines.  The ones used are 
    dependent on the argument given to the initialization of the solver 
    (defaults to python).
    
    :Initialization:
    
    Input:
     - *data* - (:class:`~pyclaw.data.Data`) An instance of a Data object whose
       parameters can be used to initialize this solver
    Output:
     - (:class:`ClawSolver1D`) - Initialized 1d clawpack solver
    """

    def __init__(self,data=None):
        r"""
        Create 1d Clawpack solver
        
        See :class:`ClawSolver1D` for more info.
        """   
        self.ndim = 1

        super(ClawSolver1D,self).__init__(data)


    # ========== Homogeneous Step =====================================
    def step_hyperbolic(self,solution):
        r"""
        Take one time step on the homogeneous hyperbolic system.

        :Input:
         - *solution* - (:class:`~pyclaw.solution.Solution`) Solution that 
           will be evolved
        """
        import numpy as np

        state = solution.states[0]
        grid = state.grid

        self.apply_q_bcs(state)
            
        meqn,mbc = state.meqn,self.mbc
          
        if(self.kernel_language == 'Fortran'):
            classic = __import__(self.so_name)

            mx = grid.ng[0]
            dx,dt = grid.d[0],self.dt
            dtdx = np.zeros( (mx+2*mbc) ) + dt/dx
            
            self.qbc,cfl = classic.step1(mbc,mx,self.qbc,self.auxbc,dx,dt,self.method,self.mthlim)
            
        elif(self.kernel_language == 'Python'):
 
            q   = self.qbc
            aux = self.auxbc
            # Limiter to use in the pth family
            limiter = np.array(self.mthlim,ndmin=1)  
        
            dtdx = np.zeros( (2*self.mbc+grid.ng[0]) )

            # Find local value for dt/dx
            if state.mcapa>=0:
                dtdx = self.dt / (grid.d[0] * state.aux[mcapa,:])
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
            UL = self.mbc + grid.ng[0] + 1 

            # Update q for Godunov update
            for m in xrange(meqn):
                q[m,LL:UL] -= dtdx[LL:UL]*apdq[m,LL-1:UL-1]
                q[m,LL-1:UL-1] -= dtdx[LL-1:UL-1]*amdq[m,LL-1:UL-1]
        
            # Compute maximum wave speed
            cfl = 0.0
            for mw in xrange(wave.shape[1]):
                smax1 = np.max(dtdx[LL:UL]*s[mw,LL-1:UL-1])
                smax2 = np.max(-dtdx[LL-1:UL-1]*s[mw,LL-1:UL-1])
                cfl = max(cfl,smax1,smax2)

            # If we are doing slope limiting we have more work to do
            if self.order == 2:
                # Initialize flux corrections
                f = np.zeros( (meqn,grid.ng[0] + 2*self.mbc) )
            
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

        self.cfl.update_global_max(cfl)
        state.set_q_from_qbc(mbc,self.qbc)
   

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
        self.dim_split = True
        self.order_trans = self.trans_inc

        self.ndim = 2

        super(ClawSolver2D,self).__init__(data)

    def check_cfl_settings(self):
        if (not self.dim_split) and (self.order_trans==0):
            cfl_recommended = 0.5
        else:
            cfl_recommended = 1.0

        if self.cfl_max > cfl_recommended:
            import warnings
            warnings.warn('cfl_max is set higher than the recommended value of %s' % cfl_recommended)
            warnings.warn(str(self.cfl_desired))


    def allocate_workspace(self,solution):
        r"""
        Pack parameters into format recognized by Clawpack (Fortran) code.

        Sets the method array and the cparam common block for the Riemann solver.
        """
        import numpy as np

        state = solution.state

        meqn,maux,mwaves,mbc,aux = state.meqn,state.maux,self.mwaves,self.mbc,state.aux

        #The following is a hack to work around an issue
        #with f2py.  It involves wastefully allocating three arrays.
        #f2py seems not able to handle multiple zero-size arrays being passed.
        # it appears the bug is related to f2py/src/fortranobject.c line 841.
        if(aux == None): maux=1

        grid  = state.grid
        maxmx,maxmy = grid.ng[0],grid.ng[1]
        maxm = max(maxmx, maxmy)

        # These work arrays really ought to live inside a fortran module
        # as is done for sharpclaw
        self.aux1 = np.empty((maux,maxm+2*mbc),order='F')
        self.aux2 = np.empty((maux,maxm+2*mbc),order='F')
        self.aux3 = np.empty((maux,maxm+2*mbc),order='F')
        mwork = (maxm+2*mbc) * (5*meqn + mwaves + meqn*mwaves)
        self.work = np.empty((mwork),order='F')


    # ========== Hyperbolic Step =====================================
    def step_hyperbolic(self,solution):
        r"""
        Take a step on the homogeneous hyperbolic system using the Clawpack
        algorithm.

        Clawpack is based on the Lax-Wendroff method, combined with Riemann
        solvers and TVD limiters applied to waves.
        """
        import numpy as np


        if(self.kernel_language == 'Fortran'):
            state = solution.states[0]
            grid = state.grid
            dx,dy = grid.d
            mx,my = grid.ng
            maxm = max(mx,my)
            
            self.apply_q_bcs(state)
            qnew = self.qbc
            qold = qnew.copy('F')
            
            classic = __import__(self.so_name)

            if self.dim_split:
                #Right now only Godunov-dimensional-splitting is implemented.
                #Strang-dimensional-splitting could be added following dimsp2.f in Clawpack.

                q, cfl_x = classic.step2ds(maxm,self.mbc,mx,my, \
                      qold,qnew,self.auxbc,dx,dy,self.dt,self.method,self.mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work,1)

                q, cfl_y = classic.step2ds(maxm,self.mbc,mx,my, \
                      q,q,self.auxbc,dx,dy,self.dt,self.method,self.mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work,2)

                cfl = max(cfl_x,cfl_y)

            else:

                q, cfl = classic.step2(maxm,self.mbc,mx,my, \
                      qold,qnew,self.auxbc,dx,dy,self.dt,self.method,self.mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work)

            self.cfl.update_global_max(cfl)
            state.set_q_from_qbc(self.mbc,self.qbc)

        else:
            raise NotImplementedError("No python implementation for step_hyperbolic in 2D.")

# ============================================================================
#  ClawPack 3d Solver Class
# ============================================================================
class ClawSolver3D(ClawSolver):
    r"""
    3D Classic (Clawpack) solver.

    Solve using the wave propagation algorithms of Randy LeVeque's
    Clawpack code (www.clawpack.org).

    See also the documentation for ClawSolver1D.
    In addition to the attributes of ClawSolver1D, ClawSolver3D
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

        ClawSolver3D.no_trans: Transverse Riemann solver
        not used.  The stable CFL for this algorithm is 0.5.  Not recommended.
        
        ClawSolver3D.trans_inc: Transverse increment waves are computed
        and propagated.

        ClawSolver3D.trans_cor: Transverse increment waves and transverse
        correction waves are computed and propagated.

    Note that only Fortran routines are supported for now in 3D --
    there is no pure-python version.
    """

    no_trans  = 0
    trans_inc = 11
    trans_cor = 22

    def __init__(self,data=None):
        r"""
        Create 3d Clawpack solver
        
        See :class:`ClawSolver3D` for more info.
        """   
        # Add the functions as required attributes
        self.dim_split = True
        self.order_trans = self.trans_cor

        self.ndim = 3

        super(ClawSolver3D,self).__init__(data)

    # ========== Setup routine =============================   
    def allocate_workspace(self,solution):
        r"""
        Allocate auxN and work arrays for use in Fortran subroutines.
        """
        import numpy as np

        state = solution.states[0]

        meqn,maux,mwaves,mbc,aux = state.meqn,state.maux,self.mwaves,self.mbc,state.aux

        #The following is a hack to work around an issue
        #with f2py.  It involves wastefully allocating three arrays.
        #f2py seems not able to handle multiple zero-size arrays being passed.
        # it appears the bug is related to f2py/src/fortranobject.c line 841.
        if(aux == None): maux=1

        grid  = state.grid
        maxmx,maxmy = grid.ng[0],grid.ng[1]
        maxm = max(maxmx, maxmy)

        # These work arrays really ought to live inside a fortran module
        # as is done for sharpclaw
        self.aux1 = np.empty((maux,maxm+2*mbc,3),order='F')
        self.aux2 = np.empty((maux,maxm+2*mbc,3),order='F')
        self.aux3 = np.empty((maux,maxm+2*mbc,3),order='F')
        mwork = (maxm+2*mbc) * (31*meqn + mwaves + meqn*mwaves)
        self.work = np.empty((mwork),order='F')


    # ========== Hyperbolic Step =====================================
    def step_hyperbolic(self,solution):
        r"""
        Take a step on the homogeneous hyperbolic system using the Clawpack
        algorithm.

        Clawpack is based on the Lax-Wendroff method, combined with Riemann
        solvers and TVD limiters applied to waves.
        """
        import numpy as np


        if(self.kernel_language == 'Fortran'):
            state = solution.states[0]
            grid = state.grid
            dx,dy,dz = grid.d
            mx,my,mz = grid.ng
            maxm = max(mx,my,mz)
            
            self.apply_q_bcs(state)
            qnew = self.qbc
            qold = qnew.copy('F')
            
            classic = __import__(self.so_name)

            if self.dim_split:
                #Right now only Godunov-dimensional-splitting is implemented.
                #Strang-dimensional-splitting could be added following dimsp2.f in Clawpack.

                q, cfl_x = classic.step3ds(maxm,self.mbc,mx,my,mz, \
                      qold,qnew,self.auxbc,dx,dy,dz,self.dt,self.method,self.mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work,1)

                q, cfl_y = classic.step3ds(maxm,self.mbc,mx,my,mz, \
                      q,q,self.auxbc,dx,dy,dz,self.dt,self.method,self.mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work,2)

                q, cfl_z = classic.step3ds(maxm,self.mbc,mx,my,mz, \
                      q,q,self.auxbc,dx,dy,dz,self.dt,self.method,self.mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work,3)

                cfl = max(cfl_x,cfl_y,cfl_z)

            else:

                q, cfl = classic.step3(maxm,self.mbc,mx,my,mz, \
                      qold,qnew,self.auxbc,dx,dy,dz,self.dt,self.method,self.mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work)

            self.cfl.update_global_max(cfl)
            state.set_q_from_qbc(self.mbc,self.qbc)

        else:
            raise NotImplementedError("No python implementation for step_hyperbolic in 3D.")
