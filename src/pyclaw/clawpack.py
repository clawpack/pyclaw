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
    
        Whether to use a source splitting method, 0 for none, 1 for first 
        order Godunov splitting and 2 for second order Strang splitting.
        ``Default = 0``
        
    .. attribute:: fwave
    
        Whether to split the flux jump (rather than the jump in Q) into waves; 
        requires that the Riemann solver performs the splitting.  
        ``Default = False``
        
    .. attribute:: src
    
        Handle for function that evaluates the source term.  
        The required signature for this function is:

        def src(solver,solution,t,dt)
    
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
        
        # Required attributes for this solver
        for attr in ['limiters','order','src_split','fwave','src','start_step']:
            self._required_attrs.append(attr)
        
        # Default required attributes
        self._default_attr_values['mbc'] = 2
        self._default_attr_values['limiters'] = limiters.tvd.minmod
        self._default_attr_values['order'] = 2
        self._default_attr_values['src_split'] = 0
        self._default_attr_values['fwave'] = False
        self._default_attr_values['src'] = None
        self._default_attr_values['start_step'] = None
        self._default_attr_values['kernel_language'] = 'Fortran'
        self._default_attr_values['verbosity'] = 0

        # Call general initialization function
        super(ClawSolver,self).__init__(data)
    
    # ========== Time stepping routines ======================================
    def step(self,solution):
        r"""
        Evolve solution one time step

        The elements of the algorithm for taking one step are:
        
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
         - *solution* - (:class:`~pyclaw.solution.Solution`) solution to be evolved
         
        :Output: 
         - (bool) - True if full step succeeded, False otherwise
        """

        if self.start_step is not None:
            self.start_step(self,solution)

        if self.src_split == 2:
            self.src(self,solution,solution.t, self.dt/2.0)
    
        self.homogeneous_step(solution)

        # Check here if we violated the CFL condition, if we did, return 
        # immediately to evolve_to_time and let it deal with picking a new
        # dt
        if self.cfl >= self.cfl_max:
            return False

        # Strang splitting
        if self.src_split == 2:
            self.src(self,solution,solution.t + self.dt/2.0, self.dt/2.0)

        # Godunov Splitting
        if self.src_split == 1:
            self.src(self,solution,solution.t,self.dt)
            
        return True
            
    def homogeneous_step(self,solution):
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
        
        self.ndim = 1

        super(ClawSolver1D,self).__init__(data)


    # ========== Setup routine =============================   
    def setup(self,solution):
        r"""
        Perform essential solver setup.  This routine must be called before
        solver.step() may be called.
        """
        self.set_mthlim()
        if(self.kernel_language == 'Fortran'):
            self.set_fortran_parameters(solution)

        self.allocate_bc_arrays(solution.states[0])

    def set_fortran_parameters(self,solution):
        r"""
        Pack parameters into format recognized by Clawpack (Fortran) code.

        Sets the method array and the cparam common block for the Riemann solver.
        """
        import numpy as np

        state = solution.state

        self.method =np.ones(7, dtype=int)
        self.method[0] = self.dt_variable
        self.method[1] = self.order 
        self.method[2] = 0  # Not used in 1D
        self.method[3] = self.verbosity
        self.method[4] = self.src_split
        self.method[5] = state.mcapa + 1
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


    # ========== Homogeneous Step =====================================
    def homogeneous_step(self,solution):
        r"""
        Take one time step on the homogeneous hyperbolic system.

        :Input:
         - *solution* - (:class:`~pyclaw.solution.Solution`) Solution that 
           will be evolved
        """
        import numpy as np

        state = solution.states[0]
        grid = state.grid

        q = self.qbc
            
        meqn,maux,mwaves,mbc = state.meqn,state.maux,self.mwaves,self.mbc
          
        if(self.kernel_language == 'Fortran'):
            if self.fwave:
                import classic1fw as classic1
            else:
                import classic1

            mx = grid.ng[0]
            dx,dt = grid.d[0],self.dt
            dtdx = np.zeros( (mx+2*mbc) ) + dt/dx
            
            q,self.cfl = classic1.step1(mx,mbc,mx,q,self.auxbc,dx,dt,self.method,self.mthlim)

        elif(self.kernel_language == 'Python'):
 
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
            self.cfl = 0.0
            for mw in xrange(wave.shape[1]):
                smax1 = np.max(dtdx[LL:UL]*s[mw,LL-1:UL-1])
                smax2 = np.max(-dtdx[LL-1:UL-1]*s[mw,LL-1:UL-1])
                self.cfl = max(self.cfl,smax1,smax2)

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
        # Amal: this line need to be replaced by set_global_q    
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
    def setup(self,solution):
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
            self.set_fortran_parameters(solution)
        else: raise Exception('Only Fortran kernels are supported in 2D.')

        self.allocate_bc_arrays(solution.states[0])

    def set_fortran_parameters(self,solution):
        r"""
        Pack parameters into format recognized by Clawpack (Fortran) code.

        Sets the method array and the cparam common block for the Riemann solver.
        """
        import numpy as np

        # Grid we will be working on
        state = solution.states[0]

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

        #We ought to put method and cflv and many other things in a Fortran
        #module and set the fortran variables directly here.
        self.method =np.empty(7, dtype=int,order='F')
        self.method[0] = self.dt_variable
        self.method[1] = self.order
        if self.dim_split:
            self.method[2] = -1  # Godunov dimensional splitting
        else:
            self.method[2] = self.order_trans
        self.method[3] = self.verbosity
        self.method[4] = self.src_split  # src term
        self.method[5] = state.mcapa + 1
        self.method[6] = state.maux
            
        if self.dim_split:
            # We should probably interface directly to step2ds (rather than dimsp2)
            # and avoid the need for cflv
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

        grid  = state.grid
        maxmx,maxmy = grid.ng[0],grid.ng[1]
        maxm = max(maxmx, maxmy)

        # These work arrays really ought to live inside a fortran module
        # as is done for sharpclaw
        self.aux1 = np.empty((maux,maxm+2*mbc),order='F')
        self.aux2 = np.empty((maux,maxm+2*mbc),order='F')
        self.aux3 = np.empty((maux,maxm+2*mbc),order='F')
        mwork = (maxm+2*mbc) * (5*meqn + mwaves + meqn*mwaves) \
              + (narray-1) * (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn
        # Amal: I think no need for the term
        # (narray-1) * (maxmx + 2*mbc) * (maxmy + 2*mbc) * meqn
        # this extra q array should be created and handled in function
        # step in case we have src term with strange splitting (Do not
        # think the fortran code will complain, but not sure)
        self.work = np.empty((mwork),order='F')


    # ========== Homogeneous Step =====================================
    def homogeneous_step(self,solution):
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
            meqn,maux,mwaves,mbc = state.meqn,state.maux,self.mwaves,self.mbc
            mx,my = grid.ng[0],grid.ng[1]
            maxm = max(mx,my)
            
            dx,dy,dt = grid.d[0],grid.d[1],self.dt
            
            self.apply_q_bcs(state)
            qnew = self.qbc #(input/output)
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
                      qold,qnew,self.auxbc,dx,dy,dt,self.method,self.mthlim,self.cfl,self.cflv, \
                      self.aux1,self.aux2,self.aux3,self.work)
            else:
                q, cfl = classic2.step2(maxm,mx,my,mbc,mx,my, \
                      qold,qnew,self.auxbc,dx,dy,dt,self.method,self.mthlim,self.cfl, \
                      self.aux1,self.aux2,self.aux3,self.work)

            self.cfl = cfl
            self.copy_local_to_global(q,state,self.mbc)

        else:
            raise NotImplementedError("No python implementation for homogeneous_step in case of 2D.")
