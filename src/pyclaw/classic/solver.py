r"""
Module containing the classic Clawpack solvers.

This module contains the pure and wrapped classic clawpack solvers.  All
clawpack solvers inherit from the :class:`ClawSolver` superclass which in turn
inherits from the :class:`~pyclaw.solver.Solver` superclass.  These
are both pure virtual classes; the only solver classes that should be instantiated
are the dimension-specific ones, :class:`ClawSolver1D` and :class:`ClawSolver2D`.
"""

from ..solver import Solver

from ..limiters import tvd

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

    .. attribute:: source_split

        Which source splitting method to use: 1 for first
        order Godunov splitting and 2 for second order Strang splitting.
        ``Default = 1``

    .. attribute:: fwave

        Whether to split the flux jump (rather than the jump in Q) into waves;
        requires that the Riemann solver performs the splitting.
        ``Default = False``

    .. attribute:: step_source

        Handle for function that evaluates the source term.
        The required signature for this function is:

        def step_source(solver,state,dt)

    .. attribute:: before_step

        Function called before each time step is taken.
        The required signature for this function is:

        def before_step(solver,solution)

    .. attribute:: kernel_language

        Specifies whether to use wrapped Fortran routines ('Fortran')
        or pure Python ('Python').  ``Default = 'Fortran'``.

    .. attribute:: verbosity

        The level of detail of logged messages from the Fortran solver.
        ``Default = 0``.

    Output:
     - (:class:`ClawSolver`) - Initialized clawpack solver
    """
    # ========== Generic Init Routine ========================================
    def __init__(self,riemann_solver=None,claw_package=None):
        r"""
        See :class:`ClawSolver` for full documentation.
        """
        self.num_ghost = 2
        self.limiters = tvd.minmod
        self.order = 2
        self.source_split = 1
        self.fwave = False
        self.step_source = None
        self.before_step = None
        self.kernel_language = 'Fortran'
        self.verbosity = 0
        self.cfl_max = 1.0
        self.cfl_desired = 0.9
        self._mthlim = self.limiters
        self._method = None

        # Call general initialization function
        super(ClawSolver,self).__init__(riemann_solver,claw_package)

    # ========== Time stepping routines ======================================
    def step(self,solution):
        r"""
        Evolve solution one time step

        The elements of the algorithm for taking one step are:

        1. The :meth:`before_step` function is called

        2. A half step on the source term :func:`step_source` if Strang splitting is
           being used (:attr:`source_split` = 2)

        3. A step on the homogeneous problem :math:`q_t + f(q)_x = 0` is taken

        4. A second half step or a full step is taken on the source term
           :func:`step_source` depending on whether Strang splitting was used
           (:attr:`source_split` = 2) or Godunov splitting
           (:attr:`source_split` = 1)

        This routine is called from the method evolve_to_time defined in the
        pyclaw.solver.Solver superclass.

        :Input:
         - *solution* - (:class:`~pyclaw.solution.Solution`) solution to be evolved

        :Output:
         - (bool) - True if full step succeeded, False otherwise
        """

        if self.before_step is not None:
            self.before_step(self,solution.states[0])

        if self.source_split == 2 and self.step_source is not None:
            self.step_source(self,solution.states[0],self.dt/2.0)

        self.step_hyperbolic(solution)

        # Check here if the CFL condition is satisfied.
        # If not, return # immediately to evolve_to_time and let it deal with
        # picking a new step size (dt).
        if self.cfl.get_cached_max() >= self.cfl_max:
            return False

        if self.step_source is not None:
            # Strang splitting
            if self.source_split == 2:
                self.step_source(self,solution.states[0],self.dt/2.0)

            # Godunov Splitting
            if self.source_split == 1:
                self.step_source(self,solution.states[0],self.dt)

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
        Convenience routine to convert users limiter specification to
        the format understood by the Fortran code (i.e., a list of length num_waves).
        """
        self._mthlim = self.limiters
        if not isinstance(self.limiters,list): self._mthlim=[self._mthlim]
        if len(self._mthlim)==1: self._mthlim = self._mthlim * self.num_waves
        if len(self._mthlim)!=self.num_waves:
            raise Exception('Length of solver.limiters is not equal to 1 or to solver.num_waves')

    def set_method(self,state):
        r"""
        Set values of the solver._method array required by the Fortran code.
        These are algorithmic parameters.
        """
        import numpy as np
        #We ought to put method and many other things in a Fortran
        #module and set the fortran variables directly here.
        self._method =np.empty(7, dtype=int,order='F')
        self._method[0] = self.dt_variable
        self._method[1] = self.order
        if self.num_dim==1:
            self._method[2] = 0  # Not used in 1D
        elif self.dimensional_split:
            self._method[2] = -1  # First-order dimensional splitting
        else:
            self._method[2] = self.transverse_waves
        self._method[3] = self.verbosity
        self._method[4] = 0  # Not used for PyClaw (would be self.source_split)
        self._method[5] = state.index_capa + 1
        self._method[6] = state.num_aux

    def setup(self,solution):
        r"""
        Perform essential solver setup.  This routine must be called before
        solver.step() may be called.
        """
        # This is a hack to deal with the fact that petsc4py
        # doesn't allow us to change the stencil_width (num_ghost)
        solution.state.set_num_ghost(self.num_ghost)
        # End hack

        self.check_cfl_settings()

        self.set_mthlim()
        if(self.kernel_language == 'Fortran'):
            if self.fmod is None:
                so_name = 'clawpack.pyclaw.classic.classic'+str(self.num_dim)
                self.fmod = __import__(so_name,fromlist=['clawpack.pyclaw.classic'])
            self.set_fortran_parameters(solution)
            self.allocate_workspace(solution)
        elif self.num_dim>1:
            raise Exception('Only Fortran kernels are supported in multi-D.')

        self.allocate_bc_arrays(solution.states[0])

        self._is_set_up = True


    def set_fortran_parameters(self,solution):
        r"""
        Pack parameters into format recognized by Clawpack (Fortran) code.

        Sets the solver._method array and the cparam common block for the Riemann solver.
        """
        self.set_method(solution.state)
        # The reload here is necessary because otherwise the common block
        # cparam in the Riemann solver doesn't get flushed between running
        # different tests in a single Python session.
        reload(self.fmod)
        solution.state.set_cparam(self.fmod)
        solution.state.set_cparam(self.rp)

    def teardown(self):
        r"""
        Delete Fortran objects, which otherwise tend to persist in Python sessions.
        """
        if(self.kernel_language == 'Fortran'):
            del self.fmod



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

    Output:
     - (:class:`ClawSolver1D`) - Initialized 1d clawpack solver
    """

    def __init__(self,riemann_solver=None,claw_package=None):
        r"""
        Create 1d Clawpack solver

        See :class:`ClawSolver1D` for more info.
        """
        self.num_dim = 1

        super(ClawSolver1D,self).__init__(riemann_solver,claw_package)


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
        if state.num_aux > 0:
            self.apply_aux_bcs(state)

        num_eqn,num_ghost = state.num_eqn,self.num_ghost

        if(self.kernel_language == 'Fortran'):
            mx = grid.num_cells[0]
            dx,dt = grid.delta[0],self.dt
            dtdx = np.zeros( (mx+2*num_ghost) ) + dt/dx
            rp1 = self.rp.rp1._cpointer

            self.qbc,cfl = self.fmod.step1(num_ghost,mx,self.qbc,self.auxbc,dx,dt,self._method,self._mthlim,self.fwave,rp1)

        elif(self.kernel_language == 'Python'):

            q   = self.qbc
            aux = self.auxbc
            # Limiter to use in the pth family
            limiter = np.array(self._mthlim,ndmin=1)

            dtdx = np.zeros( (2*self.num_ghost+grid.num_cells[0]) )

            # Find local value for dt/dx
            if state.index_capa>=0:
                dtdx = self.dt / (grid.delta[0] * state.aux[state.index_capa,:])
            else:
                dtdx += self.dt/grid.delta[0]

            # Solve Riemann problem at each interface
            q_l=q[:,:-1]
            q_r=q[:,1:]
            if state.aux is not None:
                aux_l=aux[:,:-1]
                aux_r=aux[:,1:]
            else:
                aux_l = None
                aux_r = None
            wave,s,amdq,apdq = self.rp(q_l,q_r,aux_l,aux_r,state.problem_data)

            # Update loop limits, these are the limits for the Riemann solver
            # locations, which then update a grid cell value
            # We include the Riemann problem just outside of the grid so we can
            # do proper limiting at the grid edges
            #        LL    |                               |     UL
            #  |  LL |     |     |     |  ...  |     |     |  UL  |     |
            #              |                               |

            LL = self.num_ghost - 1
            UL = self.num_ghost + grid.num_cells[0] + 1

            # Update q for Godunov update
            for m in xrange(num_eqn):
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
                f = np.zeros( (num_eqn,grid.num_cells[0] + 2*self.num_ghost) )

                # Apply Limiters to waves
                if (limiter > 0).any():
                    wave = tvd.limit(state.num_eqn,wave,s,limiter,dtdx)

                # Compute correction fluxes for second order q_{xx} terms
                dtdxave = 0.5 * (dtdx[LL-1:UL-1] + dtdx[LL:UL])
                if self.fwave:
                    for mw in xrange(wave.shape[1]):
                        sabs = np.abs(s[mw,LL-1:UL-1])
                        om = 1.0 - sabs*dtdxave[:UL-LL]
                        ssign = np.sign(s[mw,LL-1:UL-1])
                        for m in xrange(num_eqn):
                            f[m,LL:UL] += 0.5 * ssign * om * wave[m,mw,LL-1:UL-1]
                else:
                    for mw in xrange(wave.shape[1]):
                        sabs = np.abs(s[mw,LL-1:UL-1])
                        om = 1.0 - sabs*dtdxave[:UL-LL]
                        for m in xrange(num_eqn):
                            f[m,LL:UL] += 0.5 * sabs * om * wave[m,mw,LL-1:UL-1]

                # Update q by differencing correction fluxes
                for m in xrange(num_eqn):
                    q[m,LL:UL-1] -= dtdx[LL:UL-1] * (f[m,LL+1:UL] - f[m,LL:UL-1])

        else: raise Exception("Unrecognized kernel_language; choose 'Fortran' or 'Python'")

        self.cfl.update_global_max(cfl)
        state.set_q_from_qbc(num_ghost,self.qbc)
        if state.num_aux > 0:
            state.set_aux_from_auxbc(num_ghost,self.auxbc)


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

    .. attribute:: dimensional_split

        If True, use dimensional splitting (Godunov splitting).
        Dimensional splitting with Strang splitting is not supported
        at present but could easily be enabled if necessary.
        If False, use unsplit Clawpack algorithms, possibly including
        transverse Riemann solves.

    .. attribute:: transverse_waves

        If dimensional_split is True, this option has no effect.  If
        dim_plit is False, then transverse_waves should be one of
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

    def __init__(self,riemann_solver=None):
        r"""
        Create 2d Clawpack solver

        See :class:`ClawSolver2D` for more info.
        """
        self.dimensional_split = True
        self.transverse_waves = self.trans_inc

        self.num_dim = 2

        self.aux1 = None
        self.aux2 = None
        self.aux3 = None
        self.work = None

        super(ClawSolver2D,self).__init__(riemann_solver)

    def check_cfl_settings(self):
        if (not self.dimensional_split) and (self.transverse_waves==0):
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

        num_eqn,num_aux,num_waves,num_ghost,aux = state.num_eqn,state.num_aux,self.num_waves,self.num_ghost,state.aux

        #The following is a hack to work around an issue
        #with f2py.  It involves wastefully allocating three arrays.
        #f2py seems not able to handle multiple zero-size arrays being passed.
        # it appears the bug is related to f2py/src/fortranobject.c line 841.
        if(aux == None): num_aux=1

        grid  = state.grid
        maxmx,maxmy = grid.num_cells[0],grid.num_cells[1]
        maxm = max(maxmx, maxmy)

        # These work arrays really ought to live inside a fortran module
        # as is done for sharpclaw
        self.aux1 = np.empty((num_aux,maxm+2*num_ghost),order='F')
        self.aux2 = np.empty((num_aux,maxm+2*num_ghost),order='F')
        self.aux3 = np.empty((num_aux,maxm+2*num_ghost),order='F')
        mwork = (maxm+2*num_ghost) * (5*num_eqn + num_waves + num_eqn*num_waves)
        self.work = np.empty((mwork),order='F')


    # ========== Hyperbolic Step =====================================
    def step_hyperbolic(self,solution):
        r"""
        Take a step on the homogeneous hyperbolic system using the Clawpack
        algorithm.

        Clawpack is based on the Lax-Wendroff method, combined with Riemann
        solvers and TVD limiters applied to waves.
        """
        if(self.kernel_language == 'Fortran'):
            state = solution.states[0]
            grid = state.grid
            dx,dy = grid.delta
            mx,my = grid.num_cells
            maxm = max(mx,my)

            self.apply_q_bcs(state)
            if state.num_aux > 0:
                self.apply_aux_bcs(state)
            qold = self.qbc.copy('F')

            rpn2 = self.rp.rpn2._cpointer
            rpt2 = self.rp.rpt2._cpointer

            if self.dimensional_split:
                #Right now only Godunov-dimensional-splitting is implemented.
                #Strang-dimensional-splitting could be added following dimsp2.f in Clawpack.

                self.qbc, cfl_x = self.fmod.step2ds(maxm,self.num_ghost,mx,my, \
                      qold,self.qbc,self.auxbc,dx,dy,self.dt,self._method,self._mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work,1,self.fwave,rpn2,rpt2)

                self.qbc, cfl_y = self.fmod.step2ds(maxm,self.num_ghost,mx,my, \
                      self.qbc,self.qbc,self.auxbc,dx,dy,self.dt,self._method,self._mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work,2,self.fwave,rpn2,rpt2)

                cfl = max(cfl_x,cfl_y)

            else:

                self.qbc, cfl = self.fmod.step2(maxm,self.num_ghost,mx,my, \
                      qold,self.qbc,self.auxbc,dx,dy,self.dt,self._method,self._mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work,self.fwave,rpn2,rpt2)

            self.cfl.update_global_max(cfl)
            state.set_q_from_qbc(self.num_ghost,self.qbc)
            if state.num_aux > 0:
                state.set_aux_from_auxbc(self.num_ghost,self.auxbc)

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

    .. attribute:: dimensional_split

        If True, use dimensional splitting (Godunov splitting).
        Dimensional splitting with Strang splitting is not supported
        at present but could easily be enabled if necessary.
        If False, use unsplit Clawpack algorithms, possibly including
        transverse Riemann solves.

    .. attribute:: transverse_waves

        If dimensional_split is True, this option has no effect.  If
        dim_plit is False, then transverse_waves should be one of
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

    def __init__(self,riemann_solver=None):
        r"""
        Create 3d Clawpack solver

        See :class:`ClawSolver3D` for more info.
        """
        # Add the functions as required attributes
        self.dimensional_split = True
        self.transverse_waves = self.trans_cor

        self.num_dim = 3

        self.aux1 = None
        self.aux2 = None
        self.aux3 = None
        self.work = None

        super(ClawSolver3D,self).__init__(riemann_solver)

    # ========== Setup routine =============================
    def allocate_workspace(self,solution):
        r"""
        Allocate auxN and work arrays for use in Fortran subroutines.
        """
        import numpy as np

        state = solution.states[0]

        num_eqn,num_aux,num_waves,num_ghost,aux = state.num_eqn,state.num_aux,self.num_waves,self.num_ghost,state.aux

        #The following is a hack to work around an issue
        #with f2py.  It involves wastefully allocating three arrays.
        #f2py seems not able to handle multiple zero-size arrays being passed.
        # it appears the bug is related to f2py/src/fortranobject.c line 841.
        if(aux == None): num_aux=1

        grid  = state.grid
        maxmx,maxmy = grid.num_cells[0],grid.num_cells[1]
        maxm = max(maxmx, maxmy)

        # These work arrays really ought to live inside a fortran module
        # as is done for sharpclaw
        self.aux1 = np.empty((num_aux,maxm+2*num_ghost,3),order='F')
        self.aux2 = np.empty((num_aux,maxm+2*num_ghost,3),order='F')
        self.aux3 = np.empty((num_aux,maxm+2*num_ghost,3),order='F')
        mwork = (maxm+2*num_ghost) * (31*num_eqn + num_waves + num_eqn*num_waves)
        self.work = np.empty((mwork),order='F')


    # ========== Hyperbolic Step =====================================
    def step_hyperbolic(self,solution):
        r"""
        Take a step on the homogeneous hyperbolic system using the Clawpack
        algorithm.

        Clawpack is based on the Lax-Wendroff method, combined with Riemann
        solvers and TVD limiters applied to waves.
        """
        if(self.kernel_language == 'Fortran'):
            state = solution.states[0]
            grid = state.grid
            dx,dy,dz = grid.delta
            mx,my,mz = grid.num_cells
            maxm = max(mx,my,mz)

            self.apply_q_bcs(state)
            if state.num_aux > 0:
                self.apply_aux_bcs(state)
            qnew = self.qbc
            qold = qnew.copy('F')

            rpn3  = self.rp.rpn3._cpointer
            rpt3  = self.rp.rpt3._cpointer
            rptt3 = self.rp.rptt3._cpointer

            if self.dimensional_split:
                #Right now only Godunov-dimensional-splitting is implemented.
                #Strang-dimensional-splitting could be added following dimsp2.f in Clawpack.

                q, cfl_x = self.fmod.step3ds(maxm,self.num_ghost,mx,my,mz, \
                      qold,qnew,self.auxbc,dx,dy,dz,self.dt,self._method,self._mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work,1,rpn3,rpt3,rptt3)

                q, cfl_y = self.fmod.step3ds(maxm,self.num_ghost,mx,my,mz, \
                      q,q,self.auxbc,dx,dy,dz,self.dt,self._method,self._mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work,2,rpn3,rpt3,rptt3)

                q, cfl_z = self.fmod.step3ds(maxm,self.num_ghost,mx,my,mz, \
                      q,q,self.auxbc,dx,dy,dz,self.dt,self._method,self._mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work,3,rpn3,rpt3,rptt3)

                cfl = max(cfl_x,cfl_y,cfl_z)

            else:

                q, cfl = self.fmod.step3(maxm,self.num_ghost,mx,my,mz, \
                      qold,qnew,self.auxbc,dx,dy,dz,self.dt,self._method,self._mthlim,\
                      self.aux1,self.aux2,self.aux3,self.work,rpn3,rpt3,rptt3)

            self.cfl.update_global_max(cfl)
            state.set_q_from_qbc(self.num_ghost,self.qbc)
            if state.num_aux > 0:
                state.set_aux_from_auxbc(self.num_ghost,self.auxbc)

        else:
            raise NotImplementedError("No python implementation for step_hyperbolic in 3D.")
