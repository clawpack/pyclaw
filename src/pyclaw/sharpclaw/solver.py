r"""
Module containing SharpClaw solvers for PyClaw/PetClaw

#  File:        sharpclaw.py
#  Created:     2010-03-20
#  Author:      David Ketcheson
"""
# Solver superclass
from clawpack.pyclaw.solver import Solver
from clawpack.pyclaw.util import add_parent_doc

# Reconstructor
try:
    # load c-based WENO reconstructor (PyWENO)
    from clawpack.pyclaw.limiters import reconstruct as recon
except ImportError:
    # load old WENO5 reconstructor
    from clawpack.pyclaw.limiters import recon

def default_tfluct(self):
    r"""This is a dummy routine and should never be called, check Euler1D
        to learn how to pass tfluct functions to the sharpclaw solver
    """
    if self.tfluct_solver:
        raise Exception("You set solver.tfluct_solver=True, but solver.tfluct has not been set.")
    pass

class SharpClawSolver(Solver):
    r"""
    Superclass for all SharpClawND solvers.

    Implements Runge-Kutta time stepping and the basic form of a 
    semi-discrete step (the dq() function).  If another method-of-lines
    solver is implemented in the future, it should be based on this class,
    which then ought to be renamed to something like "MOLSolver".

    .. attribute:: lim_type

        Limiter(s) to be used.
            - 0: No limiting.
            - 1: TVD reconstruction.
            - 2: WENO reconstruction.

        ``Default = 2``

    .. attribute:: weno_order

        Order of the WENO reconstruction. From 1st to 17th order (PyWENO)

        ``Default = 5``

    .. attribute:: time_integrator

        Time integrator to be used. Currently implemented methods:

            - 'Euler'  : 1st-order Forward Euler integration
            - 'SSP33'  : 3rd-order strong stability preserving method of Shu & Osher
            - 'SSP104' : 4th-order strong stability preserving method Ketcheson
            - 'SSPLMM32': 2nd-order strong stability preserving 3-step linear multistep method,
                          using Euler for starting values
            - 'SSPLMM43': 3rd-order strong stability preserving 4-step linear multistep method
                          using SSPRK22 for starting values
            - 'RK'     : Arbitrary Runge-Kutta method, specified by setting `solver.a`
                         and `solver.b` to the Butcher arrays of the method.
            - 'LMM'    : Arbitrary linear multistep method, specified by setting the
                         coefficient arrays `solver.alpha` and `solver.beta`.

        ``Default = 'SSP104'``

    .. attribute:: char_decomp

        Type of WENO reconstruction.
        0: conservative variables WENO reconstruction (standard).
        1: Wave-slope reconstruction.
        2: characteristic-wise WENO reconstruction.
        3: transmission-based WENO reconstruction.
        ``Default = 0``

    .. attribute:: tfluct_solver

        Whether a total fluctuation solver have to be used. If True the function
        that calculates the total fluctuation must be provided.
        ``Default = False``

    .. attribute:: tfluct

        Pointer to Fortran routine to calculate total fluctuation
        ``Default = default_tfluct (None)``

    .. attribute:: aux_time_dep

        Whether the auxiliary array is time dependent.
        ``Default = False``
    
    .. attribute:: kernel_language

        Specifies whether to use wrapped Fortran routines ('Fortran')
        or pure Python ('Python').  
        ``Default = 'Fortran'``.

    .. attribute:: num_ghost

        Number of ghost cells.
        ``Default = 3``

    .. attribute:: fwave
    
        Whether to split the flux jump (rather than the jump in Q) into waves; 
        requires that the Riemann solver performs the splitting.  
        ``Default = False``

    .. attribute:: cfl_desired

        Desired CFL number.
        ``Default = 2.45``

    .. attribute:: cfl_max

        Maximum CFL number.
        ``Default = 2.50``

    .. attribute:: dq_src

        Whether a source term is present. If it is present the function that 
        computes its contribution must be provided.
        ``Default = None``

    .. attribute:: call_before_step_each_stage

        Whether to call the method `self.before_step` before each RK stage.
        ``Default = False``

    """
    _sspcoeff = {
       'Euler' :    1.0,
       'SSP22':     1.0,
       'SSP33':     1.0,
       'SSP104':    6.0,
       'SSPLMM32':  0.5,
       'SSPLMM43':  1./3.,
       'SSPLMM53':  0.5,
       'RK':        None,
       'LMM':       None
       }

    _cfl_default = {
        'SSP104':   [2.45, 2.5],
        'SSPLMM32': [0.24, 0.25],
        'SSPLMM43': [0.15, 1./6.],
        'SSPLMM53': [0.24, 0.25]
        }

    # ========================================================================
    #   Initialization routines
    # ========================================================================
    def __init__(self,riemann_solver=None,claw_package=None):
        r"""
        Set default options for SharpClawSolvers and call the super's __init__().
        """
        self.limiters = [1]
        self.lim_type = 2
        self.weno_order = 5
        self.time_integrator = 'SSP104'
        self.char_decomp = 0
        self.tfluct_solver = False
        self.tfluct = default_tfluct
        self.aux_time_dep = False
        self.kernel_language = 'Fortran'
        self.num_ghost = 3
        self.fwave = False
        self.cfl_desired = None
        self.cfl_max = None
        self.dq_src = None
        self.call_before_step_each_stage = False
        self._mthlim = self.limiters
        self._method = None
        self._registers = None
        self.dq_dt = None
        self.dt_old = None        

        # Used only if time integrator is 'RK'
        self.a = None
        self.b = None
        self.c = None

        # Used only if time integrator is a multistep method
        self.sspcoeff0 = None
        self.alpha = None
        self.beta = None
        self.lmm_steps = 4
        self.sspcoeff = None
        self.prev_dq_dt_values = []
        self.prev_dt_values = []
        self.prev_dtFE_values = []
        self.check_lmm_cond = False
        self.lmm_cond = True

        # Call general initialization function
        super(SharpClawSolver,self).__init__(riemann_solver,claw_package)

 
    def setup(self,solution):
        """
        Allocate RK stage arrays or previous step solutions and fortran routine work arrays.
        """
        if self.lim_type == 2:
            self.num_ghost = (self.weno_order+1)//2

        if self.lim_type == 2 and self.weno_order != 5 and self.kernel_language == 'Python':
            raise Exception('Only 5th-order WENO reconstruction is implemented in Python kernels. \
                             Use Fortran for higher-order WENO.')
        # This is a hack to deal with the fact that petsc4py
        # doesn't allow us to change the stencil_width (num_ghost)
        state = solution.state
        state.set_num_ghost(self.num_ghost)
        # End hack

        if self.time_integrator == 'LMM':
            assert self.dt_variable == False, \
                'Must set solver.dt_variable=False for LMM integrator.'
        try:
            if 'SSPLMM' in self.time_integrator:
                if self.time_integrator == 'SSPLMMk2':
                    assert 3 <= self.lmm_steps, \
                        'Must set solver.lmm_steps greater than 2 for 2nd order SSPLMM integrator.'
                    self.sspcoeff = (self.lmm_steps - 2.)/(self.lmm_steps - 1.)
                    # The choice of cfl_desired and cfl_max is intended for LMM with many steps (up to 20). 
                    # If more steps are chosen the solution may not be accurate enough.
                    # For a smaller number of steps, higher values of cfl_desired and cfl_max can be used.
                    if self.cfl_max is None:
                        self.cfl_desired = 0.14*self.sspcoeff 
                        self.cfl_max = 0.15*self.sspcoeff 
                elif self.time_integrator == 'SSPLMMk3':
                    assert 4 <= self.lmm_steps <= 5, \
                        'Must set solver.lmm_steps equal to 4 or 5 for 3rd order SSPLMM integrator.'
                    self.sspcoeff = (self.lmm_steps - 3.)/(self.lmm_steps - 1.)
                    if self.cfl_max is None:
                        self.cfl_desired = 0.48*self.sspcoeff 
                        self.cfl_max = 0.5*self.sspcoeff 
            else:
                if self.cfl_max is None:
                    self.cfl_desired = self._cfl_default[self.time_integrator][0]
                    self.cfl_max = self._cfl_default[self.time_integrator][1]
            if self.cfl_desired is None:
                self.cfl_desired = 0.9*self.cfl_max
        except KeyError:
            raise KeyError('Maximum CFL number is not provided.')

        self._allocate_registers(solution)
        self._set_mthlim()

        state = solution.states[0]
 
        if self.kernel_language=='Fortran':
            if self.fmod is None:
                so_name = 'clawpack.pyclaw.sharpclaw.sharpclaw'+str(self.num_dim)
                self.fmod = __import__(so_name,fromlist=['clawpack.pyclaw.sharpclaw'])
            state.set_cparam(self.fmod)
            state.set_cparam(self.rp)
            state.set_cparam(self.tfluct)
            self._set_fortran_parameters(state,self.fmod.clawparams,self.fmod.workspace,self.fmod.reconstruct)

        self._allocate_bc_arrays(state)

        super(SharpClawSolver,self).setup(solution)


    def __del__(self):
        r"""
        Deallocate F90 module arrays.
        Also delete Fortran objects, which otherwise tend to persist in Python sessions.
        """
        if self.kernel_language=='Fortran':
            self.fmod.clawparams.dealloc_clawparams()
            self.fmod.workspace.dealloc_workspace(self.char_decomp)
            self.fmod.reconstruct.dealloc_recon_workspace(self.fmod.clawparams.lim_type,self.fmod.clawparams.char_decomp)
            del self.fmod

        super(SharpClawSolver,self).__del__()


    # ========== Time stepping routines ======================================
    def step(self,solution,take_one_step,tstart,tend):
        """Evolve q over one time step.

        Take one step with a Runge-Kutta or multistep method as specified by
        `solver.time_integrator`.
        """
        state = solution.states[0]
        step_index = self.status['numsteps'] + 1
        if self.accept_step == True:
            self.cfl.set_global_max(0.)
            self.dq_dt = self.dq(state) / self.dt

        if 'LMM' in self.time_integrator:
            step_index = self.update_saved_values(state,step_index)

        self.get_dt(solution.t,tstart,tend,take_one_step)

        # Recompute cfl number based on current step-size
        cfl = self.cfl.get_cached_max()
        self.cfl.set_global_max(self.dt / self.dt_old * cfl)
        self.dt_old = self.dt

        ### Runge-Kutta methods ###
        if self.time_integrator == 'Euler':
            state.q += self.dt*self.dq_dt

        elif self.time_integrator == 'SSP22':
            self.ssp22(state)

        elif self.time_integrator == 'SSP33':
            self._registers[0].q = state.q + self.dt*self.dq_dt
            self._registers[0].t = state.t + self.dt

            if self.call_before_step_each_stage:
                self.before_step(self,self._registers[0])
            self._registers[0].q = 0.75*state.q + 0.25*(self._registers[0].q + self.dq(self._registers[0]))
            self._registers[0].t = state.t + 0.5*self.dt

            if self.call_before_step_each_stage:
                self.before_step(self,self._registers[0])
            state.q = 1./3.*state.q + 2./3.*(self._registers[0].q + self.dq(self._registers[0]))

        elif self.time_integrator == 'SSP104':
            self.ssp104(state)

        elif self.time_integrator == 'RK':
            # General explicit RK with specified coefficients

            # This is pulled out of the loop in order to use dq_dt
            self._registers[0].q = self.dt*self.dq_dt
            self._registers[0].t = state.t

            num_stages = len(self.b)
            for i in range(1,num_stages):
                self._registers[i].q[:] = state.q
                self._registers[i].t = state.t + self.dt*self.c[i]
                if self.call_before_step_each_stage:
                    self.before_step(self,self._registers[i])
                for j in range(i):
                    self._registers[i].q += self.a[i,j]*self._registers[j].q

                # self._registers[i].q eventually stores dt*f(y_i) after stage solution y_i is computed
                self._registers[i].q = self.dq(self._registers[i])

            for j in range(num_stages):
                state.q += self.b[j]*self._registers[j].q

        ### Linear multistep methods ###
        elif self.time_integrator in ['SSPLMMk2', 'SSPLMMk3']:
            num_steps = self.lmm_steps
            if step_index < num_steps:
                # Use SSP22 Runge-Kutta method for starting values
                self.ssp22(state)
            else:
                if self.time_integrator == 'SSPLMMk2':
                    omega_k_minus_1 = sum(self.prev_dt_values[1:])/self.dt
                    r = (omega_k_minus_1-1.)/omega_k_minus_1 # SSP coefficient

                    delta = 1./omega_k_minus_1**2
                    beta = (omega_k_minus_1+1.)/omega_k_minus_1
                    state.q = beta*(r*state.q + self.dt*self.dq_dt) + delta*self._registers[-num_steps].q
                else:
                    omega_k_minus_1 = sum(self.prev_dt_values[1:])/self.dt
                    omega_k = omega_k_minus_1 + 1.
                    r = (omega_k_minus_1-2.)/omega_k_minus_1 # SSP coefficient

                    delta0 = (4*omega_k - omega_k_minus_1**2)/omega_k_minus_1**3
                    beta0 = omega_k/omega_k_minus_1**2
                    beta_k_minus_1 = omega_k**2/omega_k_minus_1**2

                    state.q = beta_k_minus_1*(r*state.q + self.dt*self.prev_dq_dt_values[-1]) + \
                            (r*beta0 + delta0)*self._registers[-num_steps].q + \
                            beta0*self.dt*self.prev_dq_dt_values[-num_steps]

        elif self.time_integrator == 'LMM':
            if step_index < len(self._registers):
                self.ssp104(state) # Use SSP104 for starting values
            else:
                # Update solution: alpha[-1] and beta[-1] correspond to solution at the previous step
                state.q = self.alpha[-1]*self._registers[-1].q + self.beta[-1]*self.dt*self.prev_dq_dt_values[-1]
                for i in range(self.lmm_steps-1):
                    state.q += self.alpha[i]*self._registers[i].q + self.beta[i]*self.dt*self.prev_dq_dt_values[i]

        else:
            raise Exception('Unrecognized time integrator')
            return False


    def ssp22(self,state):
        self._registers[0].q = state.q + self.dt*self.dq_dt
        self._registers[0].t = state.t + self.dt

        if self.call_before_step_each_stage:
            self.before_step(self,self._registers[0])
        state.q = 0.5*(state.q + self._registers[0].q + self.dq(self._registers[0]))


    def ssp104(self,state):
        if self.time_integrator == 'SSP104':
            s1 = self._registers[0]
            s1.q[:] = state.q
        elif self.time_integrator == 'LMM':
            # Okay to copy state objects here since this only happens a few times
            import copy
            s1 = copy.deepcopy(state)

        s1.q = state.q + self.dt*self.dq_dt/6.
        s1.t = state.t + self.dt/6.

        for i in range(4):
            if self.call_before_step_each_stage:
                self.before_step(self,s1)
            s1.q = s1.q + self.dq(s1)/6.
            s1.t = s1.t + self.dt/6.

        state.q = state.q/25. + 0.36 * s1.q
        s1.q = 15. * state.q - 5. * s1.q
        s1.t = state.t + self.dt/3.

        for i in range(4):
            if self.call_before_step_each_stage:
                self.before_step(self,s1)
            s1.q = s1.q + self.dq(s1)/6.
            s1.t = s1.t + self.dt/6.

        if self.call_before_step_each_stage:
            self.before_step(self,s1)
        state.q += 0.6 * s1.q + 0.1 * self.dq(s1)


    def update_saved_values(self,state,step_index):
        r""" 
        Updates lists of saved function evaluations, solution values, dt and dtFE for LMMs.
        For 3rd-order SSPLMM additional conditions are checked if self.check_lmm_cond is set to True.
        If these conditions are violated, the step is rejected.
        """
        if (self.prev_dt_values == []):
            # This only happens at the very beginning of the computation
            self._registers[-1].q[:] = state.q
            self._registers[-1].t = state.t
            self.prev_dq_dt_values.append(self.dq_dt)
            self.prev_dt_values.append(0.) # Not used
            if 'SSPLMM' in self.time_integrator:
                cfl = self.cfl.get_cached_max()
                dtFE = self.dt / cfl * self.cfl_max / self.sspcoeff
                self.prev_dtFE_values.append(dtFE)
                self.sspcoeff0 = self._sspcoeff['SSP22']
        elif self.accept_step == True: # Previous step was accepted
            if 'SSPLMM' in self.time_integrator:
                cfl = self.cfl.get_cached_max()
                dtFE = self.dt / cfl * self.cfl_max / self.sspcoeff

                if self.time_integrator == 'SSPLMMk3' and self.check_lmm_cond:
                    self.lmm_cond = self.check_3rd_ord_cond(state,step_index,dtFE)
                    if not self.lmm_cond:
                        self.accept_step = False
                        state.q[:] = self._registers[-1].q
                        self.dq_dt = self.prev_dq_dt_values[-1]
                        state.t = self._registers[-1].t
                        self.status['numsteps'] -= 1
                        return self.status['numsteps'] + 1

                if step_index <= len(self._registers):  # Startup
                    if self.time_integrator == 'SSPLMMk3':
                        self.prev_dq_dt_values.append(self.dq_dt)
                    self.prev_dt_values.append(self.dt_old)
                    self.prev_dtFE_values.append(dtFE)
                else:
                    if self.time_integrator == 'SSPLMMk3':
                        self.prev_dq_dt_values = self.prev_dq_dt_values[1:] + self.prev_dq_dt_values[:1]
                        self.prev_dq_dt_values[-1] = self.dq_dt

                    # Roll and update prev_dt_values and prev_dtFE_values lists
                    self.prev_dt_values = self.prev_dt_values[1:] + self.prev_dt_values[:1]
                    self.prev_dt_values[-1] = self.dt_old
                    self.prev_dtFE_values = self.prev_dtFE_values[1:] + self.prev_dtFE_values[:1]
                    self.prev_dtFE_values[-1] = dtFE
            elif self.time_integrator == 'LMM':
                if step_index <= len(self._registers):
                    self.prev_dq_dt_values.append(self.dq_dt)
                else:
                    self.prev_dq_dt_values = self.prev_dq_dt_values[1:] + self.prev_dq_dt_values[:1]
                    self.prev_dq_dt_values[-1] = self.dq_dt

            # Roll and update saved solution
            self._registers = self._registers[1:] + self._registers[:1]
            self._registers[-1].q[:] = state.q
            self._registers[-1].t = state.t

        return step_index


    def check_3rd_ord_cond(self,state,step_index,dtFE):
        r"""
        This routine checks the additional conditions for the 3rd-order SSPLMMs.
        This is a posteriori check after a step is accepted.
        In particular, there is a condition on the step size for the starting values and 
        a condition on the ratio of forward Euler step sizes at very step.
        If the conditions are violated we muct retrieve the previous solution and discard
        that step; otherwise the step is accepted.
        """
        lmm_cond = True

        if step_index <= len(self._registers):
            rho = 0.6 if len(self._registers) == 4 else 0.57
            if self.dt > rho * dtFE:
                lmm_cond = False

        rhoFE = 0.9 if len(self._registers) == 4 else 0.962
        dtFEm1 = self.prev_dtFE_values[-1]
        if rhoFE * dtFEm1 > dtFE or dtFE > dtFEm1 / rhoFE:
            lmm_cond = False

        return lmm_cond


    def _set_mthlim(self):
        self._mthlim = self.limiters
        if not isinstance(self.limiters,list): self._mthlim=[self._mthlim]
        if len(self._mthlim)==1: self._mthlim = self._mthlim * self.num_waves
        if len(self._mthlim)!=self.num_waves:
            raise Exception('Length of solver.limiters is not equal to 1 or to solver.num_waves')

       
    def dq(self,state):
        """
        Evaluate dq/dt * (delta t)
        """

        deltaq = self.dq_hyperbolic(state)

        if self.dq_src is not None:
            deltaq+=self.dq_src(self,state,self.dt)

        return deltaq


    def dq_hyperbolic(self,state):
        raise NotImplementedError('You must subclass SharpClawSolver.')

         
    def dqdt(self,state):
        """
        Evaluate dq/dt.  This routine is used for implicit time stepping.
        """

        self.dt = 1
        dq_dt = self.dq_hyperbolic(state)

        if self.dq_src is not None:
            dq_dt += self.dq_src(self,state,self.dt)

        return dq_dt.flatten('f')


    def _set_fortran_parameters(self,state,clawparams,workspace,reconstruct):
        """
        Set parameters for Fortran modules used by SharpClaw.
        The modules should be imported and passed as arguments to this function.
        """
        grid = state.grid
        clawparams.num_dim       = grid.num_dim
        clawparams.lim_type      = self.lim_type
        clawparams.weno_order    = self.weno_order
        clawparams.char_decomp   = self.char_decomp
        clawparams.tfluct_solver = self.tfluct_solver
        clawparams.fwave         = self.fwave
        clawparams.index_capa         = state.index_capa+1

        clawparams.num_waves     = self.num_waves
        clawparams.alloc_clawparams()
        for idim in range(grid.num_dim):
            clawparams.xlower[idim]=grid.dimensions[idim].lower
            clawparams.xupper[idim]=grid.dimensions[idim].upper
        clawparams.dx       =grid.delta
        clawparams.mthlim   =self._mthlim

        maxnx = max(grid.num_cells)+2*self.num_ghost
        workspace.alloc_workspace(maxnx,self.num_ghost,state.num_eqn,self.num_waves,self.char_decomp)
        reconstruct.alloc_recon_workspace(maxnx,self.num_ghost,state.num_eqn,self.num_waves,
                                            clawparams.lim_type,clawparams.char_decomp)


    def _allocate_registers(self,solution):
        r"""
        Instantiate State objects for Runge--Kutta stages and Linear Multistep method steps.

        This routine is only used by method-of-lines solvers (SharpClaw),
        not by the Classic solvers.  It allocates additional State objects
        to store the intermediate stages used by Runge--Kutta and Multistep 
        time integrators.

        If we create a MethodOfLinesSolver subclass, this should be moved there.
        """
        # Generally the number of registers for the starting method should be at most 
        # equal to the number of registers of the LMM
        if self.time_integrator   == 'Euler':   nregisters=0
        elif self.time_integrator == 'SSP33':   nregisters=1
        elif self.time_integrator == 'SSP104':  nregisters=1
        elif self.time_integrator == 'RK':      nregisters=len(self.b)+1
        elif self.time_integrator == 'SSPLMMk2': nregisters=self.lmm_steps
        elif self.time_integrator == 'SSPLMMk3': nregisters=self.lmm_steps
        elif self.time_integrator == 'LMM': nregisters=len(self.alpha)
        else:
            raise Exception('Unrecognized time integrator: '+self.time_integrator)
        
        state = solution.states[0]
        # Use the same class constructor as the solution for the Runge Kutta stages
        self._registers = []
        for i in range(nregisters):
            import copy
            self._registers.append(copy.deepcopy(state))


    def accept_reject_step(self,state):
        r"""
        Decide whether to accept or not the current step.
        For Runge-Kutta methods the step is accepted if cfl <= cfl_max.
        For SSPLMM32 the choice of step-size guarantees the cfl condition is satisfied for the steps the LMM
        is used. Hence, we need to check the cfl condition only for the starting steps. 
        """
        accept_step = True
        cfl = self.cfl.get_cached_max()

        if 'LMM' in self.time_integrator:
            step_index = self.status['numsteps'] + 1

            # Condition for starting RK methods
            if step_index < len(self._registers):
                if self.time_integrator == 'LMM':
                    sspcoeff_ratio = 1.
                else:
                    sspcoeff_ratio = self.sspcoeff0/self.sspcoeff

                if cfl > sspcoeff_ratio * self.cfl_max:
                    accept_step = False

        # Check cfl condition for Runge-Kutta methods
        else:
            if cfl > self.cfl_max:
                accept_step = False

        return accept_step


    def get_dt_new(self):
        r"""
        Set size of next step depending on the time integrator and
        whether or not the current step was accepted.
        """
        self.dt_old = self.dt        
        cfl = self.cfl.get_cached_max()

        if 'SSPLMM' in self.time_integrator:
            step_index = self.status['numsteps'] + 1

            if step_index < len(self._registers):
                # Step-size update of starting methods
                sspcoeff_ratio = self.sspcoeff0/self.sspcoeff
                self.dt = sspcoeff_ratio * self.dt * self.cfl_desired / cfl
                if self.time_integrator == 'SSPLMMk3' and self.check_lmm_cond and not self.lmm_cond:
                    rho = 0.6 if len(self._registers)== 4 else 0.57
                    self.dt = rho * self.dt
            else:
                # Step size selection guarantees CFL condition is satisfied.
                # Only need to check 3rd-order LMM's condition
                if self.accept_step:
                    s = len(self._registers)
                    p = int(self.time_integrator[-1])
                    mu = min([self.prev_dtFE_values[i] for i in range(s)])
                    H = sum(self.prev_dt_values[1:])
                    self.dt = H * mu / (H + (p-1)*mu)
                elif self.time_integrator == 'SSPLMMk3' and self.check_lmm_cond:
                    self.dt = 0.5 * self.dt

        # Step-size update for RK methods
        else:
            self.dt = self.dt * self.cfl_desired / cfl



# ========================================================================
class SharpClawSolver1D(SharpClawSolver):
# ========================================================================
    """
    SharpClaw solver for one-dimensional problems.
    
    Used to solve 1D hyperbolic systems using the SharpClaw algorithms,
    which are based on WENO reconstruction and Runge-Kutta time stepping.
    """

    __doc__ += add_parent_doc(SharpClawSolver)
    
    def __init__(self,riemann_solver=None,claw_package=None):
        r"""
        See :class:`SharpClawSolver1D` for more info.
        """   
        self.num_dim = 1
        self.reflect_index = [1]
        super(SharpClawSolver1D,self).__init__(riemann_solver,claw_package)


    def dq_hyperbolic(self,state):
        r"""
        Compute dq/dt * (delta t) for the hyperbolic hyperbolic system.

        Note that the capa array, if present, should be located in the aux
        variable.

        Indexing works like this (here num_ghost=2 as an example)::

         0     1     2     3     4     mx+num_ghost-2     mx+num_ghost      mx+num_ghost+2
                     |                        mx+num_ghost-1 |  mx+num_ghost+1
         |     |     |     |     |   ...   |     |     |     |     |
            0     1  |  2     3            mx+num_ghost-2    |mx+num_ghost       
                                                  mx+num_ghost-1   mx+num_ghost+1

        The top indices represent the values that are located on the grid
        cell boundaries such as waves, s and other Riemann problem values, 
        the bottom for the cell centered values such as q.  In particular
        the ith grid cell boundary has the following related information::

                          i-1         i         i+1
                           |          |          |
                           |   i-1    |     i    |
                           |          |          |

        Again, grid cell boundary quantities are at the top, cell centered
        values are in the cell.

        """
    
        import numpy as np

        self._apply_bcs(state)
        q = self.qbc
        aux = self.auxbc 

        grid = state.grid
        mx = grid.num_cells[0]

        ixy=1

        if self.kernel_language=='Fortran':
            rp1 = self.rp.rp1._cpointer
            if self.tfluct_solver:
                tfluct1 = self.tfluct.tfluct1._cpointer
            else:
                tfluct1 = self.tfluct

            dq,cfl=self.fmod.flux1(q,self.auxbc,self.dt,state.t,ixy,mx,self.num_ghost,mx,rp1,tfluct1)

        elif self.kernel_language=='Python':

            dtdx = np.zeros( (mx+2*self.num_ghost) ,order='F')
            dq   = np.zeros( (state.num_eqn,mx+2*self.num_ghost) ,order='F')

            # Find local value for dt/dx
            if state.index_capa>=0:
                dtdx = self.dt / (grid.delta[0] * aux[state.index_capa,:])
            else:
                dtdx += self.dt/grid.delta[0]
 
            if aux.shape[0]>0:
                aux_l=aux[:,:-1]
                aux_r=aux[:,1: ]
            else:
                aux_l = None
                aux_r = None

            #Reconstruct (wave reconstruction uses a Riemann solve)
            if self.lim_type==-1: #1st-order Godunov
                ql=q; qr=q
            elif self.lim_type==0: #Unlimited reconstruction
                raise NotImplementedError('Unlimited reconstruction not implemented')
            elif self.lim_type==1: #TVD Reconstruction
                raise NotImplementedError('TVD reconstruction not implemented')
            elif self.lim_type==2: #WENO Reconstruction
                if self.char_decomp==0: #No characteristic decomposition
                    ql,qr=recon.weno(5,q)
                elif self.char_decomp==1: #Wave-based reconstruction
                    q_l=q[:,:-1]
                    q_r=q[:,1: ]
                    wave,s,amdq,apdq = self.rp(q_l,q_r,aux_l,aux_r,state.problem_data)
                    ql,qr=recon.weno5_wave(q,wave,s)
                elif self.char_decomp==2: #Characteristic-wise reconstruction
                    raise NotImplementedError

            # Solve Riemann problem at each interface
            q_l=qr[:,:-1]
            q_r=ql[:,1: ]
            wave,s,amdq,apdq = self.rp(q_l,q_r,aux_l,aux_r,state.problem_data)

            # Loop limits for local portion of grid
            # THIS WON'T WORK IN PARALLEL!
            LL = self.num_ghost - 1
            UL = grid.num_cells[0] + self.num_ghost + 1

            # Compute maximum wave speed
            cfl = 0.0
            for mw in range(self.num_waves):
                smax1 = np.max( dtdx[LL  :UL]  *s[mw,LL-1:UL-1])
                smax2 = np.max(-dtdx[LL-1:UL-1]*s[mw,LL-1:UL-1])
                cfl = max(cfl,smax1,smax2)

            #Find total fluctuation within each cell
            wave,s,amdq2,apdq2 = self.rp(ql,qr,aux,aux,state.problem_data)

            # Compute dq
            for m in range(state.num_eqn):
                dq[m,LL:UL] = -dtdx[LL:UL]*(amdq[m,LL:UL] + apdq[m,LL-1:UL-1] \
                                + apdq2[m,LL:UL] + amdq2[m,LL:UL])

        else: 
            raise Exception('Unrecognized value of solver.kernel_language.')

        self.cfl.update_global_max(cfl)
        return dq[:,self.num_ghost:-self.num_ghost]
    

# ========================================================================
class SharpClawSolver2D(SharpClawSolver):
# ========================================================================
    """ Two Dimensional SharpClawSolver
    """

    __doc__ += add_parent_doc(SharpClawSolver)

    def __init__(self,riemann_solver=None,claw_package=None):
        r"""
        Create 2D SharpClaw solver
        
        See :class:`SharpClawSolver2D` for more info.
        """   
        self.num_dim = 2
        self.reflect_index = [1,2]

        super(SharpClawSolver2D,self).__init__(riemann_solver,claw_package)


    def dq_hyperbolic(self,state):
        """Compute dq/dt * (delta t) for the hyperbolic hyperbolic system

        Note that the capa array, if present, should be located in the aux
        variable.

        Indexing works like this (here num_ghost=2 as an example)::

         0     1     2     3     4     mx+num_ghost-2     mx+num_ghost      mx+num_ghost+2
                     |                        mx+num_ghost-1 |  mx+num_ghost+1
         |     |     |     |     |   ...   |     |     |     |     |
            0     1  |  2     3            mx+num_ghost-2    |mx+num_ghost       
                                                  mx+num_ghost-1   mx+num_ghost+1

        The top indices represent the values that are located on the grid
        cell boundaries such as waves, s and other Riemann problem values, 
        the bottom for the cell centered values such as q.  In particular
        the ith grid cell boundary has the following related information::

                          i-1         i         i+1
                           |          |          |
                           |   i-1    |     i    |
                           |          |          |

        Again, grid cell boundary quantities are at the top, cell centered
        values are in the cell.

        """
        self._apply_bcs(state)
        q = self.qbc 

        grid = state.grid

        num_ghost=self.num_ghost
        mx=grid.num_cells[0]
        my=grid.num_cells[1]
        maxm = max(mx,my)

        if self.kernel_language=='Fortran':
            rpn2 = self.rp.rpn2._cpointer
            if self.tfluct_solver:
                tfluct2 = self.tfluct.tfluct2._cpointer
            else:
                tfluct2 = self.tfluct

            dq,cfl=self.fmod.flux2(q,self.auxbc,self.dt,state.t,num_ghost,maxm,mx,my,rpn2,tfluct2)

        else: raise Exception('Only Fortran kernels are supported in 2D.')

        self.cfl.update_global_max(cfl)
        return dq[:,num_ghost:-num_ghost,num_ghost:-num_ghost]

# ========================================================================
class SharpClawSolver3D(SharpClawSolver):
# ========================================================================
    """ Three Dimensional SharpClawSolver
    """

    __doc__ += add_parent_doc(SharpClawSolver)

    def __init__(self,riemann_solver=None,claw_package=None):
        r"""
        Create 3D SharpClaw solver
        
        See :class:`SharpClawSolver3D` for more info.
        """   
        self.num_dim = 3
        self.reflect_index = [1,2,3]

        super(SharpClawSolver3D,self).__init__(riemann_solver,claw_package)


    def teardown(self):
        r"""
        Deallocate F90 module arrays.
        Also delete Fortran objects, which otherwise tend to persist in Python sessions.
        """
        if self.kernel_language=='Fortran':
            self.fmod.workspace.dealloc_workspace(self.char_decomp)
            self.fmod.reconstruct.dealloc_recon_workspace(self.fmod.clawparams.lim_type,self.fmod.clawparams.char_decomp)
            self.fmod.clawparams.dealloc_clawparams()
            del self.fmod


    def dq_hyperbolic(self,state):
        """Compute dq/dt * (delta t) for the hyperbolic hyperbolic system

        Note that the capa array, if present, should be located in the aux
        variable.

        Indexing works like this (here num_ghost=2 as an example)::

         0     1     2     3     4     mx+num_ghost-2     mx+num_ghost      mx+num_ghost+2
                     |                        mx+num_ghost-1 |  mx+num_ghost+1
         |     |     |     |     |   ...   |     |     |     |     |
            0     1  |  2     3            mx+num_ghost-2    |mx+num_ghost       
                                                  mx+num_ghost-1   mx+num_ghost+1

        The top indices represent the values that are located on the grid
        cell boundaries such as waves, s and other Riemann problem values, 
        the bottom for the cell centered values such as q.  In particular
        the ith grid cell boundary has the following related information::

                          i-1         i         i+1
                           |          |          |
                           |   i-1    |     i    |
                           |          |          |

        Again, grid cell boundary quantities are at the top, cell centered
        values are in the cell.

        """
        self._apply_bcs(state)
        q = self.qbc 

        grid = state.grid

        num_ghost=self.num_ghost
        mx=grid.num_cells[0]
        my=grid.num_cells[1]
        mz=grid.num_cells[2]
        maxm = max(mx,my,mz)

        if self.kernel_language=='Fortran':
            rpn3 = self.rp.rpn3._cpointer
            if self.tfluct_solver:            
                tfluct3 = self.tfluct.tfluct3._cpointer
            else:
                tfluct3 = self.tfluct

            dq,cfl=self.fmod.flux3(q,self.auxbc,self.dt,state.t,num_ghost,maxm,mx,my,mz,rpn3,tfluct3)

        else: raise Exception('Only Fortran kernels are supported in 3D.')

        self.cfl.update_global_max(cfl)
        return dq[:,num_ghost:-num_ghost,num_ghost:-num_ghost,num_ghost:-num_ghost]
