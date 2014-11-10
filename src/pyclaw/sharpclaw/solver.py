r"""
Module containing SharpClaw solvers for PyClaw/PetClaw

#  File:        sharpclaw.py
#  Created:     2010-03-20
#  Author:      David Ketcheson
"""
# Solver superclass
from clawpack.pyclaw.solver import Solver, CFLError
from clawpack.pyclaw.util import add_parent_doc

# Reconstructor
try:
    # load c-based WENO reconstructor (PyWENO)
    from clawpack.pyclaw.limiters import reconstruct as recon
except ImportError:
    # load old WENO5 reconstructor
    from clawpack.pyclaw.limiters import recon

def before_step(solver,solution):
    r"""
    Dummy routine called before each step
    
    Replace this routine if you want to do something before each time step.
    """
    pass

def default_tfluct():
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

    .. attribute:: before_step
    
        Function called before each time step is taken.
        The required signature for this function is:
        
        def before_step(solver,solution)

    .. attribute:: lim_type

        Limiter(s) to be used.
        0: No limiting.
        1: TVD reconstruction.
        2: WENO reconstruction.
        ``Default = 2``

    .. attribute:: weno_order

        Order of the WENO reconstruction. From 1st to 17th order (PyWENO)
        ``Default = 5``

    .. attribute:: time_integrator

        Time integrator to be used. Currently implemented methods:

        'Euler'  : 1st-order Forward Euler integration
        'SSP33'  : 3rd-order strong stability preserving method of Shu & Osher
        'SSP104' : 4th-order strong stability preserving method Ketcheson
        'SSPMS32': 2nd-order strong stability preserving 3-step linear multistep method,
                   using Euler for starting values
        'SSPMS43': 3rd-order strong stability preserving 4-step linear multistep method
                   using SSPRK22 for starting values
        'RK'     : Arbitrary Runge-Kutta method, specified by setting `solver.a`
                   and `solver.b` to the Butcher arrays of the method.
        'LMM'    : Arbitrary linear multistep method, specified by setting the
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
       'SSP33':     1.0,
       'SSP104' :   6.0,
       'SSPMS32' :  0.5,
       'SSPMS43' :  1./3.,
       'RK':        None,
       'LMM':       None
       }

    _cfl_default = {
        'SSP104':   [2.45, 2.5],
        'SSPMS32':  [0.16, 0.2],
        'SSPMS43':  [0.14, 0.16]
        }

    # ========================================================================
    #   Initialization routines
    # ========================================================================
    def __init__(self,riemann_solver=None,claw_package=None):
        r"""
        Set default options for SharpClawSolvers and call the super's __init__().
        """
        self.limiters = [1]
        self.before_step = before_step
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

        # Used only if time integrator is 'RK'
        self.a = None
        self.b = None
        self.c = None

        # Used only if time integrator is a multistep method
        self.step_index = 1
        self.alpha = None
        self.beta = None

        # Call general initialization function
        super(SharpClawSolver,self).__init__(riemann_solver,claw_package)
        
    def setup(self,solution):
        """
        Allocate RK stage arrays or previous step solutions and fortran routine work arrays.
        """
        if self.lim_type == 2:
            self.num_ghost = (self.weno_order+1)/2

        if self.lim_type == 2 and self.weno_order != 5 and self.kernel_language == 'Python':
            raise Exception('Only 5th-order WENO reconstruction is implemented in Python kernels. \
                             Use Fortran for higher-order WENO.')
        # This is a hack to deal with the fact that petsc4py
        # doesn't allow us to change the stencil_width (num_ghost)
        state = solution.state
        state.set_num_ghost(self.num_ghost)
        # End hack

        self._allocate_registers(solution)
        self._set_mthlim()
        try:
            if self.cfl_max is None:
                self.cfl_desired  = self._cfl_default[self.time_integrator][0]
                self.cfl_max  = self._cfl_default[self.time_integrator][1]
            if self.cfl_desired is None:
                self.cfl_desired = 0.9*self.cfl_max
        except KeyError:
            raise KeyError('Maximum CFL number is not provided.')

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
    def step(self,solution):
        """Evolve q over one time step.

        Take one step with a Runge-Kutta or multistep method as specified by
        `solver.time_integrator`.
        """
        state = solution.states[0]

        self.before_step(self,state)

        try:
            if self.time_integrator=='Euler':
                deltaq=self.dq(state)
                state.q+=deltaq

            elif self.time_integrator=='SSP33':
                deltaq=self.dq(state)
                self._registers[0].q=state.q+deltaq
                self._registers[0].t =state.t+self.dt

                if self.call_before_step_each_stage:
                    self.before_step(self,self._registers[0])
                deltaq=self.dq(self._registers[0])
                self._registers[0].q= 0.75*state.q + 0.25*(self._registers[0].q+deltaq)
                self._registers[0].t = state.t+0.5*self.dt

                if self.call_before_step_each_stage:
                    self.before_step(self,self._registers[0])
                deltaq=self.dq(self._registers[0])
                state.q = 1./3.*state.q + 2./3.*(self._registers[0].q+deltaq)


            elif self.time_integrator=='SSP104':
                state.q = self.ssp104(state)


            elif self.time_integrator=='RK':
                # General RK with specified coefficients
                # self._registers[i].q actually stores dt*f(y_i)
                num_stages = len(self.b)
                for i in range(num_stages):
                    self._registers[i].q = state.q.copy()
                    for j in range(i):
                        self._registers[i].q += self.a[i,j]*self._registers[j].q
                    self._registers[i].t = state.t + self.dt * self.c[i]
                    self._registers[i].q = self.dq(self._registers[i])

                for j in range(num_stages):
                    state.q += self.b[j]*self._registers[j].q


            elif self.time_integrator == 'SSPMS32':
                # Store initial solution
                if self.step_index == 1:
                    for i in range(2):
                        self._registers[-2+i].dt = self.dt
                    self._registers[-1].q = state.q.copy()

                if self.step_index < 3:
                    # Use Euler method for starting values
                    deltaq = self.dq(state)
                    state.q += deltaq
                    self.step_index += 1
                
                else:
                    omega = (self._registers[-2].dt + self._registers[-1].dt)/self.dt
                    # ssp coefficient
                    r = (omega-1.)/omega
                    # method coefficients 
                    delta = 1./omega**2
                    beta = (omega+1.)/omega
                    deltaq = self.dq(state)
                    state.q = beta*(r*state.q + deltaq) + delta*self._registers[-3].q

                # Update stored solutions
                for i in range(2):
                    self._registers[-3+i].q = self._registers[-2+i].q.copy()
                    self._registers[-3+i].dt = self._registers[-2+i].dt
                self._registers[-1].q = state.q.copy()
                self._registers[-1].dt = self.dt


            elif self.time_integrator == 'SSPMS43':
                # Store initial solution
                if self.step_index == 1:
                    for i in range(3):
                        self._registers[-3+i].dt = self.dt
                    self._registers[-1].q = state.q.copy()

                if self.step_index < 4:
                    # Use SSP22 method for starting values
                    import copy
                    s1 = copy.deepcopy(state)
                    s1.set_num_ghost(self.num_ghost)

                    deltaq=self.dq(state)
                    s1.q = state.q + deltaq
                    s1.t = state.t + self.dt
                    deltaq = self.dq(s1)
                    state.q = 0.5*(state.q + s1.q + deltaq)

                    self.step_index += 1
                
                else:
                    H = self._registers[-3].dt + self._registers[-2].dt + self._registers[-1].dt
                    omega3 = H/self.dt
                    omega4 = omega3 + 1.
                    # SSP coefficient
                    r = (omega3-2.)/omega3
                    # method coefficients
                    delta0 = (4*omega4 - omega3**2)/omega3**3
                    beta0 = omega4/omega3**2
                    beta3 = omega4**2/omega3**2
                    deltaq = self.dq(state)
                    deltaqm4 = self.dq(self._registers[-4]) 
                    state.q = beta3*(r*state.q + deltaq) + \
                            (r*beta0+delta0)*self._registers[-4].q + beta0*deltaqm4

                # Update stored solutions
                for i in range(3):
                    self._registers[-4+i].q = self._registers[-3+i].q.copy()
                    self._registers[-4+i].dt = self._registers[-3+i].dt
                self._registers[-1].q = state.q.copy()
                self._registers[-1].dt = self.dt


            elif self.time_integrator == 'LMM':
                num_steps = len(self.alpha)

                # Store initial solution
                if self.step_index == 1:
                    self._registers[-num_steps].q  = state.q.copy()
                    self._registers[-num_steps].dq = self.dq(state)

                if self.step_index < num_steps:
                    # Using SSP104 for previous step values
                    state.q = self.ssp104(state)
                    self._registers[-num_steps+self.step_index].q = state.q.copy()
                    self._registers[-num_steps+self.step_index].dq = self.dq(state)
                    self.step_index += 1
                else:
                    # Update solution: alpha[-1] and beta[-1] correspond to solution at the previous step
                    state.q = self.alpha[-1]*self._registers[-1].q + self.beta[-1]*self._registers[-1].dq
                    for i in range(-num_steps,-1):
                        state.q += self.alpha[i]*self._registers[i].q + self.beta[i]*self._registers[i].dq
                        self._registers[i].q = self._registers[i+1].q.copy()
                        self._registers[i].dq = self._registers[i+1].dq.copy()
                    # Store current solution and function evaluation
                    self._registers[-1].q = state.q.copy()
                    self._registers[-1].dq = self.dq(state)


            else:
                raise Exception('Unrecognized time integrator')
        except CFLError:
            return False


    def ssp104(self,state):
        if self.time_integrator == 'SSP104':
            s1=self._registers[0]
            s2=self._registers[1]
            s1.q = state.q.copy()
        elif self.time_integrator == 'LMM':
            import copy
            s1 = copy.deepcopy(state)
            s1.set_num_ghost(self.num_ghost)
            s2 = copy.deepcopy(s1)

        deltaq=self.dq(state)
        s1.q = state.q + deltaq/6.
        s1.t = state.t + self.dt/6.

        for i in xrange(4):
            if self.call_before_step_each_stage:
                self.before_step(self,s1)
            deltaq=self.dq(s1)
            s1.q=s1.q + deltaq/6.
            s1.t =s1.t + self.dt/6.

        s2.q = state.q/25. + 9./25 * s1.q
        s1.q = 15. * s2.q - 5. * s1.q
        s1.t = state.t + self.dt/3.

        for i in xrange(4):
            if self.call_before_step_each_stage:
                self.before_step(self,s1)
            deltaq=self.dq(s1)
            s1.q=s1.q + deltaq/6.
            s1.t =s1.t + self.dt/6.

        if self.call_before_step_each_stage:
            self.before_step(self,s1)
        deltaq = self.dq(s1)
 
        return s2.q + 0.6 * s1.q + 0.1 * deltaq


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

        # Check here if we violated the CFL condition, if we did, return 
        # immediately to evolve_to_time and let it deal with picking a new
        # dt
        if self.cfl.get_cached_max() > self.cfl_max:
            raise CFLError('cfl_max exceeded')

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
        deltaq = self.dq_hyperbolic(state)

        if self.dq_src is not None:
            deltaq+=self.dq_src(self,state,self.dt)

        return deltaq.flatten('f')


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
        elif self.time_integrator == 'SSP104':  nregisters=2
        elif self.time_integrator == 'RK':      nregisters=len(self.b)+1
        elif self.time_integrator == 'SSPMS32': nregisters=3
        elif self.time_integrator == 'SSPMS43': nregisters=4
        elif self.time_integrator == 'LMM':
            nregisters=len(self.alpha)
            self.dt_variable = False
        else:
            raise Exception('Unrecognized time intergrator')
        
        state = solution.states[0]
        # use the same class constructor as the solution for the Runge Kutta stages
        State = type(state)
        self._registers = []
        for i in xrange(nregisters):
            #Maybe should use State.copy() here?
            self._registers.append(State(state.patch,state.num_eqn,state.num_aux))
            self._registers[-1].problem_data                = state.problem_data
            self._registers[-1].set_num_ghost(self.num_ghost)
            self._registers[-1].t                           = state.t
            if state.num_aux > 0: self._registers[-1].aux   = state.aux


    def get_cfl_max(self):
        """
        Set maximum CFL number for current step depending on time integrator
        """
        if self.time_integrator[:-2] == 'SSPMS' and self.step_index >= len(self._registers):
            s = len(self._registers)-2
            H = self._registers[-2].dt
            for i in range(s):
                H += self._registers[-3-i].dt
            r = (H - s*self._registers[-1].dt)/H # ssp coefficient at the current step
            sigma = r/self._sspcoeff[self.time_integrator]
        else:
            sigma = 1.0

        return sigma * self.cfl_max

    def get_dt_new(self):
        """
        Set time-step for next step depending on time integrator
        """
        # desired time-step 
        dt_des = self.dt * self.cfl_desired / self.cfl.get_cached_max()

        if self.time_integrator[:-2] == 'SSPMS' and self.step_index >= len(self._registers):
            s = len(self._registers)-2
            H = self._registers[-1].dt
            for i in range(s):
                H += self._registers[-2-i].dt
            sigma = H / (self._sspcoeff[self.time_integrator]*H + s*dt_des)
        else:
            sigma = 1.0

        return sigma * dt_des



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
                dtdx = self.dt / (grid.delta[0] * state.aux[state.index_capa,:])
            else:
                dtdx += self.dt/grid.delta[0]
 
            aux=self.auxbc
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
            for mw in xrange(self.num_waves):
                smax1 = np.max( dtdx[LL  :UL]  *s[mw,LL-1:UL-1])
                smax2 = np.max(-dtdx[LL-1:UL-1]*s[mw,LL-1:UL-1])
                cfl = max(cfl,smax1,smax2)

            #Find total fluctuation within each cell
            wave,s,amdq2,apdq2 = self.rp(ql,qr,aux,aux,state.problem_data)

            # Compute dq
            for m in xrange(state.num_eqn):
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
