r"""
Module containing SharpClaw solvers for PyClaw/PetClaw

#  File:        sharpclaw.py
#  Created:     2010-03-20
#  Author:      David Ketcheson
"""
# Solver superclass
from ..solver import Solver, CFLError

# Reconstructor
try:
    # load c-based WENO reconstructor (PyWENO)
    from ..limiters import reconstruct as recon
except ImportError:
    # load old WENO5 reconstructor
    from ..limiters import recon


def before_step(solver,solution):
    r"""
    Dummy routine called before each step
    
    Replace this routine if you want to do something before each time step.
    """
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

        Time integrator to be used.
        Euler: forward Euler method.
        SSP33: 3-stages, 3rd-order SSP Runge-Kutta method.
        SSP104: 10-stages, 4th-order SSP Runge-Kutta method.
        ``Default = 'SSP104'``

    .. attribute:: char_decomp

        Type of WENO reconstruction.
        0: conservative variables WENO reconstruction (standard).
        1: characteristic-wise WENO reconstruction.
        2: transmission-based WENO reconstruction.
        ``Default = 0``

    .. attribute:: tfluct_solver

        Whether a total fluctuation solver have to be used. If True the function
        that calculates the total fluctuation must be provided.
        ``Default = False``

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
        self.aux_time_dep = False
        self.kernel_language = 'Fortran'
        self.num_ghost = 3
        self.fwave = False
        self.cfl_desired = 2.45
        self.cfl_max = 2.5
        self.dq_src = None
        self.call_before_step_each_stage = False
        self._mthlim = self.limiters
        self._method = None
        self._rk_stages = None
        

        # Call general initialization function
        super(SharpClawSolver,self).__init__(riemann_solver,claw_package)
        
    def setup(self,solution):
        """
        Allocate RK stage arrays and fortran routine work arrays.
        """
        self.num_ghost = (self.weno_order+1)/2

        # This is a hack to deal with the fact that petsc4py
        # doesn't allow us to change the stencil_width (num_ghost)
        state = solution.state
        state.set_num_ghost(self.num_ghost)
        # End hack

        self.allocate_rk_stages(solution)
        self.set_mthlim()

        state = solution.states[0]
 
        if self.kernel_language=='Fortran':
            if self.fmod is None:
                so_name = 'clawpack.pyclaw.sharpclaw.sharpclaw'+str(self.num_dim)
                self.fmod = __import__(so_name,fromlist=['clawpack.pyclaw.sharpclaw'])
            state.set_cparam(self.fmod)
            state.set_cparam(self.rp)
            self.set_fortran_parameters(state,self.fmod.clawparams,self.fmod.workspace,self.fmod.reconstruct)

        self.allocate_bc_arrays(state)

        self._is_set_up = True

    # ========== Time stepping routines ======================================
    def step(self,solution):
        """Evolve q over one time step.

        Take on Runge-Kutta time step using the method specified by
        self..time_integrator.  Currently implemented methods:

        'Euler'  : 1st-order Forward Euler integration
        'SSP33'  : 3rd-order strong stability preserving method of Shu & Osher
        'SSP104' : 4th-order strong stability preserving method Ketcheson
        """
        state = solution.states[0]

        self.before_step(self,state)

        try:
            if self.time_integrator=='Euler':
                deltaq=self.dq(state)
                state.q+=deltaq

            elif self.time_integrator=='SSP33':
                deltaq=self.dq(state)
                self._rk_stages[0].q=state.q+deltaq
                self._rk_stages[0].t =state.t+self.dt

                if self.call_before_step_each_stage:
                    self.before_step(self,self._rk_stages[0])
                deltaq=self.dq(self._rk_stages[0])
                self._rk_stages[0].q= 0.75*state.q + 0.25*(self._rk_stages[0].q+deltaq)
                self._rk_stages[0].t = state.t+0.5*self.dt

                if self.call_before_step_each_stage:
                    self.before_step(self,self._rk_stages[0])
                deltaq=self.dq(self._rk_stages[0])
                state.q = 1./3.*state.q + 2./3.*(self._rk_stages[0].q+deltaq)


            elif self.time_integrator=='SSP104':
                s1=self._rk_stages[0]
                s2=self._rk_stages[1]
                s1.q = state.q.copy()

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
                state.q = s2.q + 0.6 * s1.q + 0.1 * deltaq
            else:
                raise Exception('Unrecognized time integrator')
        except CFLError:
            return False


    def set_mthlim(self):
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


    def set_fortran_parameters(self,state,clawparams,workspace,reconstruct):
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

    def allocate_rk_stages(self,solution):
        r"""
        Instantiate State objects for Runge--Kutta stages.

        This routine is only used by method-of-lines solvers (SharpClaw),
        not by the Classic solvers.  It allocates additional State objects
        to store the intermediate stages used by Runge--Kutta time integrators.

        If we create a MethodOfLinesSolver subclass, this should be moved there.
        """
        if self.time_integrator   == 'Euler':  nregisters=1
        elif self.time_integrator == 'SSP33':  nregisters=2
        elif self.time_integrator == 'SSP104': nregisters=3
 
        state = solution.states[0]
        # use the same class constructor as the solution for the Runge Kutta stages
        State = type(state)
        self._rk_stages = []
        for i in xrange(nregisters-1):
            #Maybe should use State.copy() here?
            self._rk_stages.append(State(state.patch,state.num_eqn,state.num_aux))
            self._rk_stages[-1].problem_data       = state.problem_data
            self._rk_stages[-1].set_num_ghost(self.num_ghost)
            self._rk_stages[-1].t                = state.t
            if state.num_aux > 0:
                self._rk_stages[-1].aux              = state.aux




# ========================================================================
class SharpClawSolver1D(SharpClawSolver):
# ========================================================================
    """
    SharpClaw solver for one-dimensional problems.
    
    Used to solve 1D hyperbolic systems using the SharpClaw algorithms,
    which are based on WENO reconstruction and Runge-Kutta time stepping.
    """
    def __init__(self,riemann_solver=None,claw_package=None):
        r"""
        See :class:`SharpClawSolver1D` for more info.
        """   
        self.num_dim = 1
        super(SharpClawSolver1D,self).__init__(riemann_solver,claw_package)


    def teardown(self):
        r"""
        Deallocate F90 module arrays.
        Also delete Fortran objects, which otherwise tend to persist in Python sessions.
        """
        if self.kernel_language=='Fortran':
            self.fmod.clawparams.dealloc_clawparams()
            self.fmod.workspace.dealloc_workspace(self.char_decomp)
            self.fmod.reconstruct.dealloc_recon_workspace(self.fmod.clawparams.lim_type,self.fmod.clawparams.char_decomp)
            del self.fmod


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

        self.apply_q_bcs(state)
        if state.num_aux > 0:
            self.apply_aux_bcs(state)
        q = self.qbc 

        grid = state.grid
        mx = grid.num_cells[0]

        ixy=1

        if self.kernel_language=='Fortran':
            rp1 = self.rp.rp1._cpointer
            dq,cfl=self.fmod.flux1(q,self.auxbc,self.dt,state.t,ixy,mx,self.num_ghost,mx,rp1)

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
    """SharpClaw evolution routine in 2D
    
    This class represents the 2D SharpClaw solver.  Note that there are 
    routines here for interfacing with the fortran time stepping routines only.
    """
    def __init__(self,riemann_solver=None,claw_package=None):
        r"""
        Create 2D SharpClaw solver
        
        See :class:`SharpClawSolver2D` for more info.
        """   
        self.num_dim = 2

        super(SharpClawSolver2D,self).__init__(riemann_solver,claw_package)


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
        self.apply_q_bcs(state)
        if state.num_aux > 0:    
            self.apply_aux_bcs(state)
        q = self.qbc 

        grid = state.grid

        num_ghost=self.num_ghost
        mx=grid.num_cells[0]
        my=grid.num_cells[1]
        maxm = max(mx,my)

        if self.kernel_language=='Fortran':
            rpn2 = self.rp.rpn2._cpointer
            dq,cfl=self.fmod.flux2(q,self.auxbc,self.dt,state.t,num_ghost,maxm,mx,my,rpn2)

        else: raise Exception('Only Fortran kernels are supported in 2D.')

        self.cfl.update_global_max(cfl)
        return dq[:,num_ghost:-num_ghost,num_ghost:-num_ghost]
