r"""
Module containing SharpClaw solvers for PyClaw/PetClaw

#  File:        sharpclaw.py
#  Created:     2010-03-20
#  Author:      David Ketcheson
"""
# Solver superclass
from pyclaw.solver import Solver, CFLError

# Reconstructor
try:
    # load c-based WENO reconstructor (PyWENO)
    from pyclaw.limiters import reconstruct as recon
except:
    # load old WENO5 reconstructor
    from pyclaw.limiters import recon


def start_step(solver,solutions):
    r"""
    Dummy routine called before each step
    
    Replace this routine if you want to do something before each time step.
    """
    pass

def src(solver,grid,q,t):
    r"""
    Dummy routine called to calculate a source term
    
    Replace this routine if you want to include a source term.
    """
    pass

 
class SharpClawSolver(Solver):
    r"""
    Superclass for all SharpClawND solvers.

    Implements Runge-Kutta time stepping and the basic form of a 
    semi-discrete step (the dq() function).  If another method-of-lines
    solver is implemented in the future, it should be based on this class,
    which then ought to be renamed to something like "MOLSolver".
    """
    
    # ========================================================================
    #   Initialization routines
    # ========================================================================
    def __init__(self, data=None):
        r"""
        Set default options for SharpClawSolvers and call the super's __init__().
        """
        
        # Required attributes for this solver
        for attr in ['limiters','start_step','lim_type','time_integrator',
                     'char_decomp','src_term','aux_time_dep','mwaves']:
            self._required_attrs.append(attr)
        
        # Defaults for required attributes
        self._default_attr_values['limiters'] = [1]
        self._default_attr_values['start_step'] = start_step
        self._default_attr_values['lim_type'] = 2
        self._default_attr_values['time_integrator'] = 'SSP104'
        self._default_attr_values['char_decomp'] = 0
        self._default_attr_values['tfluct_solver'] = False
        self._default_attr_values['aux_time_dep'] = False
        self._default_attr_values['src_term'] = False
        self._default_attr_values['kernel_language'] = 'Fortran'
        self._default_attr_values['mbc'] = 3
        self._default_attr_values['fwave'] = False
        self._default_attr_values['cfl_desired'] = 2.45
        self._default_attr_values['cfl_max'] = 2.5
        
        # Call general initialization function
        super(SharpClawSolver,self).__init__(data)
        
    # ========== Time stepping routines ======================================
    def step(self,solutions):
        """Evolve q over one time step.

        Take on Runge-Kutta time step using the method specified by
        self..time_integrator.  Currently implemented methods:

        'Euler'  : 1st-order Forward Euler integration
        'SSP33'  : 3rd-order strong stability preserving method of Shu & Osher
        'SSP104' : 4th-order strong stability preserving method Ketcheson
        """
        from pyclaw.solution import Solution
        state = solutions['n'].states[0]

        self.start_step(self,solutions)

        try:
            if self.time_integrator=='Euler':
                deltaq=self.dq(state)
                state.q+=deltaq

            elif self.time_integrator=='SSP33':
                deltaq=self.dq(state)
                self.rk_stages[0].q=state.q+deltaq
                self.rk_stages[0].t =state.t+self.dt
                deltaq=self.dq(self.rk_stages[0])
                self.rk_stages[0].q= 0.75*state.q + 0.25*(self.rk_stages[0].q+deltaq)
                self.rk_stages[0].t = state.t+0.5*self.dt
                deltaq=self.dq(self.rk_stages[0])
                state.q = 1./3.*state.q + 2./3.*(self.rk_stages[0].q+deltaq)

            elif self.time_integrator=='SSP104':
                s1=self.rk_stages[0]
                s2=self.rk_stages[1]
                s1.q = state.q.copy()

                deltaq=self.dq(state)
                s1.q = state.q + deltaq/6.
                s1.t = state.t + self.dt/6.

                for i in range(4):
                    deltaq=self.dq(s1)
                    s1.q=s1.q + deltaq/6.
                    s1.t =s1.t + self.dt/6.

                s2.q = state.q/25. + 9./25 * s1.q
                s1.q = 15. * s2.q - 5. * s1.q
                s1.t = state.t + self.dt/3.

                for i in range(4):
                    deltaq=self.dq(s1)
                    s1.q=s1.q + deltaq/6.
                    s1.t =s1.t + self.dt/6.

                deltaq = self.dq(s1)
                state.q = s2.q + 0.6 * s1.q + 0.1 * deltaq
            else:
                raise Exception('Unrecognized time integrator')
        except CFLError:
            return False


    def set_mthlim(self):
        self.mthlim = self.limiters
        if not isinstance(self.limiters,list): self.mthlim=[self.mthlim]
        if len(self.mthlim)==1: self.mthlim = self.mthlim * self.mwaves
        if len(self.mthlim)!=self.mwaves:
            raise Exception('Length of solver.limiters is not equal to 1 or to solver.mwaves')
 
       
    def dq(self,state):
        """
        Evaluate dq/dt * (delta t)
        """

        deltaq = self.dq_homogeneous(state)

        # Check here if we violated the CFL condition, if we did, return 
        # immediately to evolve_to_time and let it deal with picking a new
        # dt
        self.communicateCFL()
        if self.cfl >= self.cfl_max:
            raise CFLError('cfl_max exceeded')

        # Godunov Splitting -- really the source term should be called inside rkstep
        if self.src_term == 1:
            deltaq+=self.src(state,q,state.t)

        return deltaq

    def dq_homogeneous(state):
        raise NotImplementedError('You must subclass SharpClawSolver.')

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
    
         
    def dqdt(self,state):
        """
        Evaluate dq/dt.  This routine is used for implicit time stepping.
        """

        self.dt = 1
        deltaq = self.dq_homogeneous(state)

        # Godunov Splitting -- really the source term should be called inside rkstep
        if self.src_term == 1:
            deltaq+=self.src(state.grid,q,state.t)

        return deltaq.flatten('f')


    def set_fortran_parameters(self,state,clawparams,workspace,reconstruct):
        """
        Set parameters for Fortran modules used by SharpClaw.
        The modules should be imported and passed as arguments to this function.

        """
        grid = state.grid
        clawparams.ndim          = grid.ndim
        clawparams.lim_type      = self.lim_type
        clawparams.char_decomp   = self.char_decomp
        clawparams.tfluct_solver = self.tfluct_solver
        clawparams.fwave         = self.fwave
        if state.capa is not None:
            clawparams.mcapa         = 1
        else:
            clawparams.mcapa         = 0

        clawparams.mwaves        = self.mwaves
        clawparams.alloc_clawparams()
        for idim in range(grid.ndim):
            clawparams.xlower[idim]=grid.dimensions[idim].lower
            clawparams.xupper[idim]=grid.dimensions[idim].upper
        clawparams.dx       =grid.d
        clawparams.mthlim   =self.mthlim

        maxnx = max(grid.ng)+2*self.mbc
        workspace.alloc_workspace(maxnx,self.mbc,state.meqn,self.mwaves,self.char_decomp)
        reconstruct.alloc_recon_workspace(maxnx,self.mbc,state.meqn,self.mwaves,
                                            clawparams.lim_type,clawparams.char_decomp)


# ========================================================================
class SharpClawSolver1D(SharpClawSolver):
# ========================================================================
    """
    SharpClaw solver for one-dimensional problems.
    
    Used to solve 1D hyperbolic systems using the SharpClaw algorithms,
    which are based on WENO reconstruction and Runge-Kutta time stepping.
    """
    def __init__(self, data=None):
        r"""
        See :class:`SharpClawSolver1D` for more info.
        """   
        self.ndim = 1
        super(SharpClawSolver1D,self).__init__(data)


    def setup(self,solutions):
        """
        Allocate RK stage arrays and fortran routine work arrays.
        """
        self.allocate_rk_stages(solutions)
        self.set_mthlim()
 
        if self.kernel_language=='Fortran':
            from sharpclaw1 import clawparams, workspace, reconstruct
            import sharpclaw1
            state = solutions['n'].states[0]
            state.set_cparam(sharpclaw1)
            self.set_fortran_parameters(state,clawparams,workspace,reconstruct)

    def teardown(self):
        r"""
        Deallocate F90 module arrays.
        """
        if self.kernel_language=='Fortran':
            from sharpclaw1 import clawparams, workspace, reconstruct
            clawparams.dealloc_clawparams()
            workspace.dealloc_workspace(self.char_decomp)
            reconstruct.dealloc_recon_workspace(clawparams.lim_type,clawparams.char_decomp)


    def dq_homogeneous(self,state):
        r"""
        Compute dq/dt * (delta t) for the homogeneous hyperbolic system.

        Note that the capa array, if present, should be located in the aux
        variable.

        Indexing works like this (here mbc=2 as an example)::

         0     1     2     3     4     mx+mbc-2     mx+mbc      mx+mbc+2
                     |                        mx+mbc-1 |  mx+mbc+1
         |     |     |     |     |   ...   |     |     |     |     |
            0     1  |  2     3            mx+mbc-2    |mx+mbc       
                                                  mx+mbc-1   mx+mbc+1

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

        q = self.qbc(state) # Can we make use of state.qbc instead
                            # of calling self.qbc() again?
        grid = state.grid
        mx = grid.ng[0]

        ixy=1
        aux=state.aux
        if state.maux == 0: 
            aux = np.empty( (state.maux,mx+2*self.mbc) ,order='F')
        else:
            aux = self.auxbc(state)

        if self.kernel_language=='Fortran':
            from sharpclaw1 import flux1
            dq,self.cfl=flux1(q,aux,self.dt,state.t,ixy,mx,self.mbc,mx)

        elif self.kernel_language=='Python':

            dtdx = np.zeros( (mx+2*self.mbc) ,order='F')
            dq   = np.zeros( (state.meqn,mx+2*self.mbc) ,order='F')

            # Find local value for dt/dx
            if state.capa is not None:
                dtdx = self.dt / (grid.d[0] * state.capa)
            else:
                dtdx += self.dt/grid.d[0]
 
            if aux is not None:
                aux_l=aux[:,:-1]
                aux_r=aux[:,1: ]
            else:
                aux_l = None
                aux_r = None

            #Reconstruct (wave reconstruction uses a Riemann solve)
            if self.lim_type==-1: #1st-order Godunov
                ql=q; qr=q;
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
                    wave,s,amdq,apdq = self.rp(q_l,q_r,aux_l,aux_r,state.aux_global)
                    ql,qr=recon.weno5_wave(q,wave,s)
                elif self.char_decomp==2: #Characteristic-wise reconstruction
                    raise NotImplementedError

            # Solve Riemann problem at each interface
            q_l=qr[:,:-1]
            q_r=ql[:,1: ]
            wave,s,amdq,apdq = self.rp(q_l,q_r,aux_l,aux_r,state.aux_global)

            # Loop limits for local portion of grid
            # THIS WON'T WORK IN PARALLEL!
            LL = self.mbc - 1
            UL = grid.ng[0] + self.mbc + 1

            # Compute maximum wave speed
            self.cfl = 0.0
            for mw in xrange(self.mwaves):
                smax1 = np.max( dtdx[LL  :UL]  *s[mw,LL-1:UL-1])
                smax2 = np.max(-dtdx[LL-1:UL-1]*s[mw,LL-1:UL-1])
                self.cfl = max(self.cfl,smax1,smax2)

            #Find total fluctuation within each cell
            wave,s,amdq2,apdq2 = self.rp(ql,qr,aux,aux,state.aux_global)

            # Compute dq
            for m in xrange(state.meqn):
                dq[m,LL:UL] = -dtdx[LL:UL]*(amdq[m,LL:UL] + apdq[m,LL-1:UL-1] \
                                + apdq2[m,LL:UL] + amdq2[m,LL:UL])

        else: raise Exception('Unrecognized value of solver.kernel_language.')
        
        return dq[:,self.mbc:-self.mbc]
    

# ========================================================================
class SharpClawSolver2D(SharpClawSolver):
# ========================================================================
    """SharpClaw evolution routine in 2D
    
    This class represents the 2D SharpClaw solver.  Note that there are 
    routines here for interfacing with the fortran time stepping routines only.
    """
    def __init__(self, data=None):
        r"""
        Create 2D SharpClaw solver
        
        See :class:`SharpClawSolver2D` for more info.
        """   
        
        self.ndim = 2

        super(SharpClawSolver2D,self).__init__(data)


    def setup(self,solutions):
        """
        Allocate RK stage arrays and fortran routine work arrays.
        """
        self.allocate_rk_stages(solutions)
        self.set_mthlim()
 
        if self.kernel_language=='Fortran':
            from sharpclaw2 import clawparams, workspace, reconstruct
            import sharpclaw2
            state = solutions['n'].states[0]
            state.set_cparam(sharpclaw2)
            self.set_fortran_parameters(state,clawparams,workspace,reconstruct)


    def teardown(self):
        r"""
        Deallocate F90 module arrays.
        """
        if self.kernel_language=='Fortran':
            from sharpclaw2 import clawparams, workspace, reconstruct
            workspace.dealloc_workspace(self.char_decomp)
            reconstruct.dealloc_recon_workspace(clawparams.lim_type,clawparams.char_decomp)
            clawparams.dealloc_clawparams()


    def dq_homogeneous(self,state):
        """Compute dq/dt * (delta t) for the homogeneous hyperbolic system

        Note that the capa array, if present, should be located in the aux
        variable.

        Indexing works like this (here mbc=2 as an example)::

         0     1     2     3     4     mx+mbc-2     mx+mbc      mx+mbc+2
                     |                        mx+mbc-1 |  mx+mbc+1
         |     |     |     |     |   ...   |     |     |     |     |
            0     1  |  2     3            mx+mbc-2    |mx+mbc       
                                                  mx+mbc-1   mx+mbc+1

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

        q = self.qbc(state) # Can we make use of state.qbc instead
                            # of calling self.qbc() again?


        grid = state.grid

        mbc=self.mbc
        mx=grid.ng[0]
        my=grid.ng[1]
        maxm = max(mx,my)

        aux=state.aux
        if state.maux == 0:
            aux = np.empty( (state.maux,mx+2*mbc,my+2*mbc), order='F' )
        else:
            aux = self.auxbc(state)

        if self.kernel_language=='Fortran':
            from sharpclaw2 import flux2
            dq,self.cfl=flux2(q,aux,self.dt,state.t,mbc,maxm,mx,my)

        else: raise Exception('Only Fortran kernels are supported in 2D.')

        return dq[:,mbc:-mbc,mbc:-mbc]
