r"""
This module contains the implicit sharpclaw solvers which requires 
fsolve of scipy.optimize nonlinear solver.

All implicit sharpclaw solvers inherit from the 
:class:`ImplicitSharpClawSolverfsolve` superclass which in turn inherits from the 
:class:`~petclaw.solver.Solver` superclass. These are both pure virtual classes;
the only solver classes that should be instantiated are the dimension-specific
ones, :class:`ImplicitSharpClawSolverfsolve1D` and 
:class:`ImplicitSharpClawSolverfsolve2D`.

Authors: Matteo Parsani
         David Ketcheson
"""

# Solver superclass
import petclaw.solver

# Riemann solvers
import riemann

from scipy.optimize import fsolve
from scipy.optimize.nonlin import newton_krylov


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

# ============================================================================
#  Generic implicit SharpClaw solver class
# ============================================================================
class ImplicitSharpClawSolverfsolve(petclaw.solver.PetSolver):
    r"""
    Superclass for all ImplicitSharpClawND solvers.

    Implements implicit time stepping methods combined with the PETSc's 
    nonlinear solver SNES. If another method-of-lines solver is implemented in 
    the future, it should be based on this class,which then ought to be renamed 
    to something like "ImplicitMOLSolver".
    """
    
    # ========================================================================
    #   Initialization routines
    # ========================================================================
    def __init__(self, data=None):
        r"""
        Set default options for ImplicitSharpClawSolvers and call the 
        super's __init__().
        """
        
        # Required attributes for this solver
        for attr in ['limiters','start_step','lim_type','time_integrator',
                     'char_decomp','src_term','aux_time_dep','mwaves']:
            self._required_attrs.append(attr)
        
        # Defaults for required attributes
        self._default_attr_values['limiters'] = [1]
        self._default_attr_values['start_step'] = start_step
        self._default_attr_values['lim_type'] = 2
        self._default_attr_values['time_integrator'] = 'BEuler'
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
        super(ImplicitSharpClawSolverfsolve,self).__init__(data)


    # ========== Setup Routine ===============================================
    def initiate(self,solutions):
        r"""
        Called before any set of time steps.
        
        This routine will be called once before the solver is used via the
        :class:`~pyclaw.controller.Controller`.
        """
        
        # Import libraries
        from petsc4py import PETSc
        from numpy import empty
        from pyclaw.state import State

        # Call b4step
        self.start_step(self,solutions)

        # Get state
        state = solutions['n'].state
    

        # This is a hack to deal with the fact that petsc4py
        # doesn't allow us to change the stencil_width (mbc)
        state.set_stencil_width(self.mbc)
        # End hack

        # Set mthlim
        self.set_mthlim()


        # Create PETSc vectors in charge of containig:
        # bVec: the constant part of the nonlinear algebraic system of equations
        # fVec: nonlinear vector-valued function
        self.bVec    = state.gqVec.duplicate()
        self.fVec    = state.gqVec.duplicate()

        # Create PETSc nonlinear solver
        self.snes    = PETSc.SNES().create()
        #self.snes.setJacobian(None, self.Jac) 

        # Ought to implement a copy constructor for State
        self.impsol_stage = State(state.grid)
        self.impsol_stage.meqn             = state.meqn
        self.impsol_stage.maux             = state.maux
        self.impsol_stage.aux_global       = state.aux_global
        self.impsol_stage.t                = state.t
        if state.maux > 0:
            self.impsol_stage.aux          = state.aux



    # ========== Time stepping routines ======================================
    def step(self,solutions):
        r"""
        Evolve q for one time step using the method specified by
        self.time_integrator. Currently implemented methods:

        'BEuler': 1st-order backward Euler integration


        This routine is called from the method evolve_to_time defined in the
        pyclaw.solver.Solver superclass. 

        :Input:
         - *solutions* - (:class:`~pyclaw.solution.Solution`) Dictionary of 
           solutions to be evolved
         
        :Output: 
         - (bool) - True if full step succeeded, False otherwise
        """
        
        from pyclaw.solution import Solution

        # Get state object
        state = solutions['n'].states[0]

        
        # Call b4step, pyclaw should be subclassed if this is needed
        # NOTE: THIS IS NOT THE RIGHT PLACE WHEN WE HAVE AN IMPLICIT  
        # TIME STEPPING BECAUSE, FOR INSTANCE, THE AUX ARRAY COULD DEPEND ON
        # THE SOLUTION ITSELF.
        # THEN start_step FUNCTION SHOULD BE PLACED IN THE PART OF THE CODE 
        # WHERE THE NONLINEAR FUNCTION IS COMPUTED!!!!!!!!!
        self.start_step(self,solutions)

        self.updatesolution(solutions)
        # Check here if we violated the CFL condition, if we did, return 
        # immediately to evolve_to_time and let it pick pick up a new dt. 
        # Even for steady state calculations the control of the CFL is important
        # because in the first few pseudo time steps, where high frequency 
        # errors components must be damped out and expelled from the 
        # computational domain. 
        # TODO: implement a CFL-law
        self.communicateCFL()

        if self.cfl >= self.cfl_max:
            return False
        else:
            return True

    def updatesolution(self,solutions):
        r"""
        Updates solution.
        """
        raise Exception("Dummy routine, please override!")

    def set_mthlim(self):
        self.mthlim = self.limiters
        if not isinstance(self.limiters,list): self.mthlim=[self.mthlim]
        if len(self.mthlim)==1: self.mthlim = self.mthlim * self.mwaves
        if len(self.mthlim)!=self.mwaves:
            raise Exception('Length of solver.limiters is not equal to 1 or to solver.mwaves')


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


# ============================================================================
#  Implicit SharpClaw 1d Solver Class
# ============================================================================
class ImplicitSharpClawSolverfsolve1D(ImplicitSharpClawSolverfsolve):
    """
    Implicit SharpClaw solver for one-dimensional problems.
    
    Used to solve 1D hyperbolic systems using WENO reconstruction and implicit 
    time stepping technique.
    """

    def __init__(self,data=None):
        r"""
        Create 1d implicit SharpClawpack solver.
        
        See :class:`ImplicitSharpClawSolver1D` for more info.
        """   
        
        self.ndim = 1

        # Call superclass __init__
        super(ImplicitSharpClawSolverfsolve1D,self).__init__(data)


    # ========== Setup routine =============================   
    def setup(self,solutions):
        r"""
        Perform essential solver setup. This routine must be called before
        solver.step() may be called.

        Set Fortran data structures (for Clawpack). 
        """
        
        # Call parent's "setup" function
        self.initiate(solutions)

        # Set Fortran data structure for the 1D implicit SharpClaw solver
        if(self.kernel_language == 'Fortran'):
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


    def updatesolution(self,solutions):
        r"""
        Updates solution using fsolve.
        """

        # Import libraries
        import numpy as np
        from numpy import zeros, reshape, empty

        # Get state
        state = solutions['n'].states[0]
        
        # Store q at the current time level. It will be used later to construct 
        # the nonlinear vector-valued function
        qold = state.q.copy()
    
        # Set initial guess for the nonlinear solver
        qinit = qold.copy()


        # Define lambda function.
        nonlinearFun = lambda y : self.nonlinearFunctionBEuler(qold,y,state)
        qnew,info,ier,msg = fsolve(nonlinearFun,qinit[0,:],full_output=1)

        # 
        state.q[0,:] = qnew[:]
        


    def nonlinearFunctionBEuler(self,qold,y,state):
        r"""
        Computes the nonlinear function for the backward Euler scheme.

        :Input:
         - *qin* - Current approximation of the solution at the next time level,
         i.e. solution of the previous nonlinear solver's iteration.
        """

        # Import libraries
        import numpy as np
        from numpy import zeros, reshape, empty

        mx = state.grid.ng[0]
        dx = state.grid.d[0]
        mbc = self.mbc
        dt = self.dt

        
        # Define and set to zero the ratio between dt and dx 
        dtdx = np.zeros( (mx+2*mbc) ) + dt/dx

        # Auxbc is set here and not outside of this function because it is 
        # more general. Indeed aux could depend on q which changes at each 
        # nonlinear iteration!
        if state.maux>0:
            state.aux = self.auxbc(state)
        else:
            aux=np.empty((state.maux,mx+2*mbc), order='F')

        # Have to do this because of how qbc works...
        state.q = reshape(y,(state.meqn,mx),order='F') 
        qapprox = self.qbc(state)


        #print qapprox

        # Import library
        from sharpclaw1 import flux1

        # Call fortran routine 
        ixy = 1
        sd,self.cfl=flux1(qapprox,aux,dt,state.t,ixy,mx,mbc,mx)

        #print qold
        # Return the nonlinear vector-valued function
        return qapprox[0,mbc:-mbc]-sd[0,mbc:-mbc]-qold[0,:]
