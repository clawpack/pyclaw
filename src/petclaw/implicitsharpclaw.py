r"""
This module contains the implicit sharpclaw solvers which requires PETSc toolkit.

Currently only the bacward Euler (BE) scheme is implemented for the time 
discretization. This is the initial version and the structure of the solver is
based on BE discretization. Next step must be the generalization for "any"
implicit TS.


All implicit sharpclaw solvers inherit from the :class:`ImplicitSharpClawSolver` 
superclass which in turn inherits from the :class:`~petclaw.solver.Solver` superclass. 
These are both pure virtual classes; the only solver classes that should be instantiated
are the dimension-specific ones, :class:`ImplicitSharpClawSolver1D` and 
:class:`ImplicitSharpClawSolver2D`.

Authors: Matteo Parsani
         David Ketcheson
"""
# Solver superclass
import petclaw.solver

# Riemann solvers
import riemann


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
class ImplicitSharpClawSolver(petclaw.solver.PetSolver):
    r"""
    Superclass for all ImplicitSharpClawND solvers.

    Implements backward Euler time stepping combined with the PETSc's nonlinear
    algebraic solver library SNES. If another method-of-lines solver is 
    implemented in the future, it should be based on this class,which then ought 
    to be renamed to something like "ImplicitMOLSolver".
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
        super(ImplicitSharpClawSolver,self).__init__(data)


    # ========== Time stepping routines ======================================
    def step(self,solutions):
        r"""
        Evolve q for one time step.

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
        # TIME STEPPING BECAUSE, FOR INSTANCE, THE AUX ARRAY COULD DEPEND ON THE SOLUTION ITSELF.
        # THEN start_step FUNCTION SHOULD BE PLACED IN THE PART OF THE CODE WHERE THE
        # NONLINEAR FUNCTION IS COMPUTED!!!!!!!!!
        self.start_step(self,solutions)


        # Compute solution at the new time level
        self.updatesolution(state)

        # Check here if we violated the CFL condition, if we did, return 
        # immediately to evolve_to_time and let it deal with picking a new dt. 
        # Even for steady state calculations the control of the CFL is important, especially
        # in the first few pseudo time steps, where high frequency errors components must
        # be damped out and expelled from the computational domain. 
        # TODO: implement a CFL-law
        self.communicateCFL()


    def updatesolution(self,state):
        r"""
        Compute slution at the new time level for the implicit
        SharpClaw solver.
        
        This is a dummy routine and must be overridden.
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
#  Implicit ClawPack 1d Solver Class
# ============================================================================
class ImplicitSharpClawSolver1D(ImplicitSharpClawSolver):
    """
    Implicit SharpClaw solver for one-dimensional problems.
    
    Used to solve 1D hyperbolic systems using WENO reconstruction and implicit 
    time stepping technique.
    """

    def __init__(self,data=None):
        r"""
        Create 1d implicit SharpClawpack solver
        
        See :class:`ImplicitSharpClawSolver1D` for more info.
        """   
        
        self.ndim = 1

        super(ImplicitSharpClawSolver1D,self).__init__(data)


    # ========== Setup routine =============================   
    def setup(self,solutions):
        r"""
        Perform essential solver setup. This routine must be called before
        solver.step() may be called.

        Set Fortran data structures (for Clawpack) and set up a DA with
        the appropriate stencil width.
        """
        from petsc4py import PETSc
        from numpy import empty
        from pyclaw.state import State

        # This is a hack to deal with the fact that petsc4py
        # doesn't allow us to change the stencil_width (mbc)
        state = solutions['n'].state
        state.set_stencil_width(self.mbc)
        # End hack

        self.set_mthlim()
        if(self.kernel_language == 'Fortran'):
            from sharpclaw1 import clawparams, workspace, reconstruct
            import sharpclaw1
            state = solutions['n'].states[0]
            state.set_cparam(sharpclaw1)
            self.set_fortran_parameters(state,clawparams,workspace,reconstruct)

        self.bVec    = state.gqVec.duplicate()
        self.fVec    = state.gqVec.duplicate()
        self.snes    = PETSc.SNES().create()

        #Ought to implement a copy constructor for State
        self.impsol_stage = State(state.grid)
        self.impsol_stage.meqn             = state.meqn
        self.impsol_stage.maux             = state.maux
        self.impsol_stage.aux_global       = state.aux_global
        self.impsol_stage.t                = state.t
        if state.maux > 0:
            self.impsol_stage.aux          = state.aux

    def teardown(self):
        r"""
        Deallocate F90 module arrays.
        """
        if self.kernel_language=='Fortran':
            from sharpclaw1 import clawparams, workspace, reconstruct
            clawparams.dealloc_clawparams()
            workspace.dealloc_workspace(self.char_decomp)
            reconstruct.dealloc_recon_workspace(clawparams.lim_type,clawparams.char_decomp)



    def updatesolution(self,state):
        r"""
        Compute slution at the new time level for the implicit 1D
        Lax-Wendroff scheme, i.e.

        q^(n+1) = q^(n) + R(q^(n+1)),

        where q^(n)      = solution at the current time level
              q^(n+1)    = solution at the new time level (solution of the nonlinear system)
              R(q^(n+1)) = nonlinear function arising from the spatial/time discretization
        """

        import numpy as np
        import sys, petsc4py
        petsc4py.init(sys.argv)
        from petsc4py import PETSc
        
        # Define the constant part of the equation.
        # For the implicit LW scheme this could either zero or the solution at 
        # the current time level (q^n). In this case we set it equal to the 
        # solution at the current time level.
        self.bVec.setArray(state.q)

        # NOTE: HERE WE SHOULD PROBABLY PUT AN IF STATEMENT (?) THE ALLOWS TO SELECT
        #THE IMPLICIT TIME INTEGRATOR THE USER HAS SELECTED. ALL THE IMPLICIT 
        #TIME INTEGRATORS REQUIRE THE SOLUTION OF A NONLINEAR ALGEBRAIC SYSTEM.
        #THE DIFFERENCE AMONG THEM IS THE BUILDING BLOCKS OF THE NONLINEAR
        #FUNCTION INDICATED BY THE LETTER "nlF". FOR INSTANCE FOR THE BACKWARD 
        #EULER SCHEME, I.E.

        #(q^(n+1) - q^(n))/dt + SD(q^(n+1)) = 0 with sd(q^(n+1)) is a nonlinear
        #term arising from the spatial discretization (MOL).

        #THE NONLINEAR FUNCTION IS:

        #nlF = q^(n+1) + dt * sd(q^(n+1)) 

        #THEREFORE THE NONLINEAR SYSTEM TO BE SOLVED IS:

        #nlF = q^(n)

        #WHERE q^(n) DEFINES THE CONSTANT PART OF THE SYSTEM.

        #CAN WE DO BETTER? I MEAN COULD WE DEFINE THE TS INDEPENDENTLY OF THE 
        #SPATIAL DIMENSION OF THE PROBLEM (LOOK AT evalNonLinearFunction)? 
        #YES, WE COULD! BUT IT WOULD REQUIRE TO CHANGE evalNonLinearFunction



        #  Register the function in charge of computing the nonlinear residual
        self.snes.setFunction(self.evalNonLinearFunctionBE, self.fVec)

        # Configure the nonlinear solver to use a matrix-free Jacobian
        self.snes.setUseMF(True)
        #self.snes.getKSP().setType('cg')
        self.snes.setFromOptions()

        self.snes.appctx=(state)

        # Solve the nonlinear problem
        self.snes.solve(self.bVec, state.gqVec)


    def evalNonLinearFunctionBE(self,snes,qin,F):
        r"""
        Computes the nonlinear function for the backward Euler scheme.


        :Input:
         - *qin* - Current approximation of the solution at the next time level,
         i.e. solution of the previous nonlinear solver's iteration.
        """
        
        import numpy as np
        from numpy import zeros, reshape, empty

        state = snes.appctx

        mx = state.grid.ng[0]
        dx = state.grid.d[0]
        mbc = self.mbc
        dt = self.dt

        if state.maux>0:
            state.aux = self.auxbc(state)
        else:
            aux=np.empty((state.maux,mx+2*mbc), order='F')

        dtdx = np.zeros((mx+2*mbc))

        #Have to do this because of how qbc works...
        state.q = reshape(qin,(state.meqn,mx),order='F') 
        qapprox = self.qbc(state)

        from sharpclaw1 import flux1
        ixy = 1
        sd,self.cfl=flux1(qapprox,aux,dt,state.t,ixy,mx,mbc,mx)


        assert sd.flags['F_CONTIGUOUS']
        F.setArray(qapprox[:,mbc:-mbc]-sd[:,mbc:-mbc])













