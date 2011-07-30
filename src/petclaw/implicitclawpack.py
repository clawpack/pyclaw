r"""
Module containing the implicit classic Clawpack solvers.

This module contains the implicit clawpack solvers, i.e. implicit Lax-Wendroff
scheme.

All implicit clawpack solvers inherit from the :class:`ImplicitClawSolver` 
superclass which in turn inherits from the :class:`~pyclaw.solver.Solver` 
superclass. These are both pure virtual classes; the only solver classes that 
should be instantiated are the dimension-specific ones, 
:class:`ImplicitClawSolver1D` and :class:`ImplicitClawSolver2D`.

:Authors:
    Matteo Parsani  -- July 2011
    David Ketcheson -- July 2011
"""

# Import modules
import petclaw.solver
import pyclaw.limiters.tvd
from pyclaw.solver import CFLError
import riemann


# Define some default functions
def start_step(solver,solutions):
    r"""
    Dummy routine called before each step.
    
    Replace this routine if you want to do something before each time step.
    """
    pass

def src(solver,grid,q,t):
    r"""
    Dummy routine called to calculate a source term.
    
    Replace this routine if you want to include a source term.
    """
    pass

# ============================================================================
#  Generic implicit Clawpack solver class
# ============================================================================
class ImplicitClawSolver(petclaw.solver.PetSolver):
    r"""
    Generic implicit Clawpack solver.
    
    All implicit Clawpack solvers inherit from this base class.
    
    .. attribute:: mthlim 
    
        Limiter to be used on each wave.  ``Default = pyclaw.limiters.tvd.minmod``
    
    .. attribute:: order
    
        Order of the solver, either 1 for first order or 2 for second order 
        corrections.  ``Default = 2``
    
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
     - (:class:`ImplicitClawSolver`) - Initialized implicit clawpack solver
    
    """

    def __init__(self,data=None):
        r"""
        See :class:`ImplicitClawSolver` for full documentation.
        """
        
        # Required attributes for this solver
        for attr in ['limiters','order','src_split','fwave','src','start_step']:
            self._required_attrs.append(attr)
        
        # Default required attributes
        self._default_attr_values['mbc'] = 2
        self._default_attr_values['limiters'] = pyclaw.limiters.tvd.minmod
        self._default_attr_values['order'] = 2
        self._default_attr_values['src_split'] = 0
        self._default_attr_values['fwave'] = False
        self._default_attr_values['src'] = src
        self._default_attr_values['start_step'] = start_step
        self._default_attr_values['kernel_language'] = 'Fortran'
        self._default_attr_values['verbosity'] = 0

        # Call general initialization function
        super(ImplicitClawSolver,self).__init__(data)

    
    def initiate(self,solutions):
        r"""
        Called before any set of time steps.
        
        This routine will be called once before the solver is used via the
        :class:`~pyclaw.controller.Controller`.
        """

        # Import modules
        from petsc4py import PETSc
        from numpy import empty
        from pyclaw.state import State

        # Call b4step
        self.start_step(self,solutions)

        # Get state
        state = solutions['n'].state
    
        # Set up a DA with the appropriate stencil width.
        state.set_stencil_width(self.mbc)

        # Set mthlim
        self.set_mthlim()

        # Create PETSc vectors in charge of containig:
        # bVec: the constant part of the nonlinear algebraic system of equations
        # fVec: nonlinear vector-valued function
        self.bVec = state.gqVec.duplicate()
        self.fVec = state.gqVec.duplicate()

        # Create PETSc nonlinear solver
        self.snes = PETSc.SNES().create()

        # Ought to implement a copy constructor for State
        self.impsol_stage = State(state.grid)
        self.impsol_stage.meqn             = state.meqn
        self.impsol_stage.maux             = state.maux
        self.impsol_stage.aux_global       = state.aux_global
        self.impsol_stage.t                = state.t
        if state.maux > 0:
            self.impsol_stage.aux          = state.aux


    def step(self,solutions):
        r"""
        Evolve q for one time step.

        This routine is called from the method evolve_to_time defined in the
        pyclaw.solver.Solver superclass.

        :Input:
         - *solutions* - (:class:`~pyclaw.solution.Solution`) Dictionary of 
           solutions to be evolved.
         
        :Output: 
         - (bool) - True if full step succeeded, False otherwise.
        """

        # Import library
        from pyclaw.solution import Solution

        # Get state
        state = solutions['n'].states[0]
        
        # Call b4step, pyclaw should be subclassed if this is needed
        # NOTE: THIS IS NOT THE RIGHT PLACE WHEN WE HAVE AN IMPLICIT  
        # TIME STEPPING BECAUSE, FOR INSTANCE, THE AUX ARRAY COULD DEPEND ON 
        # THE SOLUTION ITSELF.
        # THEN start_step FUNCTION SHOULD BE PLACED IN THE PART OF THE CODE 
        # WHERE THE NONLINEAR FUNCTION IS COMPUTED!!!!!!!!!
        self.start_step(self,solutions)

        ########################################
        # Compute solution at the new time level
        ########################################
        # Import libraries
        import sys, petsc4py
        from petsc4py import PETSc
        petsc4py.init(sys.argv)

        # Set the constant part of the equation and register the function in 
        # charge of computing the nonlinear vector-valued function
        self.set_bVec(state)
        self.snes.setFunction(self.nonlinearfunction, self.fVec)

        # Configure the nonlinear solver to use a matrix-free Jacobian
        self.snes.setUseMF(True)
        self.snes.setFromOptions()

        # Pass additinal properties to SNES.
        qold = state.q.copy()
        self.snes.appctx=(state,qold)

        # Solve the nonlinear problem
        self.snes.solve(self.bVec, state.gqVec)

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


    def set_bVec(self,state):
        r"""
        Set the constant part of the nonlinear algebraic system arising from the
        implicit Lax-Wendroff discretization. 

        :Input:
         - *state* -

        """

        # Set the constant part of the nonlinear algebraic system equal to zero.
        self.bVec.setArray(0)


    def set_mthlim(self):
        self.mthlim = self.limiters
        if not isinstance(self.limiters,list): self.mthlim=[self.mthlim]
        if len(self.mthlim)==1: self.mthlim = self.mthlim * self.mwaves
        if len(self.mthlim)!=self.mwaves:
            raise Exception('Length of solver.limiters is not equal to 1 or to solver.mwaves')


# ============================================================================
#  Implicit ClawPack 1d Solver Class
# ============================================================================
class ImplicitClawSolver1D(ImplicitClawSolver):
    r"""
    Implicit Clawpack evolution routine in 1D.
    
    This class represents the 1d implicit clawpack solver on a single grid.  
    Note that there are routines here for interfacing with the fortran time 
    stepping routines and the python time stepping routines.  The ones used are 
    dependent on the argument given to the initialization of the solver 
    (defaults to python).
    
    .. attribute:: rp
    
        Riemann solver function.
        
    :Initialization:
    
    Input:
     - *data* - (:class:`~pyclaw.data.Data`) An instance of a Data object whose
       parameters can be used to initialize this solver
    Output:
     - (:class:`ImplicitClawSolver1D`) - Initialized 1d implicit clawpack solver
        
    """

    def __init__(self,data=None):

        # Set physical dimensions
        self.ndim = 1

        # Call superclass __init__
        super(ImplicitClawSolver1D,self).__init__(data)

  
    def setup(self,solutions):
        r"""
        Perform essential solver setup. This routine must be called before
        solver.step() may be called.

        Set Fortran data structures (for Clawpack).
        """

        # Call parent's "setup" function
        self.initiate(solutions)

        # Set Fortran data structure for the 1D implicit ClawPack solver
        if(self.kernel_language == 'Fortran'):
            self.set_fortran_parameters(solutions)


    def set_fortran_parameters(self,solutions):
        r"""
        Pack parameters into format recognized by implicit Clawpack (Fortran) 
        code.

        Sets the method array and the cparam common block for the Riemann 
        solver.
        """
        
        # Import module
        import numpy as np

        # Get state
        state = solutions['n'].state

        # Set parameters that will be use in the Fortran routine
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

        import classicimplicit1 as classic1
        state.set_cparam(classic1)


    def teardown(self):
        r"""
        Delete Fortran objects, which otherwise tend to persist in Python 
        sessions.
        """

        if(self.kernel_language == 'Fortran'):
            import classicimplicit1 as classic1
            del classic1


    def nonlinearfunction(self,snes,qin,nlf):
        r"""
        Computes the contribution of the spatial-temporal discretization and the 
        source term to the nonlinear function.

        :Input:
         - *qin* - Current approximation of the solution at the next time level,
         i.e. solution of the previous nonlinear solver's iteration.
        """

        # Import module
        import numpy as np
        from numpy import zeros, reshape, empty

        # Get state
        state,qold = snes.appctx

        # Get some parameters that will be used in the Fortran routine 
        mx = state.grid.ng[0]
        dx = state.grid.d[0]
        mbc = self.mbc
        dt = self.dt

        # Auxbc is set here and not outside of this function because it is 
        # more general. Indeed aux could depend on q which changes at each 
        # nonlinear iteration!
        if state.maux>0:
            state.aux = self.auxbc(state)
        else:
            aux = empty((state.maux,mx+2*mbc))

        # Define and set to zero the ratio between dt and dx 
        dtdx = zeros((mx+2*mbc)) + dt/dx

        #Have to do this because of how qbc works...
        state.q = reshape(qin,(state.meqn,mx),order='F') 
        qapprox = self.qbc(state)

        # Import module
        import classicimplicit1 as classic1

        # Call Fortran routine
        dq,self.cfl = classic1.step1imp(mx,mbc,mx,qapprox,aux,dx,dt,self.method,self.mthlim)

        # Compute the nonlinear vector-valued function
        assert dq.flags['F_CONTIGUOUS']
        nlf.setArray(qapprox[:,mbc:-mbc]-dq[:,mbc:-mbc]-qold[:,:])

