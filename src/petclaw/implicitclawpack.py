r"""
Module containing the implicit classic Clawpack solvers.

This module contains the implicit clawpack solvers (implicit Lax-Wendroff scheme).  
All implicit clawpack solvers inherit from the :class:`ImplicitClawSolver` superclass which in turn 
inherits from the :class:`~pyclaw.solver.Solver` superclass. These
are both pure virtual classes; the only solver classes that should be instantiated
are the dimension-specific ones, :class:`ImplicitClawSolver1D` and :class:`ImplicitClawSolver2D`.

:Authors:
    Matteo Parsani -- Initial version (July 2011)
"""

import petclaw.solver

import pyclaw.limiters.tvd
import riemann

from pyclaw.clawpack import start_step, src

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
    # ========== Generic Init Routine ========================================
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
    
    # ========== Setup Routine ===============================================
    def setup(self,solutions):
        r"""
        Called before any set of time steps.
        
        This routine will be called once before the solver is used via the
        :class:`~pyclaw.controller.Controller`.
        """
    pass 
    
    # ========== Riemann solver library routines =============================   
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
    
    def set_riemann_solver(self,solver_name):
        r"""
        Assigns the library solver solver_name as the Riemann solver.
        
        :Input:
         - *solver_name* - (string) Name of the solver to be used, raises a 
           NameError if the solver does not exist.
        """
        raise Exception("Cannot set a Riemann solver with this class," +
                                        " use one of the derived classes.")
         
    # ========== Time stepping routines ======================================
    def step(self,solutions):
        r"""
        Evolve q of one time step
        

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
        # TIME STEPPING BECAUSE, FOR INSTANCE, THE AUX ARRAY COULD DEPENDS ON THE SOLUTION ITSELF.
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
        if self.cfl >= self.cfl_max:
            raise CFLError('cfl_max exceeded')

        return True


    def updatesolution(self,state):
        r"""
        Compute slution at the new time level for the implicit
        Lax-Wendroff scheme.
        
        This is a dummy routine and must be overridden.
        """
        raise Exception("Dummy routine, please override!")   


    def set_mthlim(self):
        self.mthlim = self.limiters
        if not isinstance(self.limiters,list): self.mthlim=[self.mthlim]
        if len(self.mthlim)==1: self.mthlim = self.mthlim * self.mwaves
        if len(self.mthlim)!=self.mwaves:
            raise Exception('Length of solver.limiters is not equal to 1 or to solver.mwaves')



# ============================================================================
#  Application context (appc) for PETSc SNES; 1D implicit Lax-Wendroff
# ============================================================================
class ImplicitLW1D:
    r"""
    Application context for the nonlinear problem arising by the implicit 
    Lax-Wendroff scheme. It contains some parametes and knows how to compute the 
    nonlinear function.     
    """
    def __init__(self,state,mwaves,mbc,method,mthlim,dt,aux,src,kernel_language='Fortran'):
        self.state = state  
        self.mwaves = mwaves
        self.mbc = mbc
        self.method = method
        self.mthlim = mthlim
        self.dt = dt
        self.aux = aux
        self.kernel_language = kernel_language

        self.grid = self.state.grid
        self.meqn = self.state.meqn
        self.maux = self.state.maux

        self.mx = self.grid.ng[0]
        self.dx = self.grid.d[0]

        self.src = src

    # ========== Evaluation of the nonlinear function ==========================  
    def evalNonLinearFunction(self,snes,qin,F):
        r"""
        Computes the contribution of the spatial-temporal discretization and the 
        source term to the nonlinear function.

        :Input:
         - *qin* - Current approximation of the solution at the next time level,
         i.e. solution of the previous nonlinear solver's iteration.
        """
        from numpy import zeros, reshape

        mx = self.mx
        mbc = self.mbc
        aux = self.aux
        dx = self.dx
        dt = self.dt
        method = self.method
        mthlim = self.mthlim
        meqn = self.meqn

        import classicimplicit1 as classic1

        # Compute the contribution of the homegeneous PDE to the nonlinear 
        # function
        ###################################################################
        dtdx = zeros((mx+2*mbc)) + dt/dx
        
        # Reshape array X before passing it to the fortran code which works 
        # with multidimensional array
        qapprox = reshape(qin,(meqn,mx+2*mbc),order='F')
        fhomo,self.cfl = classic1.homodisc1(mx,mbc,mx,qapprox,aux,dx,dt,method,mthlim)

        # Compute the contribution of the source term to the nonlinear 
        # function
        #fsrc = self.src(self.state,qapprox,self.state.t)

        # Sum the two contribution without creating an additional array
        #fhomo += fsrc

        assert fhomo.flags['F_CONTIGUOUS']
        F.setArray(fhomo)
        

# ============================================================================
#  Implicit ClawPack 1d Solver Class
# ============================================================================
class ImplicitClawSolver1D(ImplicitClawSolver):
    r"""
    Implicit Clawpack evolution routine in 1D.
    
    This class represents the 1d implicit clawpack solver on a single grid.  Note that 
    there are routines here for interfacing with the fortran time stepping 
    routines and the python time stepping routines.  The ones used are 
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
        r"""
        Create 1d implicit Clawpack solver
        
        See :class:`ImplicitClawSolver1D` for more info.
        """   
        
        # Add the functions as required attributes
        self._required_attrs.append('rp')
        self._default_attr_values['rp'] = None
        
        # Import Riemann solvers
        exec('import riemann',globals())
            
        self.ndim = 1

        super(ImplicitClawSolver1D,self).__init__(data)


    # ========== Setup routine =============================   
    def setup(self,solutions):
        r"""
        Perform essential solver setup. This routine must be called before
        solver.step() may be called.

        Set Fortran data structures (for Clawpack) and set up a DA with
        the appropriate stencil width.
        """
        # This is a hack to deal with the fact that petsc4py
        # doesn't allow us to change the stencil_width (mbc)
        state = solutions['n'].state
        state.set_stencil_width(self.mbc)
        # End hack

        self.set_mthlim()
        if(self.kernel_language == 'Fortran'):
            self.set_fortran_parameters(solutions)

        self.bVec = state.lqVec.duplicate()
        self.qnewVec = state.lqVec.duplicate()


    def set_fortran_parameters(self,solutions):
        r"""
        Pack parameters into format recognized by implicit Clawpack (Fortran) code.

        Sets the method array and the cparam common block for the Riemann solver.
        """
        import numpy as np

        state = solutions['n'].state

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
 
        ######################################################
        # IMPORT classicimplicit1 
        #####################################################
        import classicimplicit1 as classic1

        state.set_cparam(classic1)


    def teardown(self):
        r"""
        Delete Fortran objects, which otherwise tend to persist in Python sessions.
        """
        if(self.kernel_language == 'Fortran'):
            import classicimplicit1 as classic1
            del classic1


    # ========== Riemann solver library routines =============================   
    def list_riemann_solvers(self):
        r"""
        List available Riemann solvers 
        
        This routine returns a list of available Riemann solvers which is
        constructed in the Riemann solver package (_pyclaw_rp).  In this case
        it lists only the 1D Riemann solvers.
        
        :Output:
         - (list) - List of Riemann solver names valid to be used with
           :meth:`set_riemann_solver`
        
        .. note::
            These Riemann solvers are currently only accessible to the python 
            time stepping routines.
        """
        return riemann.rp_solver_list_1d
    
    def set_riemann_solver(self,solver_name):
        r"""
        Assigns the library solver solver_name as the Riemann solver.
        
        :Input:
         - *solver_name* - (string) Name of the solver to be used, raises a 
           ``NameError`` if the solver does not exist.
        """
        import logging
        if solver_name in riemann.rp_solver_list_1d:
            exec("self.rp = riemann.rp_%s_1d" % solver_name)
        else:
            logger = logging.getLogger('solver')
            error_msg = 'Could not find Riemann solver with name %s -- perhaps you need to add the solver to riemann.__init__.py.' % solver_name
            logger.warning(error_msg)
            raise NameError(error_msg)


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
        
        
        # Define aux array
        maux = state.maux
        mx = state.grid.ng[0]
        mwaves,mbc = self.mwaves,self.mbc

        if maux>0:
            aux = self.auxbc(state)
        else:
            aux=np.empty((maux,mx+2*mbc))

        # Create application context (appc) and PETSc nonlinear solver
        appc = ImplicitLW1D(state,mwaves,mbc,self.method,self.mthlim,self.dt,aux,self.src,self.kernel_language)
        snes = PETSc.SNES().create()
        
        # Define the vector in charge of containing the solution of the 
        # nonlinear system. The initial guess is qnew = q^n, i.e. solution at 
        # the current time level t^n. 
        self.qnewVec.setArray(state.qbc)
                
        # Define the function in charge of computing the nonlinear residual.
        f = PETSc.Vec().createSeq(state.qbc.size)

        # Define the constant part of the equation.
        # For the implicit LW scheme this could either zero or the solution at 
        # the current time level (q^n). In this case we set it equal to the 
        # solution at the current time level.
        self.bVec.setArray(state.qbc)

        #  Register the function in charge of computing the nonlinear residual
        snes.setFunction(appc.evalNonLinearFunction, f)
         
        # Configure the nonlinear solver to use a matrix-free Jacobian
        snes.setUseMF(True)
        snes.getKSP().setType('cg')
        snes.setFromOptions()

        # Solve the nonlinear problem
        snes.solve(self.bVec, self.qnewVec)

        # Assign to q the new value qnew.
        state.qbc = self.qnewVec.getArray()
