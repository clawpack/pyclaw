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

from pyclaw.solver import Solver

import limiters.tvd
import riemann

# ========================================================================
#  User-defined routines
# ========================================================================
def start_step(solver,solutions):
    r"""
    Dummy routine called before each step
    
    Replace this routine if you want to do something before each time step.
    """
    pass

def src(solver,solutions,t,dt):
    r"""
    Dummy routine called to calculate a source term
    
    Replace this routine if you want to include a source term.
    """
    pass

# ============================================================================
#  Generic implicit Clawpack solver class
# ============================================================================
class ImplicitClawSolver(Solver):
    r"""
    Generic implicit Clawpack solver
    
    All implicit Clawpack solvers inherit from this base class.
    
    .. attribute:: mthlim 
    
        Limiter to be used on each wave.  ``Default = limiters.tvd.minmod``
    
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
        for attr in ['limiters','order','fwave','src','start_step']:
            self._required_attrs.append(attr)
        
        # Default required attributes
        self._default_attr_values['mbc'] = 2
        self._default_attr_values['limiters'] = limiters.tvd.minmod
        self._default_attr_values['order'] = 2
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
    b.set(0)    pass 
    
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

        # Call b4step, pyclaw should be subclassed if this is needed
        self.start_step(self,solutions)


        # Get state object
        state = solutions['n'].states[0]

        # Compute slution at the new time level
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
#  Implicit ClawPack 1d Solver Class
# ============================================================================
class ImplicitClawSolver1D(ImplicitClawSolver):
    r"""
    Implicit Clawpack evolution routine in 1D
    
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
        Perform essential solver setup.  This routine must be called before
        solver.step() may be called.
        """
        self.set_mthlim()
        if(self.kernel_language == 'Fortran'):
            self.set_fortran_parameters(solutions)


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

        # Get grid object
        grid = state.grid

        # Get dimensionality of the problem
        ndim = grid.ndim

        # Get number of equations
        meqn = state.meqn

        # Get number of ghost cells
        mbc = self.mbc

        # Number of cells
        mx = state.ng[0]

        # Length of the 1D arrays that must be passed to PETSc's SNES
        vl = mx*meqn+2*mbc
            
        # Create PETSc nonlinear solver
        snes = PETSc.SNES()
        snes.create(PETSc.COMM_SELF)

        # Define the vector in charge of containing the solution of the nonlinear system
        x = PETSc.Vec().createSeq(vl)

        # Define the function in charge of computing the nonlinear residual
        f = PETSc.Vec().createSeq(vl)

        # Define the constant part of the equation.
        # For the implicit LW scheme this could either zero or the solution at the current time level (q^n).
        # In this case we set it equal to the solution at the current time level q^(n)
        b = PETSc.Vec().createSeq(vl)
        b.set(reshape(state.qbc,(mx*meqn+2*mbc,1))


        #  Register the function in charge of computing the nonlinear residual
        self.snes.setFunction(???, f) #TODO!!!!!
         
        # Initial guess x = q^n, i.e. solution at the current time level
        x = reshape(state.qbc,(mx*meqn+2*mbc,1))

        # Configure the nonlinear solver to use a matrix-free Jacobian
        snes.setUseMF(True)
        snes.getKSP().setType('cg')
        snes.setFromOptions()


        # Solve the nonlinear problem
        snes.solve(b, x)




    def homogeneous_part(self,solutions, qnext):
        r"""
        Take one time step on the homogeneous hyperbolic system and return the 
        contribution to the nonlinear function.

        Takes one time step of size dt on the hyperbolic system defined in the
        appropriate Riemann solver rp.

        :Input:
         - *solutions* - (:class:`~pyclaw.solution.Solution`) Solution that 
           will be evolved

        """
        import numpy as np

        state = solutions['n'].states[0]
        grid = state.grid
            
        meqn,maux,mwaves,mbc = state.meqn,state.maux,self.mwaves,self.mbc
          
        if maux>0:
            aux = self.auxbc(state)

        ######################################################
        # IMPORT classicimplicit1 
        ######################################################
        if(self.kernel_language == 'Fortran'):
            import classicimplicit1 as classic1

            mx = grid.ng[0]
            dx,dt = grid.d[0],self.dt
            dtdx = np.zeros( (mx+2*mbc) ) + dt/dx
            
            if(maux == 0): aux = np.empty( (maux,mx+2*mbc) )
        
       
            f,self.cfl = classic1.homodisc1(mx,mbc,mx,qnext,aux,dx,dt,self.method,self.mthlim)

            ##################################################################################
            # NOTE:  f is a multidimensional array and not a 1D array.
            # A 1D array is the correct entity that must be used for the nonlinear function. 
            # Consequently, f must be manipulated to push correctly the 
            # contribution of the hyperbolic part to the nonlinear function: f ---> fun
            ##################################################################################
            fun = np.zeros( (mx*meqn) )

            fun = reshape(f,(mx*meqn+2*mbc,1))
            #for icell in xrange(0,mx):
            #        fun[icell*meqn:(icell+1)*meqn] = f[:,icell]


        elif(self.kernel_language == 'Python'):
            raise NotImplementedError('kernel_language = Python not implemented!')

        else: raise Exception("Unrecognized kernel_language; choose 'Fortran' or 'Python'")


        return fun[:]



 



 




