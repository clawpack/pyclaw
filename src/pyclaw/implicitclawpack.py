r"""
Module containing the classic Clawpack solvers.

This module contains the pure and wrapped implicit clawpack solvers (implicit Lax-Wendroff).  
All clawpack solvers inherit from the :class:`ImplicitClawSolver` superclass which in turn 
inherits from the :class:`~pyclaw.solver.Solver` superclass.  These
are both pure virtual classes; the only solver classes that should be instantiated
are the dimension-specific ones, :class:`ImplicitClawSolver1D` and :class:`ImplicitClawSolver2D`.

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
     - (:class:`ImplicitClawSolver`) - Initialized clawpack solver
    
    :Version: 1.0 (2009-06-01)
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
        Evolve solutions one time step

        This routine encodes the generic order in a full time step in this
        order:
        
        1. The :meth:`start_step` function is called
        
        2. The hyperbolic part of the problem (divergence of the flux function) 
           is discretized,i.e. :math:`f(q)_x` in :math:`q_t + f(q)_x = 0`. 
           This is the first term that contributes to the nonlinear function.
        
        3. The source term, if any, is discretized using :fun:`src`
        

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


        # Compute to the contribution of the hyperbolic part to the nonlinear function  
        nonlinearfun = self.implicitLW(solutions)

        # HERE WE SHOULD CALL THE NONLINEAR SOLVER

 
        return True
            

    def set_mthlim(self):
        self.mthlim = self.limiters
        if not isinstance(self.limiters,list): self.mthlim=[self.mthlim]
        if len(self.mthlim)==1: self.mthlim = self.mthlim * self.mwaves
        if len(self.mthlim)!=self.mwaves:
            raise Exception('Length of solver.limiters is not equal to 1 or to solver.mwaves')


    def implicitLW(self,state):
        """
        Evaluate dq/dt * (delta t)
        """
        
        # Compute the contribution of the hyperbolic term (divergence of the flux) 
        fun = self.dq_homogeneous(state)

        # Check here if we violated the CFL condition, if we did, return 
        # immediately to evolve_to_time and let it deal with picking a new
        # dt. 
        #
        # Even for steady state calculation the control of the CFL is important, especially
        # in the first few pseudo time steps, where high frequency errors components must
        # be damped out and expelled from the computational domain. 
        # TODO: implement a CFL-law
        self.communicateCFL()
        if self.cfl >= self.cfl_max:
            raise CFLError('cfl_max exceeded')

        # Compute the contribution of the source term, if any
        if self.src_term == 1:
            fun+=self.src(state,q,state.t)

        return fun


    def dq_homogeneous(state):
        raise NotImplementedError('You must subclass ImplicitClawSolver.')

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
    

 

