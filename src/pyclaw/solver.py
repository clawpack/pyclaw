r"""
Module specifying the interface to every solver in PyClaw.
"""
import logging

class CFLError(Exception):
    """Error raised when cfl_max is exceeded.  Is this a
       reasonable mechanism for handling that?"""
    def __init__(self,msg):
        super(CFLError,self).__init__(msg)

class BC():
    """Enumeration of boundary condition names.
       This could instead just be a static dictionary."""
    custom     = 0
    outflow    = 1
    periodic   = 2
    reflecting = 3

#################### Dummy routines ######################
def default_compute_gauge_values(q,aux):
    r"""By default, record values of q at gauges.
    """
    return q

class Solver(object):
    r"""
    Pyclaw solver superclass.

    The pyclaw.Solver.solver class is an abstract class that should
    not be instantiated; rather, all Solver classes should inherit from it.

    A Solver is typically instantiated as follows::

        >>> solver = pyclaw.ClawSolver2d()

    After which solver options may be set.  It is always necessary to set
    solver.mwaves to the number of waves used in the Riemann solver.
    Typically it is also necessary to set the boundary conditions (for q, and
    for aux if an aux array is used).  Many other options may be set
    for specific solvers; for instance the limiter to be used, whether to
    use a dimensionally-split algorithm, and so forth.

    Usually the solver is attached to a controller before being used::

        >>> claw = pyclaw.Controller()
        >>> claw.solver = solver

    .. attribute:: dt
        
        Current time step, ``default = 0.1``
        
    .. attribute:: cfl
        
        Current Courant-Freidrichs-Lewy number, ``default = 1.0``
    
    .. attribute:: status
        
        Dictionary of status values for the solver with the following keys:
         - ``cflmax`` = Maximum CFL number
         - ``dtmin`` = Minimum time step taken
         - ``dtmax`` = Maximum time step taken
         - ``numsteps`` = Current number of time steps that have been taken

        solver.status is reset each time solver.evolve_to_time is called, and
        it is also returned by solver.evolve_to_time.
    
    .. attribute:: dt_variable
    
        Whether to allow the time step to vary, ``default = True``.
        If false, the initial time step size is used for all steps.
        
    .. attribute:: max_steps
    
        The maximum number of time steps allowd to reach the end time 
        requested, ``default = 1000``.  If exceeded, an exception is
        raised.
    
    .. attribute:: logger
    
        Default logger for all solvers.  Records information about the run
        and debugging messages (if requested).

    .. attribute:: bc_lower 
    
        (list of ints) Lower boundary condition types, listed in the
        same order as the Dimensions of the Grid.  See Solver.BC for
        an enumeration.

    .. attribute:: bc_upper 
    
        (list of ints) Upper boundary condition types, listed in the
        same order as the Dimensions of the Grid.  See Solver.BC for
        an enumeration.

    .. attribute:: user_bc_lower 
        
        (func) User defined lower boundary condition.
        Fills the values of qbc with the correct boundary values.
        The appropriate signature is:

        def user_bc_lower(grid,dim,t,qbc,mbc):

    .. attribute:: user_bc_upper 
    
        (func) User defined upper boundary condition.
        Fills the values of qbc with the correct boundary values.
        The appropriate signature is:

        def user_bc_upper(grid,dim,t,qbc,mbc):
 
        
    :Initialization:
    
    Output:
     - (:class:`Solver`) - Initialized Solver object
    """

            
    def __setattr__(self, key, value):
        if not hasattr(self, '_isinitialized'):
            self.__dict__['_isinitialized'] = False
        if self._isinitialized and not hasattr(self, key):
            raise TypeError("%s has no attribute %s" % (self.__class__,key))
        object.__setattr__(self,key,value)

    #  ======================================================================
    #   Initialization routines
    #  ======================================================================
    def __init__(self,claw_package=None):
        r"""
        Initialize a Solver object
        
        See :class:`Solver` for full documentation
        """ 
        # Setup solve logger
        self.logger = logging.getLogger('evolve')

        self.dt_initial = 0.1
        self.dt_max = 1e99
        self.max_steps = 1000
        self.dt_variable = True
        self.mwaves = None #Must be set later to agree with Riemann solver
        self.so_name = None #Can remove this after merging fwaves commit
        self.qbc = None
        self.auxbc = None
        self.rp = None

        # select package to build solver objects from, by default this will be
        # the package that contains the module implementing the derived class
        # for example, if ClawSolver1D is implemented in 'petclaw.solver', then
        # the computed claw_package will be 'petclaw'
        
        import sys
        if claw_package is not None and claw_package in sys.modules:
            self.claw_package = sys.modules[claw_package]
        else:
            claw_package_name = self.__module__[0:self.__module__.rfind('.')]
            if claw_package_name in sys.modules:
                self.claw_package = sys.modules[claw_package_name]
            else:
                raise NotImplementedError("Unable to determine solver package, please provide one")

        
        # Initialize time stepper values
        self.dt = self.dt_initial
        self.cfl = self.claw_package.CFL(self.cfl_desired)
       
        # Status Dictionary
        self.status = {'cflmax':self.cfl.get_cached_max(),
                       'dtmin':self.dt, 
                       'dtmax':self.dt,
                       'numsteps':0 }
        
        # No default BCs; user must set them
        self.bc_lower =    [None]*self.ndim
        self.bc_upper =    [None]*self.ndim
        self.aux_bc_lower = [None]*self.ndim
        self.aux_bc_upper = [None]*self.ndim
        
        self.user_bc_lower = None
        self.user_bc_upper = None

        self.user_aux_bc_lower = None
        self.user_aux_bc_upper = None

        self.compute_gauge_values = None
        r"""(function) - Function that computes quantities to be recorded at gaugues"""

        self.qbc          = None
        r""" Array to hold ghost cell values.  This is the one that gets passed
        to the Fortran code.  """

        self._isinitialized = True

        super(Solver,self).__init__()


    # ========================================================================
    #  Solver setup and validation routines
    # ========================================================================
    def is_valid(self):
        r"""
        Checks that all required solver attributes are set.
        
        Checks to make sure that all the required attributes for the solver 
        have been set correctly.  All required attributes that need to be set 
        are contained in the attributes list of the class.
        
        Will post debug level logging message of which required attributes 
        have not been set.
        
        :Output:
         - *valid* - (bool) True if the solver is valid, False otherwise
        
        """
        valid = True
        if any([bcmeth == BC.custom for bcmeth in self.bc_lower]):
            if self.user_bc_lower is None:
                self.logger.debug('Lower custom BC function has not been set.')
                valid = False
        if any([bcmeth == BC.custom for bcmeth in self.bc_upper]):
            if self.user_bc_lower is None:
                self.logger.debug('Upper custom BC function has not been set.')
                valid = False
        return valid
        
    def setup(self,solution):
        r"""
        Stub for solver setup routines.
        
        This function is called before a set of time steps are taken in order 
        to reach tend.  A subclass should extend or override it if it needs to 
        perform some setup based on attributes that would be set after the 
        initialization routine.  Typically this is initialization that
        requires knowledge of the solution object.
        """

        pass

    def teardown(self):
        r"""
        Stub for solver teardown routines.
        
        This function is called at the end of a simulation.
        A subclass should override it only if it needs to 
        perform some cleanup, such as deallocating arrays in a Fortran module.
        """
        pass


    def __str__(self):
        output = "Solver Status:\n"
        for (k,v) in self.status.iteritems():
            output = "\n".join((output,"%s = %s" % (k.rjust(25),v)))
        return output


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
        for i in range(nregisters-1):
            #Maybe should use State.copy() here?
            self._rk_stages.append(State(state.grid,state.meqn,state.maux))
            self._rk_stages[-1].aux_global       = state.aux_global
            self._rk_stages[-1].set_mbc(self.mbc)
            self._rk_stages[-1].t                = state.t
            if state.maux > 0:
                self._rk_stages[-1].aux              = state.aux


    # ========================================================================
    #  Boundary Conditions
    # ========================================================================    
    def allocate_bc_arrays(self,state):
        r"""
        Create numpy arrays for q and aux with ghost cells attached.
        These arrays are referred to throughout the code as qbc and auxbc.

        This is typically called by solver.setup().
        """
        import numpy as np
        qbc_dim = [n+2*self.mbc for n in state.grid.ng]
        qbc_dim.insert(0,state.meqn)
        self.qbc = np.zeros(qbc_dim,order='F')

        auxbc_dim = [n+2*self.mbc for n in state.grid.ng]
        auxbc_dim.insert(0,state.maux)
        self.auxbc = np.empty(auxbc_dim,order='F')
        if state.maux>0:
            self.apply_aux_bcs(state)

    def apply_q_bcs(self,state):
        r"""
        Fills in solver.qbc (the local vector), including ghost cell values.
    
        This function returns an array of dimension determined by the 
        :attr:`mbc` attribute.  The type of boundary condition set is 
        determined by :attr:`bc_lower` and :attr:`bc_upper` for the 
        approprate dimension.  Valid values for :attr:`bc_lower` and 
        :attr:`bc_upper` include:
        
        - 'custom'     or 0: A user defined boundary condition will be used, the appropriate 
            Dimension method user_bc_lower or user_bc_upper will be called.
        - 'outflow'    or 1: Zero-order extrapolation.
        - 'periodic'   or 2: Periodic boundary conditions.
        - 'reflecting' or 3: Wall boundary conditions. It is assumed that the second 
            component of q represents velocity or momentum.
    
        :Input:
         -  *grid* - (:class:`Grid`) The grid being operated on.
         -  *state* - The state being operated on; this may or may not be the
                      same as *grid*.  Generally it is the same as *grid* for
                      the classic algorithms and other one-level algorithms, 
                      but different for method-of-lines algorithms like SharpClaw.

        :Output:
         - (ndarray(meqn,...)) q array with boundary ghost cells added and set
         

        .. note:: 

            Note that for user-defined boundary conditions, the array sent to
            the boundary condition has not been rolled. 
        """
        
        import numpy as np

        self.qbc = state.get_qbc_from_q(self.mbc,'q',self.qbc)
        grid = state.grid
       
        for idim,dim in enumerate(grid.dimensions):
            # First check if we are actually on the boundary
            # (in case of a parallel run)
            if dim.nstart == 0:
                # If a user defined boundary condition is being used, send it on,
                # otherwise roll the axis to front position and operate on it
                if self.bc_lower[idim] == BC.custom:
                    self.qbc_lower(state,dim,state.t,self.qbc,idim)
                elif self.bc_lower[idim] == BC.periodic:
                    if dim.nend == dim.n:
                        # This process owns the whole grid
                        self.qbc_lower(state,dim,state.t,np.rollaxis(self.qbc,idim+1,1),idim)
                    else:
                        pass #Handled automatically by PETSc
                else:
                    self.qbc_lower(state,dim,state.t,np.rollaxis(self.qbc,idim+1,1),idim)

            if dim.nend == dim.n :
                if self.bc_upper[idim] == BC.custom:
                    self.qbc_upper(state,dim,state.t,self.qbc,idim)
                elif self.bc_upper[idim] == BC.periodic:
                    if dim.nstart == 0:
                        # This process owns the whole grid
                        self.qbc_upper(state,dim,state.t,np.rollaxis(self.qbc,idim+1,1),idim)
                    else:
                        pass #Handled automatically by PETSc
                else:
                    self.qbc_upper(state,dim,state.t,np.rollaxis(self.qbc,idim+1,1),idim)


    def qbc_lower(self,state,dim,t,qbc,idim):
        r"""
        Apply lower boundary conditions to qbc
        
        Sets the lower coordinate's ghost cells of *qbc* depending on what 
        :attr:`bc_lower` is.  If :attr:`bc_lower` = 0 then the user 
        boundary condition specified by :attr:`user_bc_lower` is used.  Note 
        that in this case the function :attr:`user_bc_lower` belongs only to 
        this dimension but :attr:`user_bc_lower` could set all user boundary 
        conditions at once with the appropriate calling sequence.
        
        :Input:
         - *grid* - (:class:`Grid`) Grid that the dimension belongs to
         
        :Input/Ouput:
         - *qbc* - (ndarray(...,meqn)) Array with added ghost cells which will
           be set in this routines
        """
        if self.bc_lower[idim] == BC.custom: 
            self.user_bc_lower(state,dim,t,qbc,self.mbc)
        elif self.bc_lower[idim] == BC.outflow:
            for i in xrange(self.mbc):
                qbc[:,i,...] = qbc[:,self.mbc,...]
        elif self.bc_lower[idim] == BC.periodic:
            # This process owns the whole grid
            qbc[:,:self.mbc,...] = qbc[:,-2*self.mbc:-self.mbc,...]
        elif self.bc_lower[idim] == BC.reflecting:
            for i in xrange(self.mbc):
                qbc[:,i,...] = qbc[:,2*self.mbc-1-i,...]
                qbc[idim+1,i,...] = -qbc[idim+1,2*self.mbc-1-i,...] # Negate normal velocity
        else:
            raise NotImplementedError("Boundary condition %s not implemented" % self.bc_lower)


    def qbc_upper(self,state,dim,t,qbc,idim):
        r"""
        Apply upper boundary conditions to qbc
        
        Sets the upper coordinate's ghost cells of *qbc* depending on what 
        :attr:`bc_upper` is.  If :attr:`bc_upper` = 0 then the user 
        boundary condition specified by :attr:`user_bc_upper` is used.  Note 
        that in this case the function :attr:`user_bc_upper` belongs only to 
        this dimension but :attr:`user_bc_upper` could set all user boundary 
        conditions at once with the appropriate calling sequence.
        
        :Input:
         - *grid* - (:class:`Grid`) Grid that the dimension belongs to
         
        :Input/Ouput:
         - *qbc* - (ndarray(...,meqn)) Array with added ghost cells which will
           be set in this routines
        """
 
        if self.bc_upper[idim] == BC.custom:
            self.user_bc_upper(state,dim,t,qbc,self.mbc)
        elif self.bc_upper[idim] == BC.outflow:
            for i in xrange(self.mbc):
                qbc[:,-i-1,...] = qbc[:,-self.mbc-1,...] 
        elif self.bc_upper[idim] == BC.periodic:
            # This process owns the whole grid
            qbc[:,-self.mbc:,...] = qbc[:,self.mbc:2*self.mbc,...]
        elif self.bc_upper[idim] == BC.reflecting:
            for i in xrange(self.mbc):
                qbc[:,-i-1,...] = qbc[:,-2*self.mbc+i,...]
                qbc[idim+1,-i-1,...] = -qbc[idim+1,-2*self.mbc+i,...] # Negate normal velocity
        else:
            raise NotImplementedError("Boundary condition %s not implemented" % self.bc_lower)



    def apply_aux_bcs(self,state):
        r"""
        Appends boundary cells to aux and fills them with appropriate values.
    
        This function returns an array of dimension determined by the 
        :attr:`mbc` attribute.  The type of boundary condition set is 
        determined by :attr:`aux_bc_lower` and :attr:`aux_bc_upper` for the 
        approprate dimension.  Valid values for :attr:`aux_bc_lower` and 
        :attr:`aux_bc_upper` include:
        
        - 'custom'     or 0: A user defined boundary condition will be used, the appropriate 
            Dimension method user_aux_bc_lower or user_aux_bc_upper will be called.
        - 'outflow'    or 1: Zero-order extrapolation.
        - 'periodic'   or 2: Periodic boundary conditions.
        - 'reflecting' or 3: Wall boundary conditions. It is assumed that the second 
            component of q represents velocity or momentum.
    
        :Input:
         -  *grid* - (:class:`Grid`) The grid being operated on.
         -  *state* - The state being operated on; this may or may not be the
                      same as *grid*.  Generally it is the same as *grid* for
                      the classic algorithms and other one-level algorithms, 
                      but different for method-of-lines algorithms like SharpClaw.

        :Output:
         - (ndarray(maux,...)) q array with boundary ghost cells added and set
         

        .. note:: 

            Note that for user-defined boundary conditions, the array sent to
            the boundary condition has not been rolled. 
        """
        
        import numpy as np

        self.auxbc = state.get_qbc_from_q(self.mbc,'aux',self.auxbc)

        grid = state.grid
       
        for idim,dim in enumerate(grid.dimensions):
            # First check if we are actually on the boundary
            # (in case of a parallel run)
            if dim.nstart == 0:
                # If a user defined boundary condition is being used, send it on,
                # otherwise roll the axis to front position and operate on it
                if self.aux_bc_lower[idim] == BC.custom:
                    self.auxbc_lower(state,dim,state.t,self.auxbc,idim)
                elif self.aux_bc_lower[idim] == BC.periodic:
                    if dim.nend == dim.n:
                        # This process owns the whole grid
                        self.auxbc_lower(state,dim,state.t,np.rollaxis(self.auxbc,idim+1,1),idim)
                    else:
                        pass #Handled automatically by PETSc
                else:
                    self.auxbc_lower(state,dim,state.t,np.rollaxis(self.auxbc,idim+1,1),idim)

            if dim.nend == dim.n :
                if self.aux_bc_upper[idim] == BC.custom:
                    self.auxbc_upper(state,dim,state.t,self.auxbc,idim)
                elif self.aux_bc_upper[idim] == BC.periodic:
                    if dim.nstart == 0:
                        # This process owns the whole grid
                        self.auxbc_upper(state,dim,state.t,np.rollaxis(self.auxbc,idim+1,1),idim)
                    else:
                        pass #Handled automatically by PETSc
                else:
                    self.auxbc_upper(state,dim,state.t,np.rollaxis(self.auxbc,idim+1,1),idim)


    def auxbc_lower(self,state,dim,t,auxbc,idim):
        r"""
        Apply lower boundary conditions to auxbc
        
        Sets the lower coordinate's ghost cells of *auxbc* depending on what 
        :attr:`aux_bc_lower` is.  If :attr:`aux_bc_lower` = 0 then the user 
        boundary condition specified by :attr:`user_aux_bc_lower` is used.  Note 
        that in this case the function :attr:`user_aux_bc_lower` belongs only to 
        this dimension but :attr:`user_aux_bc_lower` could set all user boundary 
        conditions at once with the appropriate calling sequence.
        
        :Input:
         - *grid* - (:class:`Grid`) Grid that the dimension belongs to
         
        :Input/Ouput:
         - *auxbc* - (ndarray(maux,...)) Array with added ghost cells which will
           be set in this routines
        """
        if self.aux_bc_lower[idim] == BC.custom: 
            self.user_aux_bc_lower(state,dim,t,auxbc,self.mbc)
        elif self.aux_bc_lower[idim] == BC.outflow:
            for i in xrange(self.mbc):
                auxbc[:,i,...] = auxbc[:,self.mbc,...]
        elif self.aux_bc_lower[idim] == BC.periodic:
            # This process owns the whole grid
            auxbc[:,:self.mbc,...] = auxbc[:,-2*self.mbc:-self.mbc,...]
        elif self.aux_bc_lower[idim] == BC.reflecting:
            for i in xrange(self.mbc):
                auxbc[:,i,...] = auxbc[:,2*self.mbc-1-i,...]
        elif self.aux_bc_lower[idim] is None:
            raise Exception("One or more of the aux boundary conditions aux_bc_upper has not been specified.")
        else:
            raise NotImplementedError("Boundary condition %s not implemented" % self.aux_bc_lower)


    def auxbc_upper(self,state,dim,t,auxbc,idim):
        r"""
        Apply upper boundary conditions to auxbc
        
        Sets the upper coordinate's ghost cells of *auxbc* depending on what 
        :attr:`aux_bc_upper` is.  If :attr:`aux_bc_upper` = 0 then the user 
        boundary condition specified by :attr:`user_aux_bc_upper` is used.  Note 
        that in this case the function :attr:`user_aux_bc_upper` belongs only to 
        this dimension but :attr:`user_aux_bc_upper` could set all user boundary 
        conditions at once with the appropriate calling sequence.
        
        :Input:
         - *grid* - (:class:`Grid`) Grid that the dimension belongs to
         
        :Input/Ouput:
         - *auxbc* - (ndarray(maux,...)) Array with added ghost cells which will
           be set in this routines
        """
 
        if self.aux_bc_upper[idim] == BC.custom:
            self.user_aux_bc_upper(state,dim,t,auxbc,self.mbc)
        elif self.aux_bc_upper[idim] == BC.outflow:
            for i in xrange(self.mbc):
                auxbc[:,-i-1,...] = auxbc[:,-self.mbc-1,...] 
        elif self.aux_bc_upper[idim] == BC.periodic:
            # This process owns the whole grid
            auxbc[:,-self.mbc:,...] = auxbc[:,self.mbc:2*self.mbc,...]
        elif self.aux_bc_upper[idim] == BC.reflecting:
            for i in xrange(self.mbc):
                auxbc[:,-i-1,...] = auxbc[:,-2*self.mbc+i,...]
        elif self.aux_bc_lower[idim] is None:
            raise Exception("One or more of the aux boundary conditions aux_bc_lower has not been specified.")
        else:
            raise NotImplementedError("Boundary condition %s not implemented" % self.aux_bc_lower)


    # ========================================================================
    #  Evolution routines
    # ========================================================================
    def evolve_to_time(self,solution,tend=None):
        r"""
        Evolve solution from solution.t to tend.  If tend is not specified,
        take a single step.
        
        This method contains the machinery to evolve the solution object in
        ``solution`` to the requested end time tend if given, or one 
        step if not.          

        :Input:
         - *solution* - (:class:`Solution`) Solution to be evolved
         - *tend* - (float) The end time to evolve to, if not provided then 
           the method will take a single time step.
            
        :Output:
         - (dict) - Returns the status dictionary of the solver
        """
        
        if tend == None:
            take_one_step = True
        else:
            take_one_step = False
            
        # Parameters for time-stepping
        tstart = solution.t

        # Reset status dictionary
        self.status['cflmax'] = self.cfl.get_cached_max()
        self.status['dtmin'] = self.dt
        self.status['dtmax'] = self.dt
        self.status['numsteps'] = 0

        # Setup for the run
        if not self.dt_variable:
            if take_one_step:
                self.max_steps = 1
            else:
                self.max_steps = int((tend - tstart + 1e-10) / self.dt)
                if abs(self.max_steps*self.dt - (tend - tstart)) > 1e-5 * (tend-tstart):
                    raise Exception('dt does not divide (tend-tstart) and dt is fixed!')
        if self.dt_variable == 1 and self.cfl_desired > self.cfl_max:
            raise Exception('Variable time-stepping and desired CFL > maximum CFL')
        if tend <= tstart:
            self.logger.info("Already at or beyond end time: no evolution required.")
            self.max_steps = 0
                
        # Main time-stepping loop
        for n in xrange(self.max_steps):
            
            state = solution.state
            
            # Adjust dt so that we hit tend exactly if we are near tend
            if solution.t + self.dt > tend and tstart < tend and not take_one_step:
                self.dt = tend - solution.t 

            # Keep a backup in case we need to retake a time step
            if self.dt_variable:
                q_backup = state.q.copy('F')
                told = solution.t
            
            self.step(solution)

            # Check to make sure that the Courant number was not too large
            cfl = self.cfl.get_cached_max()
            if cfl <= self.cfl_max:
                # Accept this step
                self.status['cflmax'] = max(cfl, self.status['cflmax'])
                if self.dt_variable==True:
                    solution.t += self.dt 
                else:
                    #Avoid roundoff error if dt_variable=False:
                    solution.t = tstart+(n+1)*self.dt
                # Verbose messaging
                self.logger.debug("Step %i  CFL = %f   dt = %f   t = %f"
                    % (n,cfl,self.dt,solution.t))
                    
                self.write_gauge_values(solution)
                # Increment number of time steps completed
                self.status['numsteps'] += 1
                # See if we are finished yet
                if solution.t >= tend or take_one_step:
                    break
            else:
                # Reject this step
                self.logger.debug("Rejecting time step, CFL number too large")
                if self.dt_variable:
                    state.q = q_backup
                    solution.t = told
                else:
                    # Give up, we cannot adapt, abort
                    self.status['cflmax'] = \
                        max(cfl, self.status['cflmax'])
                    raise Exception('CFL too large, giving up!')
                    
            # Choose new time step
            if self.dt_variable:
                if cfl > 0.0:
                    self.dt = min(self.dt_max,self.dt * self.cfl_desired 
                                    / cfl)
                    self.status['dtmin'] = min(self.dt, self.status['dtmin'])
                    self.status['dtmax'] = max(self.dt, self.status['dtmax'])
                else:
                    self.dt = self.dt_max

      
        # End of main time-stepping loop -------------------------------------

        if self.dt_variable and solution.t < tend \
                and self.status['numsteps'] == self.max_steps:
            raise Exception("Maximum number of timesteps have been taken")

        return self.status

    def step(self):
        r"""
        Take one step
        
        This method is only a stub and should be overridden by all solvers who
        would like to use the default time-stepping in evolve_to_time.
        """
        raise NotImplementedError("No stepping routine has been defined!")

    # ========================================================================
    #  Gauges
    # ========================================================================
    def write_gauge_values(self,solution):
        r"""Write solution (or derived quantity) values at each gauge coordinate
            to file.
        """
        for i,gauge in enumerate(solution.state.grid.gauges):
            x=gauge[0]; y=gauge[1]
            aux=solution.state.aux[:,x,y]
            q=solution.state.q[:,x,y]
            p=self.compute_gauge_values(q,aux)
            t=solution.t
            solution.state.grid.gauge_files[i].write(str(t)+' '+' '.join(str(j) for j in p)+'\n')  

