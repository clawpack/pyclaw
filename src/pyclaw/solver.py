r"""
Module specifying the interface to every solver in PyClaw.
"""
from __future__ import absolute_import
import logging
import numpy as np
import six
from six.moves import range

class CFLError(Exception):
    """Error raised when cfl_max is exceeded.  Is this a
       reasonable mechanism for handling that?"""
    def __init__(self,msg):
        super(CFLError,self).__init__(msg)

class BC():
    """Enumeration of boundary condition names."""
    # This could instead just be implemented as a static dictionary.
    custom     = 0
    extrap    = 1
    periodic   = 2
    wall = 3

# =================== Dummy routines =============
def default_compute_gauge_values(q,aux):
    r"""By default, record values of q at gauges.
    """
    return q

def before_step(solver,solution):
    r"""
    Dummy routine called before each step

    Replace this routine if you want to do something before each time step.
    """
    pass

class Solver(object):
    r"""
    Pyclaw solver superclass.

    The pyclaw.Solver.solver class is an abstract class that should
    not be instantiated; rather, all Solver classes should inherit from it.

    Solver initialization takes one argument -- a Riemann solver:

        >>> from clawpack import pyclaw, riemann
        >>> solver = pyclaw.ClawSolver2D(riemann.euler_4wave_2D)

    after which solver options may be set.
    It is necessary to set the boundary conditions (for q, and
    for aux if an aux array is used):

        >>> solver.bc_lower[0] = pyclaw.BC.extrap
        >>> solver.bc_upper[0] = pyclaw.BC.wall

    Many other options may be set
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
         - ``numsteps`` = Total number of time steps that have been taken

        solver.status is returned by solver.evolve_to_time.

    .. attribute:: before_step

        Function called before each time step is taken.
        The required signature for this function is:

        def before_step(solver,solution)

    .. attribute:: dt_variable

        Whether to allow the time step to vary, ``default = True``.
        If false, the initial time step size is used for all steps.

    .. attribute:: max_steps

        The maximum number of time steps allowed to reach the end time
        requested, ``default = 10000``.  If exceeded, an exception is
        raised.

    .. attribute:: logger

        Default logger for all solvers.  Records information about the run
        and debugging messages (if requested).

    .. attribute:: bc_lower

        (list of ints) Lower boundary condition types, listed in the
        same order as the Dimensions of the Patch.  See Solver.BC for
        an enumeration.

    .. attribute:: bc_upper

        (list of ints) Upper boundary condition types, listed in the
        same order as the Dimensions of the Patch.  See Solver.BC for
        an enumeration.

    .. attribute:: user_bc_lower

        (func) User defined lower boundary condition.
        Fills the values of qbc with the correct boundary values.
        The appropriate signature is:

        def user_bc_lower(patch,dim,t,qbc,num_ghost):

    .. attribute:: user_bc_upper

        (func) User defined upper boundary condition.
        Fills the values of qbc with the correct boundary values.
        The appropriate signature is:

        def user_bc_upper(patch,dim,t,qbc,num_ghost):


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

    @property
    def all_bcs(self):
        return self.bc_lower, self.bc_upper
    @all_bcs.setter
    def all_bcs(self,all_bcs):
        for i in range(self.num_dim):
            self.bc_lower[i] = all_bcs
            self.bc_upper[i] = all_bcs

    #  ======================================================================
    #   Initialization routines
    #  ======================================================================
    def __init__(self,riemann_solver=None,claw_package=None):
        r"""
        Initialize a Solver object

        See :class:`Solver` for full documentation
        """
        # Setup solve logger
        self.logger = logging.getLogger('pyclaw.solver')

        self.dt_initial = 0.1
        self.dt_max = 1e99
        self.max_steps = 10000
        self.dt_variable = True
        self.num_waves = None  # Must be set later to agree with Riemann solver
        self.qbc = None
        self.auxbc = None
        self.rp = None
        self.fmod = None
        self._is_set_up = False
        self._use_old_bc_sig = False
        self.accept_step = True
        self.before_step = None

        # select package to build solver objects from, by default this will be
        # the package that contains the module implementing the derived class
        # for example, if ClawSolver1D is implemented in 'clawpack.petclaw.solver', then
        # the computed claw_package will be 'clawpack.petclaw'

        import sys
        if claw_package is not None and claw_package in sys.modules:
            self.claw_package = sys.modules[claw_package]
        else:
            def get_clawpack_dot_xxx(modname): return modname.rpartition('.')[0].rpartition('.')[0]
            claw_package_name = get_clawpack_dot_xxx(self.__module__)
            if claw_package_name in sys.modules:
                self.claw_package = sys.modules[claw_package_name]
            else:
                raise NotImplementedError("Unable to determine solver package, please provide one")

        # Initialize time stepper values
        self.dt = self.dt_initial
        self.cfl = self.claw_package.CFL(self.cfl_desired)

        # Status Dictionary
        self.status = {'cflmax': -np.inf,
                       'dtmin': np.inf,
                       'dtmax': -np.inf,
                       'numsteps': 0 }

        # No default BCs; user must set them
        self.bc_lower =    [None]*self.num_dim
        self.bc_upper =    [None]*self.num_dim
        self.aux_bc_lower = [None]*self.num_dim
        self.aux_bc_upper = [None]*self.num_dim

        self.user_bc_lower = None
        self.user_bc_upper = None

        self.user_aux_bc_lower = None
        self.user_aux_bc_upper = None

        self.num_eqn   = None
        self.num_waves = None

        self.compute_gauge_values = default_compute_gauge_values
        r"""(function) - Function that computes quantities to be recorded at gauges"""

        self.qbc          = None
        r""" Array to hold ghost cell values.  This is the one that gets passed
        to the Fortran code.  """

        if riemann_solver is not None:
            self.rp = riemann_solver
            rp_name = riemann_solver.__name__.split('.')[-1]
            from clawpack import riemann
            if "ptwise" in rp_name:
                rp_name = rp_name.replace("_ptwise", "")
            self.num_eqn   = riemann.static.num_eqn.get(rp_name,None)
            self.num_waves = riemann.static.num_waves.get(rp_name,None)

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
        reason = None
        if any([bcmeth == BC.custom for bcmeth in self.bc_lower]):
            if self.user_bc_lower is None:
                valid = False
                reason = 'Lower custom BC function has not been set.'
        if any([bcmeth == BC.custom for bcmeth in self.bc_upper]):
            if self.user_bc_upper is None:
                valid = False
                reason = 'Upper custom BC function has not been set.'
        if self.num_waves is None:
            valid = False
            reason = 'solver.num_waves has not been set.'
        if self.num_eqn is None:
            valid = False
            reason = 'solver.num_eqn has not been set.'
        if (None in self.bc_lower) or (None in self.bc_upper):
            valid = False
            reason = 'One of the boundary conditions has not been set.'

        if reason is not None:
            self.logger.debug(reason)
        return valid, reason

    def setup(self,solution):
        r"""
        Stub for solver setup routines.

        This function is called before a set of time steps are taken in order
        to reach tend.  A subclass should extend or override it if it needs to
        perform some setup based on attributes that would be set after the
        initialization routine.  Typically this is initialization that
        requires knowledge of the solution object.
        """

        self._is_set_up = True

    def __del__(self):
        r"""
        Stub for solver teardown routines.

        This function is called at the end of a simulation.
        A subclass should override it only if it needs to
        perform some cleanup, such as deallocating arrays in a Fortran module.
        """
        self._is_set_up = False

    def __str__(self):
        output = "Solver Status:\n"
        for (k,v) in six.iteritems(self.status):
            output = "\n".join((output,"%s = %s" % (k.rjust(25),v)))
        return output

    # ========================================================================
    #  Boundary Conditions
    # ========================================================================
    def _allocate_bc_arrays(self,state):
        r"""
        Create numpy arrays for q and aux with ghost cells attached.
        These arrays are referred to throughout the code as qbc and auxbc.

        This is typically called by solver.setup().
        """
        import inspect
        for fun in (self.user_bc_lower,self.user_bc_upper,self.user_aux_bc_lower,self.user_aux_bc_upper):
            if fun is not None:
                args = inspect.getfullargspec(fun)[0]
                if len(args) == 5:
                    self.logger.warn("""The custom boundary condition
                                        function signature has been changed.
                                        The previous signature will not be
                                        supported in Clawpack 6.0.  Please see
                                        http://www.clawpack.org/pyclaw/solvers.html#change-to-custom-bc-function-signatures
                                        for more information.""")
                    self._use_old_bc_sig = True

        qbc_dim = [n+2*self.num_ghost for n in state.grid.num_cells]
        qbc_dim.insert(0,state.num_eqn)
        self.qbc = np.zeros(qbc_dim,order='F')

        auxbc_dim = [n+2*self.num_ghost for n in state.grid.num_cells]
        auxbc_dim.insert(0,state.num_aux)
        self.auxbc = np.empty(auxbc_dim,order='F')

        self._apply_bcs(state)

    def _apply_bcs(self, state):
        r"""
        Apply boundary conditions to both q and aux arrays.

        In the case of a user-defined boundary condition, both arrays
        qbc and auxbc are passed to the user function.  Typically the
        function would only modify one or the other of them, though this
        is not enforced.

        If the user function only accepts one array argument, we warn
        that this interface has been deprecated.  In Clawpack 6, we will
        drop backward compatibility.

        For parallel runs, we check whether we're actually on a domain
        boundary.  If we are just at an inter-patch boundary, nothing needs to
        be done here.
        """

        import numpy as np

        self.qbc = state.get_qbc_from_q(self.num_ghost, self.qbc)
        if state.num_aux > 0:
            self.auxbc = state.get_auxbc_from_aux(self.num_ghost, self.auxbc)

        grid = state.grid

        for (idim, dim) in enumerate(grid.dimensions):
            # Check if we are on a true boundary
            if state.grid.on_lower_boundary[idim]:

                bcs = []
                if state.num_aux > 0:
                    bcs.append({'array'  : self.auxbc,
                                'type'   : self.aux_bc_lower,
                                'custom_fun' : self.user_aux_bc_lower,
                                'variable' : 'aux'})
                bcs.append({'array'  : self.qbc,
                            'type'   : self.bc_lower,
                            'custom_fun' : self.user_bc_lower,
                            'variable' : 'q'})
                for (i, bc) in enumerate(bcs):

                    if bc['type'][idim] == BC.custom:
                        if not self._use_old_bc_sig:
                            bc['custom_fun'](state, dim, state.t, self.qbc,
                                             self.auxbc, self.num_ghost)
                        else:
                            bc['custom_fun'](state, dim, state.t, bc['array'], self.num_ghost)

                    elif bc['type'][idim] == BC.periodic \
                            and not state.grid.on_upper_boundary[idim]:
                        pass  # In a parallel run, # PETSc handles periodic BCs.

                    else:
                        self._bc_lower(bc['type'][idim], state, dim, state.t,
                                        np.rollaxis(bc['array'], idim+1, 1), idim,
                                        bc['variable'])

            if state.grid.on_upper_boundary[idim]:

                bcs = []
                if state.num_aux > 0:
                    bcs.append({'array'  : self.auxbc,
                                'type'   : self.aux_bc_upper,
                                'custom_fun' : self.user_aux_bc_upper,
                                'variable' : 'aux'})
                bcs.append({'array'  : self.qbc,
                            'type'   : self.bc_upper,
                            'custom_fun' : self.user_bc_upper,
                            'variable' : 'q'})
                for (i, bc) in enumerate(bcs):

                    if bc['type'][idim] == BC.custom:
                        if not self._use_old_bc_sig:
                            bc['custom_fun'](state, dim, state.t, self.qbc,
                                             self.auxbc, self.num_ghost)
                        else:
                            bc['custom_fun'](state, dim, state.t, bc['array'], self.num_ghost)

                    elif bc['type'][idim] == BC.periodic \
                            and not state.grid.on_lower_boundary[idim]:
                        pass  # In a parallel run, # PETSc handles periodic BCs.

                    else:
                        self._bc_upper(bc['type'][idim], state, dim, state.t,
                                        np.rollaxis(bc['array'], idim+1, 1), idim,
                                        bc['variable'])

    def _bc_lower(self, bc_type, state, dim, t, array, idim, name):
        r"""
        Apply lower boundary conditions to array.

        Sets the lower coordinate's ghost cells of *array* depending on what
        :attr:`bc_lower` is.  If :attr:`bc_lower` = 0 then the user
        boundary condition specified by :attr:`user_bc_lower` is used.  Note
        that in this case the function :attr:`user_bc_lower` belongs only to
        this dimension but :attr:`user_bc_lower` could set all user boundary
        conditions at once with the appropriate calling sequence.

        :Input:
         - *patch* - (:class:`Patch`) Patch that the dimension belongs to.

        :Input/Output:
         - *array* - (ndarray(...,num_eqn)) Array with added ghost cells which
           will be set in this routines.
        """

        if bc_type == BC.extrap:
            for i in range(self.num_ghost):
                array[:,i,...] = array[:,self.num_ghost,...]
        elif bc_type == BC.periodic:
            # This process owns the whole patch
            array[:,:self.num_ghost,...] = array[:,-2*self.num_ghost:-self.num_ghost,...]
        elif bc_type == BC.wall:
            if name == 'q':
                for i in range(self.num_ghost):
                    array[:,i,...] = array[:,2*self.num_ghost-1-i,...]
                    # Negate normal velocity
                    array[self.reflect_index[idim],i,...] = \
                        -array[self.reflect_index[idim],2*self.num_ghost-1-i,...]
            else:
                for i in range(self.num_ghost):
                    array[:,i,...] = array[:,2*self.num_ghost-1-i,...]
        else:
            if bc_type is None:
                raise Exception('Lower boundary condition not specified for either q or aux.')
            else:
                raise NotImplementedError("Boundary condition %s not implemented" % bc_type)

    def _bc_upper(self, bc_type, state, dim, t, array, idim, name):
        r"""
        Apply upper boundary conditions to array

        Sets the upper coordinate's ghost cells of *array* depending on what
        :attr:`bc_upper` is.  If :attr:`bc_upper` = 0 then the user
        boundary condition specified by :attr:`user_bc_upper` is used.  Note
        that in this case the function :attr:`user_bc_upper` belongs only to
        this dimension but :attr:`user_bc_upper` could set all user boundary
        conditions at once with the appropriate calling sequence.

        :Input:
         - *patch* - (:class:`Patch`) Patch that the dimension belongs to

        :Input/Output:
         - *array* - (ndarray(...,num_eqn)) Array with added ghost cells which will
           be set in this routines
        """

        if bc_type == BC.extrap:
            for i in range(self.num_ghost):
                array[:,-i-1,...] = array[:,-self.num_ghost-1,...]
        elif bc_type == BC.periodic:
            # This process owns the whole patch
            array[:,-self.num_ghost:,...] = array[:,self.num_ghost:2*self.num_ghost,...]
        elif bc_type == BC.wall:
            if name == 'q':
                for i in range(self.num_ghost):
                    array[:,-i-1,...] = array[:,-2*self.num_ghost+i,...]
                    #  Negate normal velocity
                    array[self.reflect_index[idim],-i-1,...] = \
                        -array[self.reflect_index[idim],-2*self.num_ghost+i,...]
            else:
                for i in range(self.num_ghost):
                    array[:,-i-1,...] = array[:,-2*self.num_ghost+i,...]
        else:
            if bc_type is None:
                raise Exception('Upper boundary condition not specified for either q or aux.')
            else:
                raise NotImplementedError("Boundary condition %s not implemented" % bc_type)

    # ========================================================================
    #  Evolution routines
    # ========================================================================
    def accept_reject_step(self,state):
        cfl = self.cfl.get_cached_max()
        if cfl > self.cfl_max:
            return False
        else:
            return True

    def get_dt_new(self):
        cfl = self.cfl.get_cached_max()
        self.dt = min(self.dt_max,self.dt * self.cfl_desired / cfl)

    def get_dt(self,t,tstart,tend,take_one_step):
        cfl = self.cfl.get_cached_max()
        if self.dt_variable and self.dt_old is not None:
            if cfl > 0.0:
                self.get_dt_new()
                self.status['dtmin'] = min(self.dt, self.status['dtmin'])
                self.status['dtmax'] = max(self.dt, self.status['dtmax'])
            else:
                self.dt = self.dt_max
        else:
            self.dt_old = self.dt

        # Adjust dt so that we hit tend exactly if we are near tend
        if not take_one_step:
            if t + self.dt > tend and tstart < tend:
                self.dt = tend - t
            if tend - t - self.dt < 1.e-14*t:
                self.dt = tend - t

    def evolve_to_time(self, solution, tend=None):
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

        if not self._is_set_up:
            self.setup(solution)

        if tend is None:
            take_one_step = True
        else:
            take_one_step = False

        # Parameters for time-stepping
        tstart = solution.t

        num_steps = 0

        # Setup for the run
        if not self.dt_variable:
            if take_one_step:
                self.max_steps = 1
            else:
                self.max_steps = int((tend - tstart + 1e-10) / self.dt)
                if abs(self.max_steps*self.dt - (tend - tstart)) >      \
                   1e-5 * (tend - tstart):
                    raise Exception('dt does not divide (tend-tstart) and dt '
                                    'is fixed!')
        if self.dt_variable == 1 and self.cfl_desired > self.cfl_max:
            raise Exception('Variable time-stepping and desired CFL > maximum '
                            'CFL')
        if not take_one_step:
            if tend <= tstart:
                self.logger.info("Already at or beyond end time: no evolution ",
                                 "required.")
                self.max_steps = 0

        # Main time-stepping loop
        for n in range(self.max_steps):

            state = solution.state

            # Keep a backup in case we need to retake a time step
            if self.dt_variable:
                q_backup = state.q.copy('F')
                told = solution.t

            if self.before_step is not None:
                self.before_step(self, solution.states[0])

            # Note that the solver may alter dt during the step() routine
            self.step(solution, take_one_step, tstart, tend)

            # Check to make sure that the Courant number was not too large
            cfl = self.cfl.get_cached_max()
            self.accept_step = self.accept_reject_step(state)
            if self.accept_step:
                # Accept this step
                self.status['cflmax'] = max(cfl, self.status['cflmax'])
                if self.dt_variable:
                    solution.t += self.dt
                else:
                    # Avoid roundoff error if dt_variable==False:
                    solution.t = tstart+(n+1)*self.dt

                # Verbose messaging
                self.logger.debug("Step %i  CFL = %f   dt = %f   t = %f"
                                  % (n,cfl,self.dt,solution.t))

                self.write_gauge_values(solution)
                # Increment number of time steps completed
                num_steps += 1
                self.status['numsteps'] += 1

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

            # See if we are finished yet
            if take_one_step:
                break
            elif solution.t >= tend:
                break

        # End of main time-stepping loop -------------------------------------

        if not take_one_step:
            if self.dt_variable and solution.t < tend \
                    and num_steps == self.max_steps:
                raise Exception("Maximum number of timesteps have been taken")

        return self.status

    def step(self,solution):
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
        import numpy as np
        if solution.num_aux == 0:
            aux = None
        for i,gauge in enumerate(solution.state.grid.gauges):
            if self.num_dim == 1:
                ix=gauge[0]
                if solution.num_aux > 0:
                    aux = solution.state.aux[:,ix]
                q=solution.state.q[:,ix]
            elif self.num_dim == 2:
                ix, iy = gauge
                if solution.num_aux > 0:
                    aux = solution.state.aux[:,ix,iy]
                q=solution.state.q[:,ix,iy]
            p=self.compute_gauge_values(q,aux)
            if not hasattr(p,'__iter__'):
                p = [p]
            t=solution.t
            if solution.state.keep_gauges:
                gauge_data = solution.state.gauge_data
                if len(gauge_data) == len(solution.state.grid.gauges):
                    gauge_data[i]=np.vstack((gauge_data[i],np.append(t,p)))
                else:
                    gauge_data.append(np.append(t,p))

            try:
                solution.state.grid.gauge_files[i].write(str(t)+' '+' '.join(str(j)
                                                         for j in p)+'\n')
            except IOError:
                raise Exception("Gauge files are not set up correctly. You should call \
                       \nthe method `setup_gauge_files` of the Grid class object \
                       \nbefore any call for `write_gauge_values` from the Solver class.")


if __name__ == "__main__":
    import doctest
    doctest.testmod()
