#!/usr/bin/env python
# encoding: utf-8
r"""
Module specifying the interface to every solver in pyclaw.

:Authors:
    Kyle T. Mandli (2008-08-19) Initial version
"""
# ============================================================================
#      Copyright (C) 2008 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import time
import copy
import logging

# Clawpack modules
from pyclaw.data import Data

class CFLError(Exception):

    def __init__(self,msg):
        super(CFLError,self).__init__(msg)

class BC():
    """Enumeration of boundary condition names."""
    custom     = 0
    outflow    = 1
    periodic   = 2
    reflecting = 3

# Boundary condition user defined function default
def default_user_bc_lower(grid,dim,t,qbc,mbc):
    r"""
    Fills the values of qbc with the correct boundary values
    
    This is a stub function which will return an exception if called.  If you
    want to use a user defined boundary condition replace this function with
    one of your own.
    """
    raise NotImplementedError("Lower user defined boundary condition unimplemented")

def default_user_bc_upper(grid,dim,t,qbc,mbc):
    r"""
    Fills the values of qbc with the correct boundary values
    
    This is a stub function which will return an exception if called.  If you
    want to use a user defined boundary condition replace this function with
    one of your own.
    """
    raise NotImplementedError("Lower user defined boundary condition unimplemented")


class Solver(object):
    r"""
    Pyclaw solver superclass

    All solvers should inherit from this object as it defines the interface 
    that the Pyclaw expects for solvers.  This mainly means the evolve_to_time
    exists and the solver can be initialized and called correctly.

    .. attribute:: dt
        
        Current time step, ``default = 0.1``
        
    .. attribute:: cfl
        
        Current Courant/Freidrichs/Levy number, ``default = 1.0``
    
    .. attribute:: status
        
        Dictionary of status values for the solver with the following keys:
         - ``cflmax`` = Maximum CFL condition reached
         - ``dtmin`` = Minimum time step taken
         - ``dtmax`` = Maximum time step taken
         - ``numsteps`` = Current number of time steps that have been taken
    
    .. attribute:: dt_variable
    
        Whether to allow the time step to vary, ``default = True``
        
    .. attribute:: max_steps
    
        The maximum number of time steps allowd to reach the end time 
        requested, ``default = 1000``
    
    .. attribute:: times
        
        A list of run times taken by the solver, each request to evolve a 
        solution to a new time a new timing is appended to this list.
    
    .. attribute:: logger
    
        Default logger for all solvers

    - *mthbc_lower* - (int) Lower boundary condition method to be used
    - *mthbc_upper* - (int) Upper boundary condition method to be used
    - *user_bc_lower* - (func) User defined lower boundary condition
    - *user_bc_upper* - (func) User defined upper boundary condition
 
        
    :Initialization:
    
    Input:
     - *data* - (:class:`Data`) Data object to initialize the solver with
    
    Output:
     - (:class:`Solver`) - Initialized Solver object
    
    :Version: 1.0 (2008-08-19)
    """
    
    # Note that some of these names are for compatibility with the claw.data
    # files and the data objects, they are not actually required and do not
    # exist past initialization
    _required_attrs = ['dt_initial','dt_max','cfl_max','cfl_desired',
            'max_steps','dt_variable','mbc']
            
    _default_attr_values = {'dt_initial':0.1, 'dt_max':1e99, 'cfl_max':1.0, 
        'cfl_desired':0.9, 'max_steps':1000, 'dt_variable':True,'mbc':2}
    
    #  ======================================================================
    #   Initialization routines
    #  ======================================================================
    def __init__(self,data=None):
        r"""
        Initialize a Solver object
        
        See :class:`Solver` for full documentation
        """ 
        # Setup solve logger
        self.logger = logging.getLogger('evolve')

        # Set default values
        for (k,v) in self._default_attr_values.iteritems():
            self.__dict__.setdefault(k,v)

        # Set data values based on the data object
        if data is not None and isinstance(data,Data):
            for attr in self._required_attrs:
                if hasattr(data,attr):
                    setattr(self,attr,getattr(data,attr))
        
        # Initialize time stepper values
        self.dt = self._default_attr_values['dt_initial']
        self.cfl = self._default_attr_values['cfl_desired']
        
        # Status Dictionary
        self.status = {'cflmax':self.cfl,
                    'dtmin':self.dt, 
                    'dtmax':self.dt,
                    'numsteps':0 }
        
        # Profile times
        self.times = []

        self.mthbc_lower = [2]*self.ndim
        self.mthbc_upper = [2]*self.ndim
        self.mthauxbc_lower = [2]*self.ndim
        self.mthauxbc_upper = [2]*self.ndim
        
        self.user_bc_lower = default_user_bc_lower
        r"""(func) - User defined boundary condition function, lower.  
        ``default = None``
        """
        self.user_bc_upper = default_user_bc_upper
        r"""(func) - User defined boundary condition function, upper. 
        ``default = None``"""

        self.user_aux_bc_lower = default_user_bc_lower
        self.user_aux_bc_upper = default_user_bc_upper

    def __str__(self):
        output = "Solver Status:\n"
        for (k,v) in self.status.iteritems():
            output = "\n".join((output,"%s = %s" % (k.rjust(25),v)))
        return output

    # ========================================================================
    #  Solver setup and validation routines
    # ========================================================================
    def is_valid(self):
        r"""
        Checks that all required attributes are set
        
        Checks to make sure that all the required attributes for the solver 
        have been set correctly.  All required attributes that need to be set 
        are contained in the attributes list of the class.
        
        Will post debug level logging message of which required attributes 
        have not been set.
        
        :Output:
         - *valid* - (bool) True if the solver is valid, False otherwise
        
        """
        valid = True
        for key in self._required_attrs:
            if not self.__dict__.has_key(key):
                self.logger.info('%s is not present.' % key)
                valid = False
        #for i,bcmeth in enumerate(self.mthbc_lower):
        #    if bcmeth == BC.custom:
        #        try:
        #            self.user_bc_lower(self,dim,None)
        #        except NotImplementedError:
        #            logger.debug('Lower BC function for %s has not been set.' % dim.name)
        #            valid = False
        #        except:
        #            pass
        #    if dim.mthbc_upper == 0:
        #        try:
        #            dim.user_bc_upper(self,dim,None)
        #        except NotImplementedError:
        #            logger.debug('Upper BC function for %s has not been set.' % dim.name)
        #            valid = False
        #        except:
        #            pass

        return valid
        
    def setup(self,solutions):
        r"""
        Stub for solver setup routines.
        
        This function is called before a set of time steps are taken in order 
        to reach tend.  A subclass should override it only if it needs to 
        perform some setup based on attributes that would be set after the 
        initialization routine.
        
        This function is just a stub here.
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
        

    def allocate_rk_stages(self,solutions):
        from pyclaw.state import State

        if self.time_integrator   == 'Euler':  nregisters=1
        elif self.time_integrator == 'SSP33':  nregisters=2
        elif self.time_integrator == 'SSP104': nregisters=3
 
        state = solutions['n'].states[0]
        self.rk_stages = []
        for i in range(nregisters-1):
            self.rk_stages.append(State(state.grid))
            self.rk_stages[-1].meqn = state.meqn
            self.rk_stages[-1].maux = state.maux
            self.rk_stages[-1].aux_global       = state.aux_global
            self.rk_stages[-1].t                = state.t
            if state.maux > 0:
                self.rk_stages[-1].aux              = state.aux


    # ========================================================================
    #  Boundary Conditions
    # ========================================================================    
    def append_ghost_cells(self,state):
        """
        Returns q with ghost cells attached.  For Solver, this means
        just creating a copy of q with extra cells.
        """
        import numpy as np

        grid = state.grid
        mbc = self.mbc
        dims = [n + 2*mbc for n in grid.ng]
        dims.insert(0,state.meqn)
        qbc = np.zeros(dims,order = 'F')
        if grid.ndim == 1:
            qbc[:,mbc:-mbc] = state.q
        elif grid.ndim == 2:
            qbc[:,mbc:-mbc,mbc:-mbc] = state.q
        elif grid.ndim == 3:
            qbc[:,mbc:-mbc,mbc:-mbc,mbc:-mbc] = state.q
        return qbc
 
    def append_ghost_cells_to_aux(self,state):
        """
        Returns aux with ghost cells attached.  For the serial Solver, this means
        just creating a copy of aux with extra cells.
        """
        import numpy as np

        grid = state.grid
        mbc = self.mbc
        dims = [n + 2*self.mbc for n in grid.ng]
        dims.insert(0,state.maux)
        auxbc = np.zeros(dims,order = 'F')
        if grid.ndim == 1:
            auxbc[:,mbc:-mbc] = state.aux
        elif grid.ndim == 2:
            auxbc[:,mbc:-mbc,mbc:-mbc] = state.aux
        elif grid.ndim == 3:
            auxbc[:,mbc:-mbc,mbc:-mbc,mbc:-mbc] = state.aux
        return auxbc


    def qbc(self,state):
        r"""
        Appends boundary cells to q and fills them with appropriate values.
    
        This function returns an array of dimension determined by the 
        :attr:`mbc` attribute.  The type of boundary condition set is 
        determined by :attr:`mthbc_lower` and :attr:`mthbc_upper` for the 
        approprate dimension.  Valid values for :attr:`mthbc_lower` and 
        :attr:`mthbc_upper` include:
        
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

        qbc=self.append_ghost_cells(state)
        grid = state.grid
       
        for idim,dim in enumerate(grid.dimensions):
            # First check if we are actually on the boundary
            # (in case of a parallel run)
            if dim.nstart == 0:
                # If a user defined boundary condition is being used, send it on,
                # otherwise roll the axis to front position and operate on it
                if self.mthbc_lower[idim] == BC.custom:
                    self.qbc_lower(grid,dim,state.t,qbc,idim)
                elif self.mthbc_lower[idim] == BC.periodic:
                    if dim.nend == dim.n:
                        # This process owns the whole grid
                        self.qbc_lower(grid,dim,state.t,np.rollaxis(qbc,idim+1,1),idim)
                    else:
                        pass #Handled automatically by PETSc
                else:
                    self.qbc_lower(grid,dim,state.t,np.rollaxis(qbc,idim+1,1),idim)

            if dim.nend == dim.n :
                if self.mthbc_upper[idim] == BC.custom:
                    self.qbc_upper(grid,dim,state.t,qbc,idim)
                elif self.mthbc_upper[idim] == BC.periodic:
                    if dim.nstart == 0:
                        # This process owns the whole grid
                        self.qbc_upper(grid,dim,state.t,np.rollaxis(qbc,idim+1,1),idim)
                    else:
                        pass #Handled automatically by PETSc
                else:
                    self.qbc_upper(grid,dim,state.t,np.rollaxis(qbc,idim+1,1),idim)
            
        return qbc


    def qbc_lower(self,grid,dim,t,qbc,idim):
        r"""
        Apply lower boundary conditions to qbc
        
        Sets the lower coordinate's ghost cells of *qbc* depending on what 
        :attr:`mthbc_lower` is.  If :attr:`mthbc_lower` = 0 then the user 
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
        import numpy as np

        if self.mthbc_lower[idim] == BC.custom: 
            self.user_bc_lower(grid,dim,t,qbc,self.mbc)
        elif self.mthbc_lower[idim] == BC.outflow:
            for i in xrange(self.mbc):
                qbc[:,i,...] = qbc[:,self.mbc,...]
        elif self.mthbc_lower[idim] == BC.periodic:
            # This process owns the whole grid
            qbc[:,:self.mbc,...] = qbc[:,-2*self.mbc:-self.mbc,...]
        elif self.mthbc_lower[idim] == BC.reflecting:
            for i in xrange(self.mbc):
                qbc[:,i,...] = qbc[:,2*self.mbc-1-i,...]
                qbc[idim+1,i,...] = -qbc[idim+1,2*self.mbc-1-i,...] # Negate normal velocity
        else:
            raise NotImplementedError("Boundary condition %s not implemented" % x.mthbc_lower)


    def qbc_upper(self,grid,dim,t,qbc,idim):
        r"""
        Apply upper boundary conditions to qbc
        
        Sets the upper coordinate's ghost cells of *qbc* depending on what 
        :attr:`mthbc_upper` is.  If :attr:`mthbc_upper` = 0 then the user 
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
 
        if self.mthbc_upper[idim] == BC.custom:
            self.user_bc_upper(grid,dim,t,qbc,self.mbc)
        elif self.mthbc_upper[idim] == BC.outflow:
            for i in xrange(self.mbc):
                qbc[:,-i-1,...] = qbc[:,-self.mbc-1,...] 
        elif self.mthbc_upper[idim] == BC.periodic:
            # This process owns the whole grid
            qbc[:,-self.mbc:,...] = qbc[:,self.mbc:2*self.mbc,...]
        elif self.mthbc_upper[idim] == BC.reflecting:
            for i in xrange(self.mbc):
                qbc[:,-i-1,...] = qbc[:,-2*self.mbc+i,...]
                qbc[idim+1,-i-1,...] = -qbc[idim+1,-2*self.mbc+i,...] # Negate normal velocity
        else:
            raise NotImplementedError("Boundary condition %s not implemented" % x.mthbc_lower)



    def auxbc(self,state):
        r"""
        Appends boundary cells to aux and fills them with appropriate values.
    
        This function returns an array of dimension determined by the 
        :attr:`mbc` attribute.  The type of boundary condition set is 
        determined by :attr:`mthauxbc_lower` and :attr:`mthauxbc_upper` for the 
        approprate dimension.  Valid values for :attr:`mthauxbc_lower` and 
        :attr:`mthauxbc_upper` include:
        
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

        auxbc=self.append_ghost_cells_to_aux(state)
        grid = state.grid
       
        for idim,dim in enumerate(grid.dimensions):
            # First check if we are actually on the boundary
            # (in case of a parallel run)
            if dim.nstart == 0:
                # If a user defined boundary condition is being used, send it on,
                # otherwise roll the axis to front position and operate on it
                if self.mthauxbc_lower[idim] == BC.custom:
                    self.auxbc_lower(grid,dim,state.t,auxbc,idim)
                elif self.mthauxbc_lower[idim] == BC.periodic:
                    if dim.nend == dim.n:
                        # This process owns the whole grid
                        self.auxbc_lower(grid,dim,state.t,np.rollaxis(auxbc,idim+1,1),idim)
                    else:
                        pass #Handled automatically by PETSc
                else:
                    self.auxbc_lower(grid,dim,state.t,np.rollaxis(auxbc,idim+1,1),idim)

            if dim.nend == dim.n :
                if self.mthauxbc_upper[idim] == BC.custom:
                    self.auxbc_upper(grid,dim,state.t,auxbc,idim)
                elif self.mthauxbc_upper[idim] == BC.periodic:
                    if dim.nstart == 0:
                        # This process owns the whole grid
                        self.auxbc_upper(grid,dim,state.t,np.rollaxis(auxbc,idim+1,1),idim)
                    else:
                        pass #Handled automatically by PETSc
                else:
                    self.auxbc_upper(grid,dim,state.t,np.rollaxis(auxbc,idim+1,1),idim)
            
        return auxbc


    def auxbc_lower(self,grid,dim,t,auxbc,idim):
        r"""
        Apply lower boundary conditions to auxbc
        
        Sets the lower coordinate's ghost cells of *auxbc* depending on what 
        :attr:`mthauxbc_lower` is.  If :attr:`mthauxbc_lower` = 0 then the user 
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
        import numpy as np

        if self.mthauxbc_lower[idim] == BC.custom: 
            self.user_aux_bc_lower(grid,dim,t,auxbc,self.mbc)
        elif self.mthauxbc_lower[idim] == BC.outflow:
            for i in xrange(self.mbc):
                auxbc[:,i,...] = auxbc[:,self.mbc,...]
        elif self.mthauxbc_lower[idim] == BC.periodic:
            # This process owns the whole grid
            auxbc[:,:self.mbc,...] = auxbc[:,-2*self.mbc:-self.mbc,...]
        elif self.mthauxbc_lower[idim] == BC.reflecting:
            for i in xrange(self.mbc):
                auxbc[:,i,...] = auxbc[:,2*self.mbc-1-i,...]
        else:
            raise NotImplementedError("Boundary condition %s not implemented" % x.mthauxbc_lower)


    def auxbc_upper(self,grid,dim,t,auxbc,idim):
        r"""
        Apply upper boundary conditions to auxbc
        
        Sets the upper coordinate's ghost cells of *auxbc* depending on what 
        :attr:`mthauxbc_upper` is.  If :attr:`mthauxbc_upper` = 0 then the user 
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
 
        if self.mthauxbc_upper[idim] == BC.custom:
            self.user_aux_bc_upper(grid,dim,t,auxbc,self.mbc)
        elif self.mthauxbc_upper[idim] == BC.outflow:
            for i in xrange(self.mbc):
                auxbc[:,-i-1,...] = auxbc[:,-self.mbc-1,...] 
        elif self.mthauxbc_upper[idim] == BC.periodic:
            # This process owns the whole grid
            auxbc[:,-self.mbc:,...] = auxbc[:,self.mbc:2*self.mbc,...]
        elif self.mthauxbc_upper[idim] == BC.reflecting:
            for i in xrange(self.mbc):
                auxbc[:,-i-1,...] = auxbc[:,-2*self.mbc+i,...]
        else:
            raise NotImplementedError("Boundary condition %s not implemented" % x.mthauxbc_lower)


    # ========================================================================
    #  Evolution routines
    # ========================================================================
    def evolve_to_time(self,*args):
        r"""
        Evolve solutions['n'] to tend
        
        This method contains the machinery to evolve the solution object in
        ``solutions['n']`` to the requested end time tend if given, or one 
        step if not.  The solutions dictionary is provided as a generic 
        interface to multistep methods that may require more than one solution
        at multiple times to evolve where it is understood that the solution 
        at ``solutions['n']`` is the solution at the current time step that is
        to be evolved.
        
        :Input:
         - *solutions* - (:class:`Solution`) Dictionary of Solutions to be 
           evolved
         - *tend* - (float) The end time to evolve to, if not provided then 
           the method will take a single time step.
            
        :Output:
         - (dict) - Returns the status dictionary of the solver
        """
        
        # Parse arguments
        solutions = args[0]
        if len(args) > 1:
            tend = args[1]
            take_one_step = False
        else:
            tend = None
            take_one_step = True
            
        # Profile time
        timing = time.time()

        # Parameters for time stepping
        retake_step = False
        tstart = solutions['n'].t

        # Reset status dictionary
        self.status['cflmax'] = self.cfl
        self.status['dtmin'] = self.dt
        self.status['dtmax'] = self.dt
        self.status['numsteps'] = 0

        # Setup for the run
        if not self.dt_variable:
            if tend < tstart or take_one_step:
                self.max_steps = 1
            else:
                self.max_steps = int((tend - tstart + 1e-10) / self.dt)
                if abs(self.max_steps*self.dt - (tend - tstart)) > 1e-5 * (tend-tstart):
                    raise Exception('dt does not divide (tend-tstart) and dt is fixed!')
        if self.dt_variable == 1 and self.cfl_desired > self.cfl_max:
            raise Exception('Variable time stepping and desired CFL > maximum CFL')
        if tend == tstart:
            self.logger.info("Already at end time, no evolution required.")
            self.max_steps = 0
                
        # Main time stepping loop
        for n in xrange(self.max_steps):
        
            # Adjust dt so that we hit tend exactly if we are near tend
            if solutions['n'].t + self.dt > tend and tstart < tend and not take_one_step:
                self.dt = tend - solutions['n'].t 

            # In case we need to retake a time step
            if self.dt_variable:
                # pass
                #Temporarily HACKed to avoid slowdown!
                #old_solution = copy.deepcopy(solutions["n"])
                qold = copy.copy(solutions["n"].state.q)
                told = solutions["n"].t
            retake_step = False  # Reset flag
            
            # Take one time step defined by the subclass
            self.step(solutions)

            # Check to make sure that the Courant number was not too large
            if self.cfl <= self.cfl_max:
                # Accept this step
                self.status['cflmax'] = max(self.cfl, self.status['cflmax'])
                if self.dt_variable==True:
                    solutions['n'].t += self.dt 
                else:
                    #Avoid roundoff error if dt_variable=False:
                    solutions['n'].t = tstart+(n+1)*self.dt
                # Verbose messaging
                self.logger.debug("Step %i  CFL = %f   dt = %f   t = %f"
                    % (n,self.cfl,self.dt,solutions['n'].t))
                    
                # Increment number of time steps completed
                self.status['numsteps'] += 1
                # See if we are finished yet
                if solutions['n'].t >= tend or take_one_step:
                    break
            else:
                # Reject this step
                self.logger.debug("Rejecting time step, CFL number too large")
                if self.dt_variable:
                    solutions['n'].state.q = qold
                    solutions['n'].t = told
                    # Retake step
                    retake_step = True
                else:
                    # Give up, we cannot adapt, abort
                    self.status['cflmax'] = \
                        max(self.cfl, self.status['cflmax'])
                    raise Exception('CFL to large, giving up!')
                    
            # Choose new time step
            if self.dt_variable:
                if self.cfl > 0.0:
                    self.dt = min(self.dt_max,self.dt * self.cfl_desired 
                                    / self.cfl)
                    self.status['dtmin'] = min(self.dt, self.status['dtmin'])
                    self.status['dtmax'] = max(self.dt, self.status['dtmax'])
                else:
                    self.dt = self.dt_max
        
        # End of main time stepping loop -------------------------------------

        if self.dt_variable and solutions['n'].t < tend \
                and self.status['numsteps'] == self.max_steps:
            raise Exception("Maximum number of timesteps have been taken")

        # Profile time
        self.times.append(time.time() - timing)
        
        return self.status

    def step(self):
        r"""
        Take one step
        
        This method is only a stub and should be overridden by all solvers who
        would like to use the default time stepping in evolve_to_time.
        """
        raise NotImplementedError("No stepping routine has been defined!")

    def communicateCFL(self):
        r"""
        Dummy function, only here for PetClaw to override.
        """
        pass
            


