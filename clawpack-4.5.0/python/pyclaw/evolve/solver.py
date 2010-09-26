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

import sys
import time
import copy
import logging

import numpy as np

# Clawpack modules
from pyclaw.data import Data

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
            'max_steps','dt_variable']
            
    _default_attr_values = {'dt_initial':0.1, 'dt_max':1e99, 'cfl_max':1.0, 
        'cfl_desired':0.9, 'max_steps':1000, 'dt_variable':True}
    
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
        return valid
        
    def setup(self):
        r"""
        Stub for solver setup routines.
        
        This function is called before a set of time steps are taken in order 
        to reach tend.  A subclass should override it only if it needs to 
        perform some setup based on attributes that would be set after the 
        initialization routine.
        
        This function is just a stub here.
        """
        pass
        
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
                self.max_steps = (tend - tstart + 1e-10) / self.dt
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
                old_solution = copy.deepcopy(solutions["n"])
            retake_step = False  # Reset flag
            
            # Take one time step defined by the subclass
            self.step(solutions)

            # Check to make sure that the Courant number was too large
            if self.cfl <= self.cfl_max:
                # Accept this step
                self.status['cflmax'] = max(self.cfl, self.status['cflmax'])
                solutions['n'].t += self.dt 
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
                solutions['n'] = old_solution
                self.logger.debug("Rejecting time step, CFL number too large")
                if self.dt_variable:
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