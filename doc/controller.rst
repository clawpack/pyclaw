.. _pyclaw_controller:

***********************
Pyclaw Controller Class
***********************

The pyclaw controller object is a convenience class for running simulations
based on the classic clawpack formats and output specifications.  It allows 
for a variety of output time specifications, output styles and other ways to 
keep a simulation organized.

The main way to use a Controller object then is to provide it with an
appropriate :class:`~pyclaw.evolve.solver.Solver` and initial 
:class:`~pyclaw.solution.Solution` object.  Then specify what kind of output
you would like different than the defaults (see 
:class:`~pyclaw.controller.Controller` for
details on what those are).  Then simply call 
:meth:`~pyclaw.controller.Controller.run` in order to run the desired 
simulation.  

::

    >>> import pyclaw.controller as controller
    >>> claw = controller.Controller()            # Instantiate a new controller
    >>> claw.solver = my_solver                   # Assign a solver
    >>> claw.solutions['n'] = my_initial_solution # Assign an initial condition

Here we would set a variety of run parameters such as ``tfinal``, 
``keep_copy`` if we wanted to plot the solutions immediately, or 
``output_format`` to specify a format other than ``ascii`` or no output files 
if we are going to use ``keep_copy = True``.  After we are all set up we just
need to call the controller's :meth:`run` method and off we go.
    
::    

    >>> claw.run()

Please see the :ref:`pyclaw_tutorial` for a detailed example of how this would 
work in its entirety.
   
:class:`pyclaw.controller.Controller`
=====================================

.. autoclass:: pyclaw.controller.Controller
   :members:
   :member-order: groupwise
   