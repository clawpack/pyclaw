.. _pyclaw_solver:

==========================
Pyclaw Solver Base Classes
==========================

The solver module provides a template for a solver that Pyclaw knows how to 
interact with.  It is expected that a subclass will override the class method
:meth:`step` in order to provide a complete solver.

The main method of interest defined in this module is the 
:meth:`evolve_to_time` method which evolves the solutions given to the end 
time provided.  Details of how this works can be found below.

:mod:`pyclaw.evolve.solver`
===========================

.. automodule:: pyclaw.evolve.solver
   :members: