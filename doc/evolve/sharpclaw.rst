.. _sharpclaw_solvers:

===============================
SharpClaw Solvers
===============================

The SharpClaw solvers are a collection of solvers that contain the
functionality of the Fortran code SharpClaw, developed in David Ketcheson's
thesis.  The 1D SharpClaw solver contains a pure Python implementation as
well as a wrapped Fortran version.  The 2D solver is in progress but not
available yet.  The SharpClaw solvers provide an interface similar to that
of the classic Clawpack solvers, but with a few different options.
The superclass solvers are not meant
to be used separately but are there to provide common routines for all the
Clawpack solvers.  Please refer to each of the inherited classes for more info
about the methods and attributes they provide each class.  The inheritance
structure is:

.. inheritance-diagram:: pyclaw.sharpclaw.sharpclaw.SharpClawSolver1D

:Example:

    This is a simple example of how to instantiate and evolve a solution to a
    later time :math:`\text{t_end}` using the linearized 1d acoustics Riemann solver

.. doctest::
    
        >>> from pyclaw import SharpClawSolver1D
        >>> solver = SharpClawSolver1D()           # Instantiate a default, 1d solver
        
        >>> solver.evolve_to_time(solution,t_end)  # Evolve the solution to t_end # doctest: +SKIP


:mod:`pyclaw.sharpclaw`
===============================

.. automodule:: pyclaw.sharpclaw
   :members:
   :show-inheritance:
