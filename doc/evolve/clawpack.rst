.. _pyclaw_clawpack_solvers:

===============================
Pyclaw Classic Clawpack Solvers
===============================

The pyclaw classic clawpack solvers are a collection of solvers that represent
the functionality of the older versions of clawpack.  It comes in two forms, a
pure python version and a python wrapping the fortran libraries.  All of the
solvers available provide the same basic interface and provide the same 
options as the old versions of clawpack.  The superclass solvers are not meant
to be used separately but there to provide common routines for all the
Clawpack solvers.  Please refer to each of the inherited classes for more info
about the methods and attributes they provide each class.  The inheritance
structure is:

.. inheritance-diagram:: pyclaw.clawpack.clawpack.ClawSolver1D pyclaw.clawpack.clawpack.ClawSolver2D

:Example:

    This is a simple example of how to instantiate and evolve a solution to a
    later time :math:`\text{t_end}` using the linearized 1d acoustics Riemann solver
    
.. doctest::

    >>> import pyclaw
    >>> solver = pyclaw.ClawSolver1D()                   # Instantiate a default, 1d solver
    >>> solver.set_mthlim = pyclaw.limiters.tvd.vanleer  # Use the van Leer limiter
    >>> solver.dt = 0.0001                               # Set the initial time step
    >>> solver.max_steps = 500                           # Set the maximum number of time steps

.. doctest::

    >>> solver.evolve_to_time(solution,t_end)  # Evolve the solution to t_end  # doctest: +SKIP


:mod:`pyclaw.clawpack`
=============================

.. automodule:: pyclaw.clawpack
   :members:
   :show-inheritance:
