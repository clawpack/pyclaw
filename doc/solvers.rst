.. contents::

.. _solvers:

**********************************************
Using PyClaw's solvers: Classic and SharpClaw
**********************************************

At present, PyClaw includes two types of solvers:

    * Classic: the original Clawpack algorithms, in 1/2/3D
    * SharpClaw: higher-order wave propagation using WENO reconstruction and
      Runge-Kutta integration, in 1/2D


Typically, all that is needed to select a different solver is to specify
it in the problem script, e.g.

    >>> solver = pyclaw.ClawSolver2D()

for the Classic solver or

    >>> solver = pyclaw.SharpClawSolver2D()

for the SharpClaw solver.  Most of the applications distributed with PyClaw
are set up to use either solver, depending on the value of the command line option
`solver_type`, which should be set to `classic` or `sharpclaw`.

Typically, for a given grid resolution, the SharpClaw solvers are more accurate
but also more computationally expensive.
For typical problems involving shocks, the Classic solvers are recommended.
For problems involving high-frequency waves, turbulence, or smooth solutions,
the SharpClaw solvers may give more accurate solutions at less cost.  This
is an active area of research and you may wish to experiment with both solvers.

Future plans include incorporation of finite-difference and discontinuous Galerkin
solvers.

Key differences between the Classic and SharpClaw solvers are:

    * The source term routine for the Classic solver should return the integral of
      the source term over a step, while the source term routine for SharpClaw
      should return the instantaneous value of the source term.

    * The solvers have different options.  For a list of options and possible
      values, see the documentation of the :class:`~pyclaw.classic.solver.ClawSolver` and 
      :class:`~pyclaw.sharpclaw.solver.SharpClawSolver` classes.


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

.. inheritance-diagram:: pyclaw.sharpclaw.solver.SharpClawSolver1D pyclaw.sharpclaw.solver.SharpClawSolver2D

:Example:

    This is a simple example of how to instantiate and evolve a solution to a
    later time :math:`\text{t_end}` using the 1D acoustics Riemann solver.

.. doctest::
    
        >>> from pyclaw import SharpClawSolver1D
        >>> solver = SharpClawSolver1D()           # Instantiate a default, 1d solver
        
        >>> solver.evolve_to_time(solution,t_end)  # Evolve the solution to t_end # doctest: +SKIP


:mod:`pyclaw.sharpclaw`
===============================

.. autoclass:: pyclaw.sharpclaw.solver.SharpClawSolver
   :members:



.. _pyclaw_clawpack_solvers:

===============================
Pyclaw Classic Clawpack Solvers
===============================

The pyclaw classic clawpack solvers are a collection of solvers that represent
the functionality of the older versions of clawpack.  It comes in two forms, a
pure python version and a python wrapping of the fortran libraries.  All of the
solvers available provide the same basic interface and provide the same 
options as the old versions of clawpack.  The superclass solvers are not meant
to be used separately but there to provide common routines for all the
Clawpack solvers.  Please refer to each of the inherited classes for more info
about the methods and attributes they provide each class.  The inheritance
structure is:

.. inheritance-diagram:: clawpack.pyclaw.classic.solver.ClawSolver1D clawpack.pyclaw.classic.solver.ClawSolver2D clawpack.pyclaw.classic.solver.ClawSolver3D

:Example:

    This is a simple example of how to instantiate and evolve a solution to a
    later time :math:`\text{t_end}` using the linearized 1d acoustics Riemann solver
    
.. doctest::

    >>> from clawpack import pyclaw
    >>> solver = pyclaw.ClawSolver1D()                   # Instantiate a default, 1d solver
    >>> solver.set_mthlim = pyclaw.limiters.tvd.vanleer  # Use the van Leer limiter
    >>> solver.dt = 0.0001                               # Set the initial time step
    >>> solver.max_steps = 500                           # Set the maximum number of time steps

.. doctest::

    >>> solver.evolve_to_time(solution,t_end)  # Evolve the solution to t_end  # doctest: +SKIP


:mod:`pyclaw.classic.solver`
=============================

.. autoclass:: clawpack.pyclaw.classic.solver.ClawSolver
   :members:
