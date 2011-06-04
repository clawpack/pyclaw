.. _problem_setup:

=============================
Setting up your own problem
=============================
The best way to set up a new problem is to find an existing problem setup that
is similar.  The main steps are:

    * Write a Riemann solver (if solving a new system of equations)
    * Set up the Makefile
    * Write the initialization script
    * Write routines for source terms, custom boundary conditions, or other customizations
    * Write a setplot.py file for visualization


Writing a Riemann solver
=============================
The Riemann package has solvers for many hyperbolic systems.  If your problem
involves a new system, you will need to write your own Riemann solver.  Please
then contribute your solver to the package by sending a pull request on Github
or e-mailing one of the developers.

For very simple problems in one dimension, it may be worthwhile to write the
Riemann solver in Python, especially if you are more comfortable with Python
than with Fortran.  For two-dimensional problems, or one-dimensional problems
requiring fine grids (or if you are impatient) the solver should be written
in Fortran.  The best approach is generally to find a similar solver in the
Riemann package and modify it to solve your system.


Setting up the Makefile
===============================
Generally you can just copy the Makefile from an example in pyclaw/apps and
replace the value of `RP_SOURCES`.  Make sure the example you choose has the
same dimensionality.  Also be sure to use the f-wave targets if your Riemann
solver is an f-wave solver.

Writing the initialization script
===================================
This script should:

    * Import the appropriate package (pyclaw or petclaw)
    * Instantiate a :class:`~pyclaw.evolve.solver.Solver` 
    * Set the Riemann solver if using a Python Riemann solver
    * Set solver.mwaves
    * Set the boundary conditions
    * Instantiate some :class:`~pyclaw.grid.Dimension` object(s) and a :class:`~pyclaw.grid.Grid`
    * Set any required global values in aux_global
    * Set grid.meqn and grid.mbc
    * Set the initial condition (grid.q)

Usually the script then instantiates a :class:`~pyclaw.controller.Controller`, sets the
initial solution and solver, and calls :meth:`~pyclaw.controller.Controller.run`.
