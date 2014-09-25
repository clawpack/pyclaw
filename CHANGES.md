#5.2.1 release

- **New 3D Euler examples**: Sedov blast, shock tube, and shock-bubble interaction, with IPython notebooks showing how to visualize.  The Sedov example has a test using HDF5 files.
- **New interface to boundary condition routines**: the routines for setting q and aux both get both q and aux now.  The old interface is supported for backward compatibility until 5.3.  All examples have been updated to use the new interface.
- New 2D Euler example: flow over a forward-facing step.  It runs, but the solution blows up with the Roe solvers; this is a known issue.
- In SharpClaw, the total fluctuation solver can now be specified in a manner similar to how the Riemann solver is specified.  The 1D Euler examples demonstrate this.
- The Riemann solvers in clawpack.riemann now have defined constants for each conserved quantity.  Many of the PyClaw examples have been updated to use those for indexing into the q array.

# 5.2.0 release

- **Linear multistep methods for timestepping in SharpClaw**: 
  Previously, only Runge-Kutta time stepping was included in SharpClaw.
  Now, you can use a linear multistep method by setting `solver.time_integrator`
  to an appropriate value.  See the docstring of the SharpClawSolver class for
  details.
- **Characteristic-wise WENO and TVD reconstructions in SharpClaw** are now available.
  The user must provide a routine that computes the eigenvectors.  See 
  `clawpack/pyclaw/examples/euler_1d/`.
- **Logging control**: it is now easy to modify the logging levels interactively,
  without modifying the logger files.  All logging configuration is
  set by default with `pyclaw/log.config`; the file `petclaw/log.config` is
  no longer used.
- It is no longer necessary to compile dummy transverse Riemann solvers when using
  an algorithm with no transverse wave propagation.
- Add functions to compute cell centers of ghost cells.
- Tests run in parallel on Travis.
- All tests now run in under two minutes on most systems.
- PyWENO-generated code updated to match latest PyWENO release.
- Writing and reading of HDF5 and netcdf files now works in serial and parallel.
- Improvements to examples.
- Miscellaneous bug fixes.


# 5.1.0 release

- SharpClaw: now includes a 3D solver.  See the 3D variable acoustics example.
  Note that the 3D solver, like the 2D solver, simply does dimension-by-dimension
  reconstruction and is thus formally 2nd order accurate, regardless of the
  order of WENO reconstruction actually used.  Nevertheless, practical accuracy for
  most multidimensional problems is substantially better than with typical 2nd-order methods.
- I/O: problem_data and mapc2p are now saved to the pickle file and read back in
- pyclaw.examples is now importable (previously it only worked if you did a dev install)
- Controller class: the status dictionary returned by controller.run() now contains
  information relevant to the whole simulation, not merely the last call of evolve_to_time().

# Clawpack 5.0 release: major changes since PyClaw 1.0

- The Riemann solver is selected (by setting solver.rp) at run-time instead of compile-time
- Syntax: many variable names have changed and simpler syntax has been
  introduced for initialization of the main object classes
- Tighter integration with the rest of Clawpack:
    - PyClaw is now part of the larger Clawpack package (and lives in the clawpack namespace)
    - PyClaw docs are part of Clawpack docs (see http://clawpack.github.io/doc/pyclaw/index.html)
- Example applications:
    - The pyclaw/apps/ directory has been renamed to pyclaw/examples
    - examples are importable, so you can get to them even if installed to site-packages
    - New examples have been added: LWR traffic flow and the Shu-Osher problem (1D Euler equations)
- Internal refactoring:
    - reorganization of pyclaw/src directory, so that Fortran and Python files are together
    - All Fortran code has been converted to free-format and uses some F90 features
    - No multiple inheritance in solver classes
- Geometry:
    - The new Domain object allows one-line generation of geometry for
      rectangular domains with uniform grids, without the need to first create
      Dimension objects
- Other features:
    - Arbitrary Runge-Kutta methods can be used for time integration in SharpClaw
    - Option to turn off output to disk (set controller.output_format=None)
- Installation:
    - PyClaw and its dependencies are installed with one line using pip: 'pip install clawpack'
    - All Fortran code is compiled at installation and placed on the PYTHONPATH
    - Makefiles and environment variables are no longer used
- Testing:
    - New test framework; tests can be run in parallel using mpiexec nosetests
    - Testing with Travis
