# 5.5.0 release
- Some improvements to what happens if a user specifies the wrong file format when reading simulation data.

# 5.4.1 release
- Fixed a few bugs involving things that didn't work in Python 3.
- Added gauges.compare_gauges() function to easily plot gauge data from different directories on a single plot.
- Enabled reading of gauges for 1D simulations.
- Bug fix: loggers no longer try to access pyclaw.io (changed to pyclaw.fileio).
- A PetClaw Patch object can be initialized using a single Dimension (no need for a list).
- Support for ForestClaw moved to src/forestclaw.  Added tests for ForestClaw.
- Added a 2D shallow water example with bathymetry in examples/shallow_2d/sill.py.
- Inform user when boundary conditions are not set.

# 5.4.0 release
- PyClaw is now compatible with both Python 2.7 and Python 3.5.
- There is a new 1D Euler shocktube example that uses the new HLLC Riemann solver.

# 5.3.1 release
- Added a new example that shows how to use a custom Riemann solver.
  It solves a two-species advection-reaction problem in 2D and is in 
  examples/advection_reaction_2d/.
- Added parallel HDF5 file reading and writing (requires installation of parallel HDF library).
- Changed some Python limiter code to give better performance.
- Fixed a test that was mistakenly not being run (euler_1d/shocksine).


# 5.3.0 release
- New time integration methods in SharpClaw: variable-step-size linear multistep methods.  See http://arxiv.org/abs/1504.04107.
- HDF5 file I/O now works.
- Capability to use pointwise Riemann solvers.  Examples in pyclaw/examples/acoustics_1d_homogeneous
  and pyclaw/examples/acoustics_2d_homogeneous.
- New examples:
  - shallow water flow over a sill
  - 2D acoustics on mapped grid with inclusions
  - 3D atmospheric Euler with gravity
- Better enforcement of SSP timestep restriction for Runge-Kutta methods (SharpClaw).
- Improvements to geometry:
    - physical coordinates with ghost cells
    - all grid attributes update automatically when grid is modified
    - better output from `print grid`
    - grid mappings can be modified on the fly
    - grid.c_center and grid.p_center methods give coordinates of a single cell center
    - grid.mapc2p functions the same way it does in visclaw
    - Several new doctests.
- Testing: more doctests, I/O tests; doctests are run on Travis; scipy and hdf5 tests run on Travis.
- Various bug fixes (10 issues closed).

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
