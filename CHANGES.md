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
