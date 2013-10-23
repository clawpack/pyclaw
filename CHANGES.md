# Major changes since PyClaw 1.0

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
