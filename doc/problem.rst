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

In some cases the user might want to reuse a complete set of Fortran routines 
(source term, initialization of the auxiliary variables, custom BC, etc.) that 
were setup in a Clawpack example. This will require a direct use of the Fortran
to Python interface `f2py <http://www.scipy.org/F2py>`_. More details can be 
found at :ref:`port_Example`.


Writing the initialization script
===================================
This script should:

    * Import the appropriate package (pyclaw or petclaw)
    * Instantiate a :class:`~pyclaw.solver.Solver` 
    * Set the Riemann solver if using a Python Riemann solver
    * Set solver.mwaves to the number of waves used in the Riemann solver
    * Set the boundary conditions
    * Instantiate some :class:`~pyclaw.grid.Dimension` object(s) and a :class:`~pyclaw.grid.Grid`
    * Set any required global values in aux_global
    * Set grid.meqn and grid.mbc
    * Set the initial condition (grid.q)

Usually the script then instantiates a :class:`~pyclaw.controller.Controller`, sets the
initial solution and solver, and calls :meth:`~pyclaw.controller.Controller.run`.

Setting initial conditions
----------------------------
Once you have initialize a State object, it contains a member state.q
whose first dimension is meqn and whose remaining dimensions are those
of the grid.  Now you must set the initial condition.  For instance::

    >>> Y,X = np.meshgrid(grid.y.center,grid.x.center)
    >>> r = np.sqrt(X**2 + Y**2)
    >>> width=0.2
    >>> state.q[0,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))
    >>> state.q[1,:,:] = 0.


Setting auxiliary variables
----------------------------
If the problem involves coefficients that vary in space or a mapped grid,
the required fields are stored in state.aux.  In order to use such fields,
you must pass the maux argument to the State initialization::

    >>> state = pyclaw.State(grid,meqn,maux)

The number of fields in state.aux (i.e., the length of its first dimension)
is set equal to maux.  The values of state.aux are set in the same way
as those of state.q.

Setting boundary conditions
----------------------------
The boundary conditions are specified through solver.bc_lower and solver.bc_upper,
each of which is a list of length solver.ndim.  The ordering of the boundary conditions
in each list is the same as the ordering of the Dimensions in the Grid; typically :math:`(x,y)`.
Thus solver.bc_lower[0] specifies the boundary condition at the left boundary and
solver.bc_upper[0] specifies the condition at the right boundary.  Similarly,
solver.bc_lower[1] and solver.bc_upper[1] specify the boundary conditions at the
top and bottom of the domain.

PyClaw includes the following built-in boundary condition implementations:

    * pyclaw.BC.periodic - periodic

    * pyclaw.BC.outflow - zero-order extrapolation

    * pyclaw.BC.reflecting - solid wall conditions, assuming that the 2nd/3rd component
                             of q is the normal velocity in x/y.

Other boundary conditions can be implemented by using pyclaw.BC.custom, and
providing a custom BC function.  The attribute solver.user_bc_lower/upper must
be set to the corresponding function handle.  For instance::

    >>> def custombc(state,dim,t,qbc,mbc):
    >>>     for i in xrange(mbc):
    >>>         qbc[0,i,:] = q0
    >>>
    >>> solver.bc_lower[0]=pyclaw.BC.custom
    >>> solver.user_bc_lower=shockbc

If the state.aux array is used, boundary conditions must be set for it
in a similar way, using solver.aux_bc_lower and solver.aux_bc_upper.
Note that although state is passed to the BC routines, they should
NEVER modify state.  Rather, they should modify qbc/auxbc.

Setting solver options
----------------------------

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

Adding source terms
==============================
Non-hyperbolic terms (representing, e.g., reaction or diffusion) can be included
in a PyClaw simulation by providing an appropriate function handle to 

    * solver.step_src if using Classic Clawpack.  In this case, the function
      specified should modify q by taking a step on the equation :math:`q_t = \psi(q)`.

    * solver.dq_src if using SharpClaw.  In this case, the function should
      return :math:`\Delta t \cdot \psi(q)`.

For an example, see pyclaw/apps/euler/2d/shockbubble/shockbubble.py.

Setting up the Makefile
===============================
Generally you can just copy the Makefile from an example in pyclaw/apps and
replace the value of `RP_SOURCES`.  Make sure the example you choose has the
same dimensionality.  Also be sure to use the f-wave targets if your Riemann
solver is an f-wave solver.


