  .. _pyclaw_tutorial:
  
***************************************
Tutorial: Solve the acoustics equations
***************************************

PyClaw is designed to solve general systems of hyperbolic PDEs of the form

.. math::
   \begin{equation}
        \kappa(x) q_t + A(q,x) q_x = 0.
    \end{equation}

As an example, in this tutorial we'll set up a simulation that solves 
the acoustics equations in one dimension:

.. math::
   \begin{eqnarray}
        &p_t + K u_x = 0\\
        &u_t + \frac{1}{\rho} p_x = 0
    \end{eqnarray}

The key to solving a particular system of equations with PyClaw or other similar 
codes is a Riemann solver.  Riemann solvers for many systems are available as part 
of the clawpack/riemann package. 

We'll assume that you've already followed the :ref:`installation` instructions.

.. note::
   The following instructions show how to set up a problem step-by-step in an
   interactive shell.  If you don't want to type all these commands, you can
   instead::
   
    $ cd $PYCLAW/apps/acoustics/1d/homogeneous 
    $ python acoustics.py iplot=1

Now launch an iPython session and import pyclaw

.. doctest::

    >>> import pyclaw

The Solver
===========
PyClaw includes various algorithms for solving hyperbolic PDEs; each is implemented
in a :class:`~pyclaw.solver.Solver` object.  So the first step is to create a solver

.. doctest::

    >>> solver = pyclaw.ClawSolver1D()

In order to avoid the complication of compiling Fortran code, we'll use a
Riemann solver implemented in Python

.. doctest::

    >>> solver.kernel_language = 'Python'

Now we import the appropriate solver from the `riemann` package and set the 
``solver.rp`` attribute, which is a function handle

.. doctest::

    >>> from riemann import rp_acoustics
    >>> solver.rp = rp_acoustics.rp_acoustics_1d
    >>> solver.num_waves = 2

The ``num_waves`` property indicates the number of waves used in the Riemann solver.

Finally, we set the boundary conditions.  We'll use a wall (wall)
condition at the left boundary and a non-wall (zero-order extrapolation)
condition at the right boundary

.. doctest::

    >>> solver.bc_lower[0] = pyclaw.BC.wall
    >>> solver.bc_upper[0] = pyclaw.BC.extrap

Dimension, Domain, and State
============================
Next we need to set up the grid.  A PyClaw grid is built from dimensions;
a 1D grid requires only 1 dimension

.. doctest::

    >>> x = pyclaw.Dimension('x', -1.0, 1.0, 200)
    
This creates a :class:`~pyclaw.geometry.Dimension` object named ``x``  on the interval ``[-1.0, 1.0]`` with ``200``
cells.  Notice that the calling sequence is similar to numpy's ``linspace``
command, except that the first argument is the name of the dimension.

.. doctest::

    >>> domain = pyclaw.Domain(x)

This creates a :class:`~pyclaw.geometry.Domain` object, which holds information about the cell center
and edge coordinates.  Finally, we set up a :class:`~pyclaw.state.State`
object, which will hold the solution itself

.. doctest::

    >>> state = pyclaw.State(domain,2)

The second argument indicates the number of equations in the hyperbolic
system we're solving: in this case, two.

Initial condition
=================
Now we will set the initial value of the solution

.. doctest::

    >>> xc = domain.grid.x.centers
    >>> from numpy import exp
    >>> state.q[0,:] = exp(-100 * (xc-0.75)**2)
    >>> state.q[1,:] = 0.

The pressure (``state.q[0,:]``) is set to a Gaussian centered at :math:`x=0.75`.
The velocity (``state.q[1,:]``) is set to zero everywhere.

Finally, we put the state into a Solution object

.. doctest::

    >>> solution = pyclaw.Solution(state,domain)

Problem-specific parameters
===========================
The acoustics equations above have some coefficients -- namely, the
bulk modulus :math:`K` and density :math:`\rho` -- that must be defined.
Furthermore, checking the code for the Riemann solver we've chosen
reveals that it expects us to provide values for the impedance :math:`Z`
and sound speed :math:`c`.  These values are stored in a Python dictionary
called problem_data that is a member of the :class:`~pyclaw.state.State`

.. doctest::

    >>> from math import sqrt
    >>> rho = 1.0
    >>> bulk = 1.0
    >>> state.problem_data['rho'] = rho
    >>> state.problem_data['bulk'] = bulk
    >>> state.problem_data['zz'] = sqrt(rho*bulk)
    >>> state.problem_data['cc'] = sqrt(bulk/rho)

The controller
===================
The most convenient way to run a PyClaw simulation is by using a
:class:`~pyclaw.controller.Controller` object.  The controller
directs the solver in advancing the solution and handles output.

.. doctest::

    >>> controller = pyclaw.Controller()
    >>> controller.solution = solution
    >>> controller.solver = solver
    >>> controller.tfinal = 1.0

At last everything is set up!  Now run the simulation

.. doctest::

    >>> controller.run()
    {'dtmin': 0.0010000000000000009, 'dtmax': 0.0090000000000000011, 'numsteps': 12, 'cflmax': 0.90000000000000013}	

This should print out a few lines indicating the output times. It also prints the minimum and maximum tipe-step used, the number of steps required for the computation and the maximum CFL number. The simplest way to plot the solution is

.. doctest::

    >>> from pyclaw import plot
    >>> plot.interactive_plot() # doctest: +SKIP
    

That's it!  Your first PyClaw simulation.  Of course, we've only
scratched the surface of what PyClaw can do, and there are many
important options that haven't been discussed here.  To get an
idea, take a look through the pyclaw/apps directory and try running
some other examples.  It's also a good idea to get more deeply
acquainted with the main :ref:`pyclaw_classes`.
