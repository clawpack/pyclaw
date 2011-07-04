  .. _pyclaw_tutorial:
  
***************
Pyclaw Tutorial
***************

PyClaw is designed to solve general systems of hyperbolic PDEs of the form

.. math::

   \kappa(x) q_t + A(q,x) q_x = 0.

As an example, in this tutorial we'll set up a simulation that solves 
the acoustics equations in one dimension:

.. math::

   p_t + K u_x = 0

   u_t + \frac{1}{\rho} p_x = 0

.. The key to solving a particular system of equations with PyClaw or other similar codes is a Riemann solver.  Riemann solvers for many systems are available as part of the clawpack/riemann package. 

We'll assume that you've already followed the :ref:`installation` instructions.

Now launch an iPython session and import pyclaw::

    >>> import pyclaw

The Solver
===========
PyClaw includes various algorithms for solving hyperbolic PDEs; each is implemented
in a :class:`~pyclaw.solver.Solver` object.  So the first step is to create a solver::

    >>> solver = pyclaw.ClawSolver1D()

In order to avoid the complication of compiling Fortran code, we'll use a
Riemann solver implemented in Python::

    >>> solver.kernel_language = 'Python'

Now we import the appropriate solver from the riemann package and set the 
solver.rp attribute::

    >>> from riemann import rp_acoustics
    >>> solver.rp = rp_acoustics.rp_acoustics_1d

Finally, we set the boundary conditions.  We'll use a reflecting (wall)
condition at the left boundary and a non-reflecting (zero-order extrapolation)
condition at the right boundary::

    >>> solver.mthbc_lower[0] = pyclaw.BC.reflecting
    >>> solver.mthbc_upper[0] = pyclaw.BC.outflow

Dimension, Grid, and State
===========================
Next we need to set up the grid.  A PyClaw grid is built from dimensions;
in our case, we only need 1 dimension::

    >>> x = Dimension('x', -1.0, 1.0, 200)
    
This creates a dimension ``x``  on the interval ``[-1.0, 1.0]`` with ``200``
grid points.  Notice that the calling sequence is similar numpy's linspace
command, except that the first argument is the name of the dimension.

::

    >>> grid = pyclaw.Grid(x)

This creates a grid object, which holds information about the cell center
and edge coordinates.  Finally, we set up a :class:`~pyclaw.state.State`
object, which will hold the solution itself::

    >>> state = pyclaw.State(grid)
    >>> state.meqn = 2

The attribute meqn indicates the number of equations in the hyperbolic
system we're solving.

Initial condition
======================
Now we will set the initial value of the solution::

    >>> xc = grid.x.center
    >>> from numpy import exp
    >>> state.q[0,:] = exp(-100 * (xc-0.75)**2)
    >>> state.q[1,:] = 0.

The pressure (grid.q[0,:]) is set to a Gaussian centered at $x=0.75$.
The velocity (grid.q[1,:]) is set to zero everywhere.

Finally, we put the state into a Solution object::

    >>> solution = pyclaw.Solution(state)


Problem-specific parameters
============================
The acoustics equations above have some coefficients -- namely, the
bulk modulus $K$ and density $\rho$ -- that must be defined.
Furthermore, checking the code for the Riemann solver we've chosen
reveals that it expects us to provide values for the impedance $Z$
and sound speed $c$.  These values are stored in a Python dictionary
called aux_global that is a member of the :class:`~pyclaw.state.State`::

    >>> from math import sqrt
    >>> rho = 1.0
    >>> bulk = 1.0
    >>> state.aux_global['rho'] = rho
    >>> state.aux_global['bulk'] = bulk
    >>> state.aux_global['zz'] = sqrt(rho*bulk)
    >>> state.aux_global['cc'] = sqrt(rho/bulk)

The controller
===================
The most convenient way to run a PyClaw simulation is by using a
:class:`~pyclaw.controller.Controller` object.  The controller
directs the solver in advancing the solution and handles output.

::

    >>> controller = pyclaw.Controller()
    >>> controller.solution = solution
    >>> controller.tfinal = 1.0

At last everything is set up!  Now run the simulation::

    >>> controller.run()

This should print out a few lines indicating the output times.
The simplest way to plot the solution is::

    >>> from pyclaw import plot
    >>> plot.plotInteractive()

That's it!  Your first PyClaw simulation.  Of course, we've only
scratched the surface of what PyClaw can do, and there are many
important options that haven't been discussed here.  To get an
idea, take a look through the pyclaw/apps directory and try running
some other examples.  It's also a good idea to get more deeply
acquainted with the main :ref:`pyclaw_classes`.
