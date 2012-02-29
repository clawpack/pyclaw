.. _output:

***********************
Advanced output options
***********************
PyClaw supports options to output more
than just the solution :math:`q`.  It can provide:

    * Output of derived quantities computed from :math:`q`; for instance,
      pressure (not a conserved quantity) could be computed from density
      and energy.
    * Output of scalar functionals, such as the total mass summed over the whole grid.
    * Output of gauge values, which are time traces of the solution at a
      single point.

Derived quantities and functionals are written out at the same times that the solution
:math:`q` is written.  While these could be computed in postprocessing, it is more efficient
to compute them at run-time for large parallel runs.  

Gauge output is written at every timestep.  In order to get this data without a
gauge, one would otherwise have to write the full solution out at every
timestep, which might be very slow.


Outputting derived quantities
===============================
It is sometimes desirable to output quantities other than those
in the vector q.  To do so, just add a function `p_function` to 
the controller that accepts the state and sets the derived quantities
in state.p

.. testsetup:: *

    import pyclaw
    import numpy as np
    solver = pyclaw.ClawSolver2D()
    x = pyclaw.Dimension('x',-1.0,1.0,100)
    y = pyclaw.Dimension('y',-1.0,1.0,100)
    domain = pyclaw.Domain([x,y])
    num_eqn = 3; num_aux = 2 
    state = pyclaw.State(domain,num_eqn,num_aux)
    grid = state.grid
    Y,X = np.meshgrid(grid.y.centers,grid.x.centers)
    r = np.sqrt(X**2 + Y**2)
    width = 0.2
    state.q[0,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.
    claw = pyclaw.Controller()

.. doctest::

    >>> def stress(state):
    ...    state.p[0,:,:] = np.exp(state.q[0,:,:]*state.aux[1,:,:]) - 1.

    >>> state.mp = 1
    >>> claw.p_function = stress

Outputting functionals
===============================
In PyClaw a functional is a scalar quantity computed from :math:`q` that is written
to file at each output time.  For now, only functionals of the form

.. math::
   \begin{equation}
	F(q) = \int |f(q)| dx dy
   \end{equation}	

are supported.  In other words, the functional must be the absolute
integral of some function of :math:`q`.  To enable writing functionals, simply
set state.mf to the number of functionals and point the controller to a 
function that computes :math:`f(q)`

.. doctest::

    >>> def compute_f(state):
    ...    state.F[0,:,:] = state.q[0,:,:]*state.q[1,:,:]

    >>> state.mf = 1

Using gauges
===================
A gauge in PyClaw is a single grid location for which output is written at
every time step.  This can be very useful in some applications, like comparing
with data from tidal gauges (from whence the name is derived) in tsunami modeling.
The gauges are managed by the grid object, and a grid at location :math:`(x,y)` 
may be added simply by calling `grid.add_gauges((x,y))`.  Multiple gauges
can be set at once by providing a list of coordinate tuples

.. testsetup::

    x1 = 0; x2 = 1; x3 = 2
    y1 = 1; y2 = 2; y3 = 0

.. doctest::

    >>> state.grid.add_gauges([(x1,y1),(x2,y2),(x3,y3)])

By default, the solution values are written out at each gauge location.
To write some other quantity, simply provide a function 
:math:`f(q,aux)` and point the solver to it

.. doctest::

    >>> def f(q,aux):
    ...    return q[1,:,:]/q[0,:,:]

    >>> solver.compute_gauge_values = f
