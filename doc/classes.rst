.. _pyclaw_classes:
  
*****************************************
Understanding Pyclaw Classes
*****************************************
.. contents::

Flow of a Pyclaw Simulation
===========================

The basic idea of a pyclaw simulation is to construct a
:class:`~pyclaw.solution.Solution` object, hand it to a
:class:`~pyclaw.solver.Solver` object, and request a solution at a new
time.  The solver will take whatever steps are necessary to evolve the solution
to the requested time.

.. image:: images/pyclaw_architecture_flow.*

The bulk of the work in order to run a simulation then is the creation and
setup of the appropriate :class:`~pyclaw.geometry.Domain`, :class:`~pyclaw.state.State`,
:class:`~pyclaw.solution.Solution`, and :class:`~pyclaw.solver.Solver`
objects needed to evolve the solution to the requested time.

Creation of a Pyclaw :class:`~pyclaw.solution.Solution`
=======================================================

.. image:: images/pyclaw_solution_structure.*

A Pyclaw :class:`~pyclaw.solution.Solution` is a container for a collection of
:class:`~pyclaw.geometry.Domain` and :class:`~pyclaw.state.State` designed with a 
view to future support of adaptive mesh refinement and multi-block simulations. The :class:`~pyclaw.solution.Solution` 
object keeps track of a list of :class:`~pyclaw.state.State` objects
and controls the overall input and output of the entire collection of 
:class:`~pyclaw.state.State` objects.  Each
:class:`~pyclaw.state.State` object inhabits a :class:`~pyclaw.geometry.Grid`, composed of
:class:`~pyclaw.geometry.Dimension` objects that define the extents 
of the :class:`~pyclaw.grid.Domain`.  Multiple states can inhabit the same
grid, but each :class:`~pyclaw.state.State` inhabits a single grid.

The process needed to create a :class:`~pyclaw.solution.Solution` object then
follows from the bottom up.

.. doctest::

    >>> from pyclaw import Solution, State, Dimension, Domain
    >>> x = Dimension('x', -1.0, 1.0, 200)
    >>> y = Dimension('y', 0.0, 1.0, 100)
    
This code creates two dimensions, a dimension ``x``  on the interval 
``[-1.0, 1.0]`` with :math:`200` grid points and a dimension ``y`` on the interval
``[0.0, 1.0]`` with :math:`100` grid points.  

.. note:: 

    Many of the attributes of a :class:`~pyclaw.geometry.Dimension`
    object are set automatically so make sure that the values you want are set
    by default.  Please refer to the :class:`~pyclaw.geometry.Dimension`
    classes definition for what the default values are.

Next we have to create a :class:`~pyclaw.geometry.Domain` object that will
contain our :class:`~pyclaw.geometry.Domain.dimensions` objects.

.. doctest::

    >>> grid = Domain([x,y])
    >>> num_eqn = 2
    >>> state = State(grid,num_eqn)


Here we create a ``grid`` with the dimensions we created earlier to make a single 2D 
:class:`~pyclaw.geometry.Domain` object. Then we set the number of equations the State 
will represent to 2. Finally, we create a :class:`~pyclaw.state.State` that inhabits 
this grid. As before, many of the attributes of the :class:`~pyclaw.geometry.Domain` 
and State objects are set automatically.

We now need to set the initial condition ``q`` and possibly ``aux`` to the correct
values.

.. doctest::

    >>> import numpy as np
    >>> sigma = 0.2
    >>> omega = np.pi
    >>> Y,X = np.meshgrid(state.grid.y.centers,state.grid.x.centers)
    >>> r = np.sqrt(X**2 + Y**2)
    >>> state.q[0,:] = np.cos(omega * r)
    >>> state.q[1,:] = np.exp(-r**2 / sigma**2)
    
We now have initialized the first entry of ``q`` to a cosine function 
evaluated at the cell centers and the second entry of ``q`` to a gaussian, again
evaluated at the grid cell centers.

Many Riemann solvers also require information about the problem we are going
to run which happen to be grid properties such as the impedence :math:`Z` and 
speed of sound :math:`c` for linear acoustics.  We can set these values in the 
``problem_data`` dictionary in one of two ways.  The first way is to set them
directly as in:

.. doctest::

    >>> state.problem_data['c'] = 1.0
    >>> state.problem_data['Z'] = 0.25
    
If you're using a Fortran Riemann solver, these values will automatically get
copied to the corresponding variables in the cparam common block of the
Riemann solver.  This is done in solver.setup(), which calls state.set_cparam().

Last we have to put our :class:`~pyclaw.state.State` object into a 
:class:`~pyclaw.solution.Solution` object to complete the process.  In this
case, since we are not using adaptive mesh refinement or a multi-block
algorithm, we do not have multiple grids.

.. doctest::

    >>> sol = Solution(state,grid)
    
We now have a solution ready to be evolved in a 
:class:`~pyclaw.solver.Solver` object.


Creation of a Pyclaw :class:`~pyclaw.solver.Solver`
==========================================================

A Pyclaw :class:`~pyclaw.solver.Solver` can represent many different
types of solvers; here we will use a 1D, classic Clawpack type of
solver.  This solver is defined in the :mod:`~pyclaw.classic.solver` module.

First we import the particular solver we want and create it with the default 
configuration.

.. doctest::

    >>> from pyclaw import ClawSolver1D, BC
    >>> solver = ClawSolver1D()
    >>> solver.bc_lower[0] = BC.periodic
    >>> solver.bc_upper[0] = BC.periodic

Next we need to tell the solver which Riemann solver to use from the
:ref:`pyclaw_rp`. We can always 
check what Riemann solvers are available to use via the :mod:`~pyclaw.riemann` 
module. Once we have picked one out, we pass it to the solver via:

.. doctest::

    >>> from pyclaw import riemann 
    >>> solver.rp = riemann.rp_acoustics.rp_acoustics_1d

In this case we have decided to use the 1D linear acoustics Riemann solver.  You 
can also set your own solver by importing the module that contains it and 
setting it directly to the `rp` attribute of the particular object in the class 
:class:`~pyclaw.classic.solver.ClawSolver1D`.

.. doctest::

    >>> import my_rp_module # doctest: +SKIP
    >>> solver.rp = my_rp_module.my_acoustics_rp # doctest: +SKIP

Last we finish up by specifying solver options, if we want to override the
defaults.  For instance, we might want to specify a particular limiter

.. doctest::

    >>> from pyclaw import limiters
    >>> solver.limiters = limiters.tvd.vanleer
    
If we wanted to control the simulation we could at this point by issuing the 
following commands:

.. doctest::

    >>> solver.evolve_to_time(sol,1.0) # doctest: +SKIP

This would evolve our solution ``sol`` to ``t = 1.0`` but we are then
responsible for all output and other setup considerations.

Creating and Running a Simulation with :class:`~pyclaw.controller.Controller`
=============================================================================

The :class:`~pyclaw.controller.Controller` coordinates the output and setup of
a run with the same parameters as the classic Clawpack.  In order to have it 
control a run, we need only to create the controller, assign it a solver and
initial condition, and call the :meth:`~pyclaw.controller.Controller.run`
method.

.. testsetup::

    import pyclaw
    x = pyclaw.Dimension('x',0.0,1.0,100)
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,2)
    sol = pyclaw.Solution(state,domain)

.. doctest::

    >>> from pyclaw.controller import Controller

    >>> claw = Controller()
    >>> claw.solver = solver
    >>> claw.solutions = sol
    
Here we have imported and created the :class:`~pyclaw.controller.Controller` 
class, assigned the :class:`~pyclaw.solver.Solver` and 
:class:`~pyclaw.solution.Solution`.

These next commands setup the type of output the controller will output.  The
parameters are similar to the ones found in the classic clawpack claw.data 
format.

.. doctest::

    >>> claw.output_style = 1
    >>> claw.num_output_times = 10
    >>> claw.tfinal = 1.0
    
When we are ready to run the simulation, we can call the 
:meth:`~pyclaw.controller.Controller.run` method.  It will then run the
simulation and output the appropriate time points.  If the 
:attr:`~pyclaw.controller.Controller.keep_copy` is set to *True* the 
controller will keep a copy of each solution output in memory in the frames array. 
For instance, you can then immediately plot the solutions output into the *frames*
array.


Restarting a simulation
=========================
To restart a simulation, simply initialize a Solution object using an output
frame from a previous run; for example, to restart from frame 3

.. doctest::

    >>> claw.solution = Solution(3, file_format='petsc')

By default, the :class:`~pyclaw.controller.Controller` will number your
output frames starting from the frame number used for initializing
the :class:`~pyclaw.solution.Solution` object.
If you want to change the default behaviour and start counting frames
from zero, you will need to pass the keyword argument
``count_from_zero=True`` to the solution initializer.


.. note::
    
    It is necessary to specify the output format ('petsc' or 'ascii').
    

If your simulation includes aux variables, you will need to either recompute them or
output the aux values at every step, following the instructions below.


Outputting aux values
===============================
To write aux values to disk at the initial time::

    >>> claw.write_aux_init = True

To write aux values at every step::

    >>> claw.write_aux_always = True

Outputting derived quantities
===============================
It is sometimes desirable to output quantities other than those
in the vector q.  To do so, just add a function `compute_p` to 
the controller that accepts the state and sets the derived quantities
in state.p

.. doctest::

    >>> def stress(state):
    ...     state.p[0,:,:] = np.exp(state.q[0,:,:]*state.aux[1,:,:]) - 1.
 
    >>> state.mp = 1
    >>> claw.compute_p = stress

