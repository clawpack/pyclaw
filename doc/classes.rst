.. _pyclaw_classes:
  
*****************************************
Pyclaw Classes
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
setup of the appropriate :class:`~pyclaw.grid.Grid`, :class:`~pyclaw.state.State`,
:class:`~pyclaw.solution.Solution`, and :class:`~pyclaw.solver.Solver`
objects needed to evolve the solution to the requested time.

Creation of a Pyclaw :class:`~pyclaw.solution.Solution`
=======================================================

.. image:: images/pyclaw_solution_structure.pdf

A Pyclaw :class:`~pyclaw.solution.Solution` is a container for a collection of
:class:`~pyclaw.grid.Grid` and :class:`~pyclaw.state.State` designed with a 
view to future support of adaptive mesh 
refinement and multi-block simulations. The :class:`~pyclaw.solution.Solution` 
object keeps track of a list of :class:`~pyclaw.state.State` objects
and controls the overall input and output of the entire collection of 
:class:`~pyclaw.state.State` objects.  Each
:class:`~pyclaw.state.State` object inhabits a `~pyclaw.grid.Grid`, composed of
:class:`~pyclaw.grid.Dimension` objects that define the extents 
of the :class:`~pyclaw.grid.Grid`.  Multiple states can inhabit the same
grid, but each :class:`~pyclaw.state.State` inhabits a single grid.

The process needed to create a :class:`~pyclaw.solution.Solution` object then
follows from the bottom up.

::

    >>> from pyclaw import Solution, State, Grid, Dimension
    
    >>> x = Dimension('x', -1.0, 1.0, 200)
    >>> y = Dimension('y', 0.0, 1.0, 100)
    
This code creates two dimensions, a dimension ``x``  on the interval 
``[-1.0, 1.0]`` with ``200`` grid points and a dimension ``y`` on the interval
``[0.0, 1.0]`` with ``100`` grid points.  

.. note:: 

    Many of the attributes of a :class:`~pyclaw.grid.Dimension`
    object are set automatically so make sure that the values you want are set
    by default.  Please refer to the :class:`~pyclaw.grid.Dimension`
    classes definition for what the default values are.

Next we have to create a :class:`~pyclaw.grid.Grid` object that will
contain our :class:`~pyclaw.grid.Dimension` objects.

::

    >>> grid = Grid([x,y])
    >>> state = State(grid)
    >>> state.num_eqn = 2

Here we create a grid with the dimensions we created earlier to make a single
2D :class:`~pyclaw.grid.Grid` object.  Then we create a `~pyclaw.state.State`
that inhabits this Grid. Finally, we set the number of equations the State
will represent to 2.  As before, many of the attributes of the
:class:`~pyclaw.grid.Grid` and State objects are set automatically.

We now need to set the initial condition ``q`` and possibly ``aux`` to the correct
values.  

::

    >>> import numpy as np
    >>> sigma = 0.2
    >>> omega = np.pi
    >>> state.q[:,0] = np.cos(omega * grid.x.center)
    >>> state.q[:,1] = np.exp(-grid.x.center**2 / sigma**2)
    
We now have initialized the first entry of q to a cosine function 
evaluated at the cell centers and the second entry of q to a gaussian, again
evaluated at the grid cell centers.

Many Riemann solvers also require information about the problem we are going
to run which happen to be grid properties such as the impedence ``Z`` and 
speed of sound ``c`` for linear acoustics.  We can set these values in the 
``aux_global`` dictionary in one of two ways.  The first way is to set them
directly as in:

::

    >>> state.aux_global['c'] = 1.0
    >>> state.aux_global[`Z`] = 0.25
    
If you're using a Fortran Riemann solver, these values will automatically get
copied to the corresponding variables in the cparam common block of the
Riemann solver.  This is done in solver.setup(), which calls grid.set_cparam().

Last we have to put our :class:`~pyclaw.state.State` object into a 
:class:`~pyclaw.solution.Solution` object to complete the process.  In this
case, since we are not using adaptive mesh refinement or a multi-block
algorithm, we do not have multiple grids.

::

    >>> sol = Solution(state)
    
We now have a solution ready to be evolved in a 
:class:`~pyclaw.solver.Solver` object.


Creation of a Pyclaw :class:`~pyclaw.solver.Solver`
==========================================================

A Pyclaw :class:`~pyclaw.solver.Solver` can represent many different
types of solvers; here we will use a 1D, classic Clawpack type of
solver.  This solver is defined in the :mod:`~pyclaw.clawpack` module.

First we import the particular solver we want and create it with the default 
configuration.

::

    >>> from pyclaw.clawpack import ClawSolver1D
    >>> solver = ClawSolver1D()
    >>> solver.bc_lower[0] = pyclaw.BC.periodic
    >>> solver.bc_upper[0] = pyclaw.BC.periodic

Next we need to tell the solver which Riemann solver to use from the
:doc:`Riemann solver package </pyclaw/evolve/rp>` .  We can always check what 
Riemann solvers are available to use via the 
:meth:`~pyclaw.ClawSolver1D.list_riemann_solvers` method.  Once we have
picked one out, we let the solver pick it out for us via:

::

    >>> solver.set_riemann_solver('acoustics')

In this case we have decided to use the linear acoustics Riemann solver.  You 
can also set your own solver by importing the module that contains it and 
setting it directly to the :attr:`~pyclaw.clawpack.ClawSolver1D.rp`
attribute to the particular function.

::

    >>> import my_rp_module
    >>> solver.rp = my_rp_module.my_acoustics_rp

Last we finish up by specifying solver options, if we want to override the
defaults.  For instance, we might want to specify a particular limiter::

    >>> solver.limiters = pyclaw.limiters.vanleer
    
If we wanted to control the simulation we could at this point by issuing the 
following commands:

::

    >>> solver.evolve_to_time(sol,1.0)
    
This would evolve our solution ``sol`` to ``t = 1.0`` but we are then
responsible for all output and other setup considerations.

Creating and Running a Simulation with :class:`~pyclaw.controller.Controller`
=============================================================================

The :class:`~pyclaw.controller.Controller` coordinates the output and setup of
a run with the same parameters as the classic Clawpack.  In order to have it 
control a run, we need only to create the controller, assign it a solver and
initial condition, and call the :meth:`~pyclaw.controller.Controller.run`
method.

::

    >>> from pyclaw.controller import Controller

    >>> claw = Controller()
    >>> claw.solver = solver
    >>> claw.solutions['n'] = sol
    
Here we have imported and created the :class:`~pyclaw.controller.Controller` 
class, assigned the :class:`~pyclaw.solver.Solver` and 
:class:`~pyclaw.solution.Solution`.

These next commands setup the type of output the controller will output.  The
parameters are similar to the ones found in the classic clawpack claw.data 
format.

::

    >>> claw.outstyle = 1
    >>> claw.nout = 10
    >>> claw.tfinal = 1.0
    
When we are ready to run the simulation, we can call the 
:meth:`~pyclaw.controller.Controller.run` method.  It will then run the
simulation and output the appropriate time points.  If the 
:attr:`~pyclaw.controller.Controller.keep_copy` is set to *True* the 
controller will keep a copy of each solution output in memory in the frames array.  For
instance, you can then immediately plot the solutions output into the *frames*
array.


Restarting a simulation
=========================
To restart a simulation, simply initialize a Solution object using an output
frame from a previous run; for example, to restart from frame 3::

    >>> claw.solution = pyclaw.Solution(3,format='petsc')

.. note::
    
    It is necessary to specify the output format ('petsc' or 'ascii').


Outputting derived quantities
===============================
It is sometimes desirable to output quantities other than those
in the vector q.  To do so, just add a function `compute_p` to 
the controller that accepts the state and sets the derived quantities
in state.p::

    >>> state.mp = 1
    >>> claw.compute_p = stress
    >>> def stress(state):
    >>>     state.p[0,:,:] = np.exp(state.q[0,:,:]*state.aux[1,:,:]) - 1.

