

.. _setrun:

*****************************************************************
Specifying run-time parameters in `setrun.py`
*****************************************************************

This section explains the parameters needed for the classic single-grid
Clawpack code.  Additional parameters are needed by extensions of the code.
For these, see:

 * AMRClaw (adaptive mesh refinement): :ref:`setrun_amrclaw`

 * GeoClaw (geophysical flows): :ref:`setrun_geoclaw`

*Describe setrun in general here!*

It may be useful to look at a specific example, e.g. 
:ref:`setrun_sample`.

Run-time parameters
-------------------

The parameters needed in 1 space dimension (*ndim=1*) are described.  In 2d
and 3d there are analogous parameters in y and z required, as mentioned
below.

.. attribute:: ndim : integer from [1,2,3]

   number of space dimensions.  

.. attribute:: xlower : float

   lower limit in the x direction.   
   For 2d or 3d, there are analogous parameters *ylower* and *zlower*.

.. attribute:: xupper : float

   upper limit in the x direction.   
   For 2d or 3d, there are analogous parameters *yupper* and *zupper*.

.. attribute:: mx : integer

   The number of grid cells in the x direction.
   For 2d or 3d, there are analogous parameters *my* and *mz*.

.. attribute:: meqn : integer

   Number of equations in the system (e.g. *meqn=1* for a scalar problem).

.. attribute:: maux : integer

   Number of auxiliary variables in the aux array (initialized in setaux.f)

.. attribute:: mcapa : integer

   Index of aux array corresponding to capacity function, if there is one.

.. attribute:: t0 : float

   Initial time, often *t0 = 0.*

.. attribute:: dt_initial: float

   Initial time step to try in first step.

.. attribute:: dt_variable: boolean

   If True, time steps are adjusted automatically based on the desired
   Courant number *cfl_desired*.  

   If False, fixed time steps of lenght *dt_initial* are used.

.. attribute:: dt_max: float

   If *dt_variable = True* then this is an upper bound on the allowable time
   step regardless of the Courant number.  Useful if there are other reasons
   to limit the time step (e.g. stiff source terms).

.. attribute:: cfl_desired: float

   If *dt_variable = True* then this is the desired Courant number.  Time
   steps will be adjusted based on the maximum wave speed seen in the *last*
   time step taken.  For a nonlinear problem this may not result in the
   Courant number being exactly the desired value in the next step, which is
   where the next attribute comes in...

   Usually *cfl_desired = 0.9* or less.

.. attribute:: cfl_max: float

   If *dt_variable = True* then this is the maximum Courant number that can
   be allowed.  If a time step results in a Courant number that's greater
   than *cfl_desired* but less than or equale to *cfl_max*, the step is
   accepted.  If the Courant number is greater than *cfl_max* then the step
   is rejected and a smaller step is taken.  (At this point the maximum wave
   speed from Riemann solutions is known, so the step can be adjusted to
   exactly hit the desired value *cfl_desired*.)

   Usually *cfl_max = 1.0* is fine.
   
.. attribute:: max_step: int

   Maximum number of time steps allowed between output times.  This is just
   to avoid infinite loops and generally a large value is fine.


.. attribute:: outstyle : integer

   There are three possible ways to specify the output
   times.   This parameter selects the desired manner to specify the times,
   and affects what other attributes are required.

     * *outstyle = 1* Requires *nout* and *tfinal*, the number of output
       times to produce up to time *tfinal*.  They will be equally spaced.
       The time steps will be adjusted to hit these times exactly. (Provided
       *dt_variable = True*.  Otherwise *dt_initial* must divide
       *tfinal/nout* an integer number of times.)

     * *outstyle = 2* Requires *nout* and *tout*, where *tout* is a list
       of *nout* desired output times.

     * *outstyle = 3* Requires *nout* and *iout*, and the solution is output 
       every *iout* time steps for a total of *nout* steps.  
       

To be continued...
