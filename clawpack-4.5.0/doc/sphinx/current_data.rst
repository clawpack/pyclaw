

.. _current_data:

***************
current_data
***************

The Python plotting routines often allow the user to specify a call back
function as plotting parameters, for example :ref:`ClawPlotAxes` has an
attribute *afteraxes* that can be set to a function to be executed after
performing all plots on the specified axes.  This is useful for setting
parameters for these axes that are not covered by the provided attributes,
for annotating the plots on these axes, for adding a plot of the true
solution, etc.

Call back functions include:

 * *beforeframe* and *afterframe* attributes of :ref:`ClawPlotData` 
 * *afteraxes* attribute of :ref:`ClawPlotAxes` 
 * *afteritem*, *aftergrid*, *plot_var*, *map2d_to_1d* attributes of :ref:`ClawPlotItem` 


All of these functions are designed to take a single argument
*current_data*, an object with attributes that may be useful to the user in
evaluating the function.  

.. warning:: :ref:`mapc2p` is one exception that does not take argument *current_data*.



Attributes of *current_data*:
-----------------------------

Some of these may be unavailable because they don't make sense in the
current context, e.g. in a *beforeframe* function.

.. attribute:: plotdata : 

    parent :ref:`ClawPlotData` object.  From this object virtually any
    relevant data can be accessed.  Other attributes are defined for
    convenience.  If you find you frequently
    need something else, you could modify :mod:`pyclaw.plotters.frametools`
    to add this to *current_data*.


.. attribute:: frameno : 

   The current frame number

.. attribute:: grid : 

   Object of :class:`pyclaw.solution.grid` with data for the last grid
   plotted.

.. attribute:: gridno : 

   Grid number of this grid, of interest only in AMR calculations.

.. attribute:: q : 

    q array for current frame, so for example the in a scalar 2d problem the
    value in the (i,j) cell would be *current_data.q[i,j,0]* (remember that
    Python always indexes starting at 0).

    In an AMR calculation q will be from the last grid plotted.  

.. attribute:: aux : 

    aux array for current frame, provided these have been output by the
    Fortran code.  At the moment this requires modifying the library routine 
    `outN.f` to set outaux = .true.  so that fort.a files are produced along
    with fort.q files.  [This should be an input parameter!]

    If fort.a files are not found then current_data.aux will be None.

    In an AMR calculation aux will be from the last grid plotted.  

.. attribute:: var : 

    array of the variable actually plotted most recently, e.g. if
    *plotitem.plot_var == 0* then in 2d *current_data.var[i,j] ==
    current_data.q[i,j,0]*.

.. attribute:: level : 

   For AMR computations, where *current_data.grid* is for the last grid plotted,
   *current_data.level* is the AMR level of this grid.  Particularly useful
   in `aftergrid` functions.

.. attribute:: t : 

    the current time t.

.. attribute:: x : 

    x array of cell centers corresponding to *current_data.q*.  

.. attribute:: y : 

    y array of cell centers corresponding to *current_data.q* (only in 2d).  

.. attribute:: xlower :

    left edge of current grid.

.. attribute:: ylower :

    left edge of current grid in y (only in 2d).



