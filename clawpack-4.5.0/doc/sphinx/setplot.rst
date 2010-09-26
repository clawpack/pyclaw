.. _setplot:

*********************************************
Using setplot.py to specify the desired plots 
*********************************************


The desired plots are specified by creating an object
of class ClawPlotData that contains specifications of what *figures* should be
created, and within each figure what sets of *axes* should be drawn, and
within each axes what *items* should be plotted (lines, contour plots,
etc.).



Plotting Data Objects
=====================

More details about each class of objects can be found
on these pages:

.. toctree::
   :maxdepth: 2

   ClawPlotData
   ClawPlotFigure
   ClawPlotAxes
   ClawPlotItem


For examples, see :ref:`plotexamples`.

.. _setplot_overview:

Overview
========

The approach outlined below may seem more complicated than necessary, and it
would be if all you ever want to do is plot one set of data at each output
time.  However, when adaptive mesh refinement is used each frame of data may
contain several grids and so creating the desired plot requires looping over
all grids.  This is done by the plotting utilities described in :ref:`plotting`,
but for this to work it is necessary to specify what plot(s) are desired.

Most example directories contain a file setplot.py that contains a
function setplot(plotdata). This function
sets various attributes of plotdata
in order to specify how plotting is to be done.

The object plotdata is of class ClawPlotData.  The way to set up the plot
structure is to follow this outline:


 * Specify some attributes of setplot that determine  what sort of plots
   will be produced and where they will be stored, e.g.::

      plotdata.plotdir = '_plots'

   will cause hardcopy to go to subdirectory _plots of the current directory.
   (Attributes like plotdir that are only used for hardcopy are often set in
   the script plotclaw.py rather than in setplot.  See :ref:`plot_files`.)

   There are many other :ref:`ClawPlotData` attributes and methods.

 * Specify one or more Figures to be created for each frame, e.g.::
 
      plotfigure = plotdata.new_plotfigure(name='Solution', figno=1)
 
   plotfigure is now an object of class ClawPlotFigure and various attributes
   can be set, e.g.:: 
   
      plotfigure.kwargs = {'figsize':[8,12], 'facecolor':'#ff9999'}
   
   to specify any keyword arguments 
   that should be used when creating this figure in matplotlib.
   The above would create a figure that is
   8 inches by 12 inches with a pink background.  


   For more options,
   see the `matplotlib documentation <http://matplotlib.sourceforge.net/>`_
   for the `figure command
   <http://matplotlib.sourceforge.net/api/figure_api.html#matplotlib.figure.Figure>`_.

   There are many other :ref:`plotfigure` attributes and methods.

 * Specify one or more Axes to be created within each figure, e.g.::

    plotaxes = plotfigure.new_plotaxes()

   Note that new_plotaxes is a method of class ClawPlotFigure and creates a set
   of axes specific to the particular object plotfigure.

   plotaxes is now an object of class ClawPlotAxes  and various attributes
   can be set, e.g.::

      plotfigure.ylimits = [-1, 1]


   There are many other :ref:`ClawPlotAxes` attributes and methods.


 *  Specify one or more Items to be created within these axes, e.g.::

      plotitem = plotaxes.new_plotitem(plot_type='1d_plot')

    Note that new_plotitem is a method of class ClawPlotAxes and creates a plot
    object (e.g. line, contour plot, etc)
    specific to the particular object plotaxes.  

    plotitem is now an object of class ClawPlotItem  and various attributes
    can be set, e.g.::

      plotitem.plotstyle = '-'
      plotitem.color = 'r'

    for a solid line that is red.

    There are many other :ref:`ClawPlotItem` attributes and methods.


