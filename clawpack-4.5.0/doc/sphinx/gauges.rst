

.. _gauges:

***************
Gauges
***************

With AMRCLAW in two space dimensions
it is possible to specify gauge locations as points (x,y) where the values of all
components of q should be output every time step during the computation over some
time range (t1,t2).  The addition of a time range is new in svn revision 328.

Gauges are useful in several ways, e.g.:

 1. To compare computational results to measurements from 
    physical gauges such as a pressure gauge or tide gauge that
    record data as a function of time at a single point,

 2. To better visualize how the solution behaves at a single point,

 3. To better compare results obtained with different methods or grid resolutions.
    Comparing two-dimensional pcolor or contour plots can be difficult whereas
    comparing to curves that give the solution as a function of time often reveals
    more clearly differences in accuracy or nonphysical oscillations.

To use gauges, include the line::

    call setgauges()

in your file setprob.f.  This reads in the gauge locations from the file 
`setgauges.data`.  This file should have one line giving the number of gauges and
the following lines specify information for each gauge in the format::

    gaugeno  x  y  t1  t2

Rather than creating this file by hand, it is easiest to do this in your
`setrun.py` file by including lines of the form::

    clawdata.gauges = []
    clawdata.gauges.append([gaugeno, x, y, t1, t2])

where the second line is repeated as many times as desired with different values of
the parameters.  From this information the `setgauges.data` file will be
automatically created when you do "make .data" or "make .output".

During the computation the value of all components of q at all gauge locations will
be output to a single file `fort.gauge` in the output directory.  Lines of this
file have the form::

   gaugeno  level  t  q[0]  q[1] ...  q[meqn-1]

where level is the AMR level used to determine the q values at this time.
Internally the finest level available at each gauge is used, with bilinear
interpolation to the gauge locations from the 4 nearest cell centers.

If you wish to change what is output at these points, you should copy the library
routine `dumpgauge.f` to your own directory and modify it appropriately.

Plotting tools
--------------

Several Python plotting tools are available to plot the gauge data, so you do not
have to parse the file `fort.gauge` yourself.  In the `setplot` Python script you
can specify plots that are to be done for each gauge, similar to the manner in
which you can specify plots that are to be done for each time frame.  For example,
to plot the component q[0] at each gauge, include in `setplot` lines of this form::

    plotfigure = plotdata.new_plotfigure(name='q[0] at gauges', figno=300, \
                    type='each_gauge')

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-1.5, 1.5]
    plotaxes.title = 'q[0]'

    # Plot q[0] as blue line:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'b-'

Note that `plotdata.new_plotfigure` is called with `type='each_gauge'` which
denotes that this plot is to be produced for each gauge found in `setgauges.data`.
(When type is not specified, the default is `type='each_frame'` for time frame data).

If you type::

    $ make .plots

then html files will be created for the gauge plots along with the time frame plots
and will all show up in the index (usually in `_plots/_PlotIndex.html`).

When using Iplotclaw to interactively view plots, try::

    PLOTCLAW> plotgauge 1

to produce the plot for gauge 1, or simply::

    PLOTCLAW> plotgauge 

to loop through all gauges.

You can of course specify more than one plotitem on each plotaxes if you want.  For
example to plot the each gauge from the current run as a blue line and the same
gauge from some previous run (perhaps with a different grid resolution)
as a red line, you could add the following lines to the above example::

    # Plot q[0] from previous run as red line:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'r-'
    plotitem.outdir = '_output_from_previous_run'


Plotting gauge locations
------------------------

It is often convenient to plot the locations of the gauges on pcolor or contour
plots each time frame.  You can do this as follows, for example::

    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    # set other attributes as desired

    def addgauges(current_data):
        from pyclaw.plotters import gaugetools
        gaugetools.plot_gauge_locations(current_data.plotdata, \
             gaugenos='all', format_string='ko', add_labels=True)

    plotaxes.afteraxes = addgauges

You can replace `gaugenos='all'` by `gaugenos=[1,2]` or other list of specific
gauges to plot.  The `format_string` above specifies a black dot at each gauge
location and `add_labels=True` means that the gauge number will appear next to each
gauge.

If you want more control over this plotting you can of course copy the function
`plot_gauge_locations` from `pyclaw.plotters.gaugetools.py` 
to your setplot.py file and modify at will.

Examples
--------

To see an example of the use of gauges see:

   * `[$CLAW/apps/advection/2d/example1_gauges/amr]
     <claw/apps/advection/2d/example1_gauges/amr/README.html>`_ 

