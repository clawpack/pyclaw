
""" 
Automatically compute ylimits from data
=======================================

This example illustrates how to determine the ylimits that should be set on
each axis based on the data.  

One could use::
 
  plotaxes.ylimits = 'auto'

for each variable to have them automatically set by matplotlib, but then
they would be reset for each frame as the solution evolves.

Instead, this example uses a function `pyclaw.plotters.frametools.var_limits`
that reads in all the data files and determines limits for a set of
specified variables.
 
""" 

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 
    from pyclaw.plotters.frametools import var_limits


    plotdata.clearfigures()  # clear any old figures,axes,items data

    vars = [0,1]   # variables to be plotted
    use_global_limits = True
    if use_global_limits:
        # scan all frames and to determine limits
        # padding parameter adds a bit to varmin and varmax for plot limits.
        varmin,varmax,varlim = var_limits(plotdata,vars,padding=0.1)
	plimits = varlim[0]
	ulimits = varlim[1]
    else:
	plimits = 'auto'   # or set to known limits
	ulimits = 'auto'


    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Solution', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = plimits
    plotaxes.title = 'Pressure'

    # Set up item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'o-'
    plotitem.color = 'b'

    
    # Figure for q[1]
    plotfigure = plotdata.new_plotfigure(name='Velocity', figno=2)

    # Set axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = ulimits
    plotaxes.title = 'Velocity'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 1
    plotitem.plotstyle = 'o-'
    plotitem.color = 'r'

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'# pointer for index page
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = True            # also run pdflatex?

    return plotdata
    
