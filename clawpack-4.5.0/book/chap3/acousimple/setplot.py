
"""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.

"""


#--------------------------
def setplot(plotdata):
#--------------------------
    """   
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """

    from pyclaw.plotters.frametools import var_limits


    plimits = [-1., 1.]
    ulimits = [-1., 1.]
    xlimits = 'auto'          # choose automatically


    # Pressure:
    # ---------
    
    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=1)
    plotaxes = plotfigure.new_plotaxes(name='Pressure')
    plotaxes.axescmd = 'subplot(1,1,1)' 
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = plimits
    plotitem = plotaxes.new_plotitem(name='Pressure',plot_type='1d')
    plotitem.plot_var = 0       # q[0] is the pressure
    plotitem.plotstyle = '-'
    plotitem.color = 'b'


    # Velocity:
    # ---------
    plotfigure = plotdata.new_plotfigure(name='Velocity', figno=2)
    plotaxes = plotfigure.new_plotaxes(name='Velocity')
    plotaxes.axescmd = 'subplot(1,1,1)' 
    plotaxes.xlimits = xlimits
    plotaxes.ylimits = ulimits
    plotitem = plotaxes.new_plotitem(name='Velocity',plot_type='1d')
    plotitem.plot_var = 1       # q[1] is the velocity
    plotitem.plotstyle = '-'
    plotitem.color = 'b'


    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

