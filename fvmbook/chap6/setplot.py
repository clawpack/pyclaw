
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
"""
from advection import beta, x0, IC

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    plotdata.clearfigures()  # clear any old figures,axes,items data



    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes(name='Solution')
    plotaxes.xlimits = 'auto'
    if IC=='gauss_square':
        plotaxes.ylimits = [-0.5, 1.5]
    elif IC=='wavepacket':
        plotaxes.ylimits = [-1.0, 1.5]
    plotaxes.title = 'q[0]'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(name='solution', plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    plotaxes.afteraxes = plot_true_soln
    plotitem.show = True       # show on plot?
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via visclaw.frametools.printframes:

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


#-------------------
def plot_true_soln(current_data):
#-------------------
    from numpy import linspace, mod, exp, sin
    from pylab import plot
    xtrue = linspace(0.,1.,1000)
    t = current_data.t
    xshift = xtrue - t
    # periodic boundary conditions
    xshift = mod(xshift, 1.0)
    if IC=='gauss_square':
        x1 = 0.6; x2 = 0.8
        qtrue = exp(-beta * (xshift-x0)**2) + (xshift>0.6)*(xshift<0.8)
    elif IC=='wavepacket':
        qtrue = exp(-beta * (xshift-x0)**2) * sin(80.*xshift)
    plot(xtrue, qtrue, 'r')
