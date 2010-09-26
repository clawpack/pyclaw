
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from pyclaw.data import Data
userdata = Data(['claw.data','setprob.data'])

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    plotdata.clearfigures()  # clear any old figures,axes,items data

    

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0., 1.]
    plotaxes.ylimits = [-1., 1.5]
    plotaxes.title = 'q[0]'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    plotaxes.afteraxes = plot_true_soln
    

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    

#-------------------
def plot_true_soln(current_data):
#-------------------
    from numpy import linspace, mod, exp, where, sin
    from pylab import plot
    xtrue = linspace(userdata.xlower,userdata.xupper,1000)
    t = current_data.t
    xshift = xtrue - userdata.u*t
    if userdata.mthbc_xupper==2:
        # for periodic boundary conditions
        xshift = userdata.xlower + mod(xshift-userdata.xlower, \
                 userdata.xupper-userdata.xlower)
    x0 = 0.5
    qtrue = exp(-userdata.beta*(xshift - x0)**2) * sin(userdata.freq*xshift)
    plot(xtrue, qtrue, 'r')

