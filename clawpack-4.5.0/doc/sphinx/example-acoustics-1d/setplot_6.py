
""" 
Specifying a function plot_var for the plotting variable
========================================================

The pressure q[0] and q[1] are plotted on two sets of axes in a single
figure. 

The Riemann invariants R1 = (-p+Z*u)/2*Z and R2 = (p-Z*u)/2*Z
are plotted in a second figure.

This illustrates how to specify plot_var as a function.

""" 

from numpy import sqrt
from pyclaw.data import Data
probdata = Data('setprob.data')      # read problem data values
Z = sqrt(probdata.rho * probdata.K)  # impedance used in R1 and R2


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    plotdata.clearfigures()  # clear any old figures,axes,items data

    plotfigure = plotdata.new_plotfigure(name='Solution', figno=1)

    # Pressure:

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-.5,1.1]
    plotaxes.title = 'Pressure'

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    
    # Velocity:

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-.5,.5]
    plotaxes.title = 'Velocity'

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 1
    plotitem.plotstyle = 'o-'
    plotitem.color = 'b'

    # Riemann invariants:
    # -------------------

    plotfigure = plotdata.new_plotfigure(name='Riemann invariants', figno=2)

    # R1

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-.5,0.6]
    plotaxes.title = 'Left-going Riemann invariant'

    plotitem = plotaxes.new_plotitem(plot_type='1d')

    def R1(current_data):
        q = current_data.q   # solution when this function called
        p = q[:,0]           # pressure
        u = q[:,1]           # velocity
        R1 = (-p + Z*u) / (2.*Z)  # uses impedance Z set above
        return R1

    plotitem.plot_var = R1  # defined above
    plotitem.plotstyle = '-o'
    plotitem.color = 'r'
    
    # R2

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-.5,0.6]
    plotaxes.title = 'Right-going Riemann invariant'

    plotitem = plotaxes.new_plotitem(plot_type='1d')

    def R2(current_data):
        q = current_data.q   # solution when this function called
        p = q[:,0]           # pressure
        u = q[:,1]           # velocity
        R2 = (p + Z*u) / (2.*Z)  # uses impedance Z set above
        return R2

    plotitem.plot_var = R2  # defined above
    plotitem.plotstyle = '-o'
    plotitem.color = 'g'
    


    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'# pointer for index page
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 1           # layout of plots
    plotdata.latex_framesperline = 2         # layout of plots
    plotdata.latex_makepdf = True            # also run pdflatex?

    return plotdata
    
