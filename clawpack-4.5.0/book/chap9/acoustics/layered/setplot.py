
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
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    plotdata.clearfigures()  # clear any old figures,axes,items data


    def shade_layers(current_data):
        import pylab
        for i in range(0,121,2):
            pylab.fill([i,i,i+1,i+1],[-.1,.8,.8,-.1],fc=[.8,.8,.8])

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes(name='Solution')
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.xlimits = [0,120]
    plotaxes.ylimits = [-.1, .8]
    plotaxes.title = 'Pressure'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(name='solution', plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    plotaxes.afteraxes = shade_layers

    # Zoomed view
    plotaxes = plotfigure.new_plotaxes(name='q0zoom')
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = [20,40]
    plotaxes.ylimits = [-.1, .8]
    plotaxes.title = 'zoom of pressure'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(name='solution', plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    plotaxes.afteraxes = shade_layers
    

    # Figure for q[1]
    plotfigure = plotdata.new_plotfigure(name='Velocity', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes(name='Solution')
    plotaxes.xlimits = [0,120]
    plotaxes.ylimits = [-.1, .8]
    plotaxes.title = 'Velocity'
    plotaxes.axescmd = 'subplot(211)'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(name='solution', plot_type='1d')
    plotitem.plot_var = 1
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    plotaxes.afteraxes = shade_layers

    # Zoomed view
    plotaxes = plotfigure.new_plotaxes(name='q1zoom')
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.xlimits = [20,40]
    plotaxes.ylimits = [-.1, .8]
    plotaxes.title = 'zoom of velocity'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(name='solution', plot_type='1d')
    plotitem.plot_var = 1
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    plotaxes.afteraxes = shade_layers

    
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

    
