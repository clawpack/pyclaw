
""" 
Pseudo-color (pcolor) plot
==========================

Produce a pcolor plot of pressure from 2d acoustics example.

With AMR data, show the grid lines on the two coarsest levels only.
    
""" 

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from pyclaw.plotters import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    

    plotfigure = plotdata.new_plotfigure(name='Pressure', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0., 1.]
    plotaxes.ylimits = [0., 1.]
    plotaxes.title = 'Pressure'
    plotaxes.scaled = True        # so aspect ratio is 1

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.pcolor_cmin = 0.   # p=0 mapped to low end (yellow)
    plotitem.pcolor_cmax = 2.   # p=2 mapped to high end (blue)
    plotitem.add_colorbar = True
    plotitem.amr_gridlines_show = [1,1,0]   # show grid lines on 2 coarsest levels
    plotitem.amr_gridedges_show = [1]       # show grid edges on all levels
    

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 1           # layout of plots
    plotdata.latex_framesperline = 2         # layout of plots
    plotdata.latex_makepdf = True            # also run pdflatex?

    return plotdata

    
