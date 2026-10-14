""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 
import numpy as np


#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    
    # Figure for pcolor plot
    plotfigure = plotdata.new_plotfigure(name='Water height', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-3.0, 1.0]
    plotaxes.ylimits = [-1.0, 1.0]
    plotaxes.title = 'Water height'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    plotitem.plot_var = 0
    plotitem.contour_levels = np.linspace(-1.59729945e-03, 1.59729945e-03, 200)
    plotitem.contour_colors = 'k'
    plotitem.patchedges_show = 1

    plotitem.show = True       # show on plot?

    return plotdata
