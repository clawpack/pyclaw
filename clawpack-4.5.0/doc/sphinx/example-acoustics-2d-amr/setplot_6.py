
""" 
Scatter plot together with 1d reference solution
================================================

Produce a scatter plot of the pressure in cell (i,j) vs. r(i,j) = distance
from origin.  Add to this a plot of a 1d reference solution obtained by
solving the 1d radially symmetric problem.  This output is obtained from a
subdirectory 1drad/_output of the 2d directory.  See
`<claw/doc/sphinx/example-acoustics-2d/1drad/README.html>`_
    
This example adapted from book/chap21/acoustics.
""" 

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 

    import os

    plotdata.clearfigures()  # clear any old figures,axes,items data
    
    plotfigure = plotdata.new_plotfigure(name='scatter', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = [-0.2,2.5]
    plotaxes.title = 'Scatter plot'

    # Scatter of 2d data
    # ------------------
    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')
    
    def p_vs_r(current_data):
        # Return radius of each grid cell and p value in the cell
        from pylab import sqrt
        x = current_data.x
        y = current_data.y
        r = sqrt(x**2 + y**2)
        q = current_data.q
        p = q[:,:,0]
        return r,p

    plotitem.map_2d_to_1d = p_vs_r
    plotitem.plot_var = 0
    plotitem.amr_plotstyle = ['o', '^', '+']   # symbol for each AMR level
    plotitem.amr_color = ['r','k','b']    # color for each AMR level
    plotitem.show = True       # show on plot?
    
    
    # 1d reference solution
    # ---------------------
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.outdir = os.path.join(os.getcwd(),'1drad/_output')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-'
    plotitem.color = 'r'
    plotitem.kwargs = {'linewidth': 2}
    

    # ---------------------
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

    
