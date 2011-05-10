

"""
Module frametools for plotting frames of time-dependent data.

"""

import os,sys,shutil,glob
import string,re
import time
import traceback


from pyclaw.data import Data
from pyclaw.plotters import plotpages
from matplotlib.colors import Normalize 

plotter = 'matplotlib'
if plotter == 'matplotlib':
    if not sys.modules.has_key('matplotlib'):
        try:
            import matplotlib
            matplotlib.use('Agg')  # Use an image backend
        except:
            print "*** Error: problem importing matplotlib"

try:
    import pylab
except:
    print "*** Error: problem importing pylab"


try:
    from numpy import ma
except:
    print "*** Error: problem importing masked array module ma"
    


#==========================================
def plotframe(frameno, plotdata, verbose=False):
#==========================================

    """
    Plot a single frame from the computation.

    This routine checks input and then calls plotframeN for the
    proper space dimension.

    """

    if verbose:  print '    Plotting frame %s ... '  % frameno

    if plotdata.mode() == 'iplotclaw':
        pylab.ion()
        
    try:
        plotfigure_dict = plotdata.plotfigure_dict
    except:
        print '*** Error in plotframe: plotdata missing plotfigure_dict'
        print '*** This should not happen'
        return None

    if len(plotfigure_dict) == 0:
        print '*** Warning in plotframe: plotdata has empty plotfigure_dict'
        print '*** Apparently no figures to plot'


    try:
        framesoln = plotdata.getframe(frameno, plotdata.outdir)
    except:
        print '*** Cannot find frame number ',frameno
        print '*** looking in directory ', plotdata.outdir
        return None

    t = framesoln.t

    # initialize current_data containing data that will be passed
    # to afterframe, afteraxes, aftergrid commands
    current_data = Data()
    current_data.user = Data()   # for user specified attributes
                                 # to avoid potential conflicts
    current_data.plotdata = plotdata
    current_data.frameno = frameno
    current_data.t = t


    # call beforeframe if present, which might define additional 
    # attributes in current_data or otherwise set up plotting for this
    # frame.

    beforeframe =  getattr(plotdata, 'beforeframe', None)
    if beforeframe:
        if isinstance(beforeframe, str):
            # a string to be executed
            exec(beforeframe)
        else:
            # assume it's a function
            try:
                output = beforeframe(current_data)
                if output: current_data = output
            except:
                print '*** Error in beforeframe ***'
                raise



    # iterate over each single plot that makes up this frame:
    # -------------------------------------------------------
 
    if plotdata._mode == 'iplotclaw':
        print '    Plotting Frame %s at t = %s' % (frameno,t)
        requested_fignos = plotdata.iplotclaw_fignos
    else:
        requested_fignos = plotdata.print_fignos
    plotted_fignos = []

    plotdata = set_show(plotdata)   # set _show attributes for which figures
                                    # and axes should be shown.

    plotdata = set_outdirs(plotdata)  # set _outdirs attribute to be list of
                                      # all outdirs for all items

    # loop over figures to appear for this frame: 
    # -------------------------------------------

    for figname in plotdata._fignames:
        plotfigure = plotdata.plotfigure_dict[figname]
        if (not plotfigure._show) or (plotfigure.type != 'each_frame'):
            continue  # skip to next figure 

        figno = plotfigure.figno
        #print '+++ Figure: ',figname,figno
        if requested_fignos != 'all':
            if figno not in requested_fignos:
                continue # skip to next figure

        plotted_fignos.append(figno)


        if not plotfigure.kwargs.has_key('facecolor'):
            # use Clawpack's default bg color (tan)
            plotfigure.kwargs['facecolor'] = '#ffeebb'   

        # create figure and set handle:
        plotfigure._handle = pylab.figure(num=figno, **plotfigure.kwargs)

        pylab.ioff()
        if plotfigure.clf_each_frame:
            pylab.clf()


        try:
            plotaxes_dict = plotfigure.plotaxes_dict
        except:
            print '*** Error in plotframe: plotdata missing plotaxes_dict'
            print '*** This should not happen'
            return  None

        if (len(plotaxes_dict) == 0) or (len(plotfigure._axesnames) == 0):
            print '*** Warning in plotframe: plotdata has empty plotaxes_dict'
            print '*** Apparently no axes to plot in figno ',figno

        # loop over axes to appear on this figure:
        # ----------------------------------------

        for axesname in plotfigure._axesnames:
            plotaxes = plotaxes_dict[axesname]
            if not plotaxes._show:
                continue   # skip this axes if no items show

            # create the axes:
            axescmd = getattr(plotaxes,'axescmd','subplot(1,1,1)')
            axescmd = 'plotaxes._handle = pylab.%s' % axescmd
            exec(axescmd)
            pylab.hold(True)


            # NOTE: This was rearranged December 2009 to 
            # loop over grids first and then over plotitems so that 
            # finer grids will plot over coarser grids of other items.
            # Needed in particular if masked arrays are used, e.g. if
            # two different items do pcolor plots with different color maps
            # for different parts of the domain (e.g. water and land).

            # Added loop over outdirs to prevent problems if different
            # items use data from different outdirs since loop over items
            # is now inside loop on grids.


            # loop over all outdirs:

            for outdir in plotdata._outdirs:
                try:
                    framesoln = plotdata.getframe(frameno, outdir)
                except:
                    print '*** Cannot find frame number ',frameno
                    print '*** looking in directory ', outdir
                    return current_data


                if framesoln.t != t:
                    print '*** Warning: t values do not agree for frame ',frameno
                    print '*** t = %g for outdir = %s' % (t,plotdata.outdir)
                    print '*** t = %g for outdir = %s' % (framesoln.t,outdir)

                current_data.framesoln = framesoln

                #print "+++ Looping over grids in outdir = ",outdir

                # loop over grids:
                # ----------------
    
                for gridno in range(len(framesoln.grids)):
                    #print '+++ gridno = ',gridno
                    grid = framesoln.grids[gridno]
                    current_data.grid = grid
                    current_data.q = grid.q
                    current_data.aux = grid.aux
                    current_data.xlower = grid.dimensions[0].lower
                    current_data.xupper = grid.dimensions[0].upper
    
                    # loop over items:
                    # ----------------
    
                    for itemname in plotaxes._itemnames:
                        
                        plotitem = plotaxes.plotitem_dict[itemname]
                        #print '+++ %s: %s' % (itemname,plotitem.outdir)

                        #import pdb; pdb.set_trace()
                        item_outdir = plotitem.outdir
                        if not plotitem.outdir:
                            item_outdir = plotdata.outdir
                        if item_outdir != outdir:
                            # skip to next item
                            #print '+++ skipping, plotitem.outdir=',item_outdir
                            #print '+++           plotdata.outdir=',plotdata.outdir
                            #print '+++                    outdir=',outdir
                            continue
                            

                        ndim = plotitem.ndim

                        # option to suppress printing some levels:
                        try:
                            pp_amr_data_show = plotitem.amr_data_show
                            i = min(len(pp_amr_data_show), grid.level) - 1
                            show_this_level = pp_amr_data_show[i]
                        except:
                            show_this_level = True

                        if plotitem._show and show_this_level:
                            cmd = 'output = plotitem%s(framesoln,plotitem,\
                                    current_data,gridno)'  % ndim
    
                            try:
                                exec(cmd)
                                if output: current_data = output
                                if verbose:  
                                        print '      Plotted  plotitem ', itemname
                            except:
                                print '*** Error in plotframe: problem calling plotitem%s' % ndim
                                traceback.print_exc()
                                return None
    
                    # end of loop over plotitems
                # end of loop over grids
            # end of loop over outdirs


            for itemname in plotaxes._itemnames:
                plotitem = plotaxes.plotitem_dict[itemname]
                if plotitem.afteritem:
                    print "*** ClawPlotItem.afteritem is deprecated"
                    print "*** use ClawPlotAxes.afteraxes "
                    print "*** or  ClawPlotItem.aftergrid instead"
                try:
                    if plotitem.add_colorbar:
                        pobj = plotitem._current_pobj # most recent plot object
                        pylab.colorbar(pobj)
                except:
                    pass


            if plotaxes.title_with_t:
                if (t==0.) | ((t>=0.001) & (t<1000.)):
                    pylab.title("%s at time t = %14.8f" % (plotaxes.title,t))
                else:
                    pylab.title("%s at time t = %14.8e" % (plotaxes.title,t))
            else:
                pylab.title(plotaxes.title)


            # call an afteraxes function if present:
            afteraxes =  getattr(plotaxes, 'afteraxes', None)
            if afteraxes:
                if isinstance(afteraxes, str):
                    # a string to be executed
                    exec(afteraxes)
                else:
                    # assume it's a function
                    try:
                        current_data.plotaxes = plotaxes
                        current_data.plotfigure = plotaxes._plotfigure
                        output = afteraxes(current_data)
                        if output: current_data = output
                    except:
                        print '*** Error in afteraxes ***'
                        raise

            if plotaxes.scaled:
                pylab.axis('scaled')

            # set axes limits:
            if (plotaxes.xlimits is not None) & (type(plotaxes.xlimits) is not str):
                try:
                    pylab.xlim(plotaxes.xlimits[0], plotaxes.xlimits[1])
                except:
                    pass  # let axis be set automatically
            if (plotaxes.ylimits is not None) & (type(plotaxes.ylimits) is not str):
                try:
                    pylab.ylim(plotaxes.ylimits[0], plotaxes.ylimits[1])
                except:
                    pass  # let axis be set automatically



            # end of loop over plotaxes
            
        # end of loop over plotfigures


    # call an afterframe function if present:
    afterframe =  getattr(plotdata, 'afterframe', None)
    if afterframe:
        if isinstance(afterframe, str):
            # a string to be executed
            exec(afterframe)
        else:
            # assume it's a function
            try:
                output = afterframe(current_data)
                if output: current_data = output
            except:
                print '*** Error in afterframe ***'
                raise


    if plotdata.mode() == 'iplotclaw':
        pylab.ion()
    for figno in plotted_fignos:
        pylab.figure(figno)
        pylab.draw()

    if verbose:
        print '    Done with plotframe for frame %i at time %g' % (frameno,t)

    
    # print the figure(s) to file(s) if requested:
    if (plotdata.mode() != 'iplotclaw') & plotdata.printfigs:
        # iterate over all figures that are to be printed:
        for figno in plotted_fignos:
            printfig(frameno=frameno, figno=figno, \
                    format=plotdata.print_format, plotdir=plotdata.plotdir,\
                    verbose=verbose)

    return current_data

    # end of plotframe



    
#==================================================================
def plotitem1(framesoln, plotitem, current_data, gridno):
#==================================================================
    """
    Make a 1d plot for a single plot item for the solution in framesoln.

    The current_data object holds data that should be passed into
    aftergrid or afteraxes if these functions are defined.  The functions
    may add to this object, so this function should return the possibly
    modified current_data for use in other plotitems or in afteraxes or
    afterframe.

    """

    plotdata = plotitem._plotdata
    plotfigure = plotitem._plotfigure
    plotaxes = plotitem._plotaxes

    # the following plot parameters should be set and independent of 
    # which AMR level a grid is on:

    plot_params = """
             plot_type  afteritem  mapc2p  MappedGrid gaugeno
             """.split()

    # Convert each plot parameter into a local variable starting with 'pp_'.

    for plot_param in plot_params:
        cmd = "pp_%s = getattr(plotitem, '%s',None)" \
             % (plot_param, plot_param)
        exec(cmd)


    if pp_mapc2p is None:
        # if this item does not have a mapping, check for a global mapping:
        pp_mapc2p = getattr(plotdata, 'mapc2p', None)

    grid = framesoln.grids[gridno]
    current_data.grid = grid
    current_data.q = grid.q
    current_data.aux = grid.aux

    t = framesoln.t

    # the following plot parameters may be set, depending on what
    # plot_type was requested:

    plot_params = """
             plot_var  aftergrid  plotstyle color kwargs 
             plot_var2 fill_where map_2d_to_1d plot_show
             """.split()

    # For each plot_param check if there is an amr_ list for this
    # parameter and if so select the value corresponding to the level of
    # this grid.  Otherwise, use plotitem attribute of this name.
    # The resulting variable starts with 'pp_'.

    for plot_param in plot_params:
        amr_plot_param = "amr_%s" % plot_param
        amr_list = getattr(plotitem, amr_plot_param, [])
        if len(amr_list) > 0:
            index = min(grid.level, len(amr_list)) - 1
            exec("pp_%s = amr_list[%i]" % (plot_param, index))
        else:
            exec("pp_%s = getattr(plotitem, '%s', None)" \
                 % (plot_param, plot_param))

    if pp_plot_type == '1d':
        pp_plot_type = '1d_plot'  # '1d' is deprecated
    
    if pp_plot_type in ['1d_plot', '1d_semilogy']:
        thisgridvar = get_gridvar(grid,pp_plot_var,1,current_data)
        xc_center = thisgridvar.xc_center   # cell centers
        xc_edge = thisgridvar.xc_edge       # cell edges
        var = thisgridvar.var               # variable to be plotted
        current_data.x = xc_center
        current_data.var = var

    elif pp_plot_type == '1d_fill_between':
        try: pylab.fill_between
        except: 
            print "*** This version of pylab is missing fill_between"
            print "*** Reverting to 1d_plot"
            pp_plot_type = '1d_plot'
        thisgridvar = get_gridvar(grid,pp_plot_var,1,current_data)
        thisgridvar2 = get_gridvar(grid,pp_plot_var2,1,current_data)
        xc_center = thisgridvar.xc_center   # cell centers
        xc_edge = thisgridvar.xc_edge       # cell edges
        var = thisgridvar.var               # variable to be plotted
        var2 = thisgridvar2.var             # variable to be plotted
	current_data.x = xc_center
	current_data.var = var
	current_data.var2 = var2


    elif pp_plot_type == '1d_from_2d_data':
        if not pp_map_2d_to_1d:
            print '*** Error, plot_type = 1d_from_2d_data requires '
            print '*** map_2d_to_1d function as plotitem attribute'
            raise
            return
        try:
            thisgridvar = get_gridvar(grid,pp_plot_var,2,current_data)
            xc_center = thisgridvar.xc_center   # cell centers
            yc_center = thisgridvar.yc_center
            xc_edge = thisgridvar.xc_edge   # cell edge
            yc_edge = thisgridvar.yc_edge
            var = thisgridvar.var             # variable to be plotted
            current_data.x = xc_center
            current_data.y = yc_center
            current_data.var = var
            # TODO: fix for xp, yp??
            xc_center, var = pp_map_2d_to_1d(current_data)
            xc_center = xc_center.flatten()  # convert to 1d
            var = var.flatten()  # convert to 1d

            #xc_center, var = pp_map_2d_to_1d(var,xc_center,yc_center,t)
            #xc_edge, var = pp_map_2d_to_1d(var,xc_edge,yc_edge,t)
        except:
            print '*** Error with map_2d_to_1d function'
            print 'map_2d_to_1d = ',pp_map_2d_to_1d
            raise
            return

    elif pp_plot_type == '1d_gauge_trace':
        gaugesoln = plotdata.getgauge(pp_gaugeno)
        xc_center = None
        xc_edge = None
        

    elif pp_plot_type == '1d_empty':
        pp_plot_var = 0  # shouldn't be used but needed below *FIX*
        thisgridvar = get_gridvar(grid,pp_plot_var,1,current_data)
        xc_center = thisgridvar.xc_center   # cell centers
        xc_edge = thisgridvar.xc_edge       # cell edges

    else:
        raise ValueError("Unrecognized plot_type: %s" % pp_plot_type)
        return None


    # Grid mapping:

    if pp_MappedGrid is None:
        pp_MappedGrid = (pp_mapc2p is not None)

    if (pp_MappedGrid & (pp_mapc2p is None)):
        print "*** Warning: MappedGrid == True but no mapc2p specified"
    elif pp_MappedGrid:
        X_center= pp_mapc2p(xc_center)
        X_edge = pp_mapc2p(xc_edge)
    else:
        X_center = xc_center
        X_edge = xc_edge


    # The plot commands using matplotlib:

    pylab.hold(True)


    if (pp_plot_type in ['1d_plot','1d_from_2d_data']):
        if pp_color:
            pp_kwargs['color'] = pp_color

        plotcommand = "pobj=pylab.plot(X_center,var,'%s', **pp_kwargs)"  \
                      % pp_plotstyle
        if pp_plot_show:
            exec(plotcommand)

    elif pp_plot_type == '1d_semilogy':
        if pp_color:
            pp_kwargs['color'] = pp_color

        plotcommand = "pobj=pylab.semilogy(X_center,var,'%s', **pp_kwargs)"  \
                      % pp_plotstyle
        if pp_plot_show:
            exec(plotcommand)

    elif pp_plot_type == '1d_fill_between':
        if pp_color:
            pp_kwargs['color'] = pp_color
        if pp_fill_where:
            pp_fill_where = pp_fill_where.replace('plot_var','var')
            pp_fill_where = pp_fill_where.replace('plot_var2','var2')
            plotcommand = \
              "pobj=pylab.fill_between(X_center,var,var2,%s,**pp_kwargs)"  \
                  % pp_fill_where

        else:
            plotcommand = "pobj=pylab.fill_between(X_center,var,var2,**pp_kwargs)" 
        if pp_plot_show:
            exec(plotcommand)

    elif pp_plot_type == '1d_gauge_trace':

        gauget = gaugesoln.t
        gaugeq = gaugesoln.q[3,:]
        plotcommand = "pobj=pylab.plot(gauget, gaugeq)"
        if pp_plot_show:
            exec(plotcommand)

        # interpolate to the current time t:
        try:
            i1 = pylab.find(gauget < t)[-1]
            i1 = min(i1,len(gauget)-2)
            slope = (gaugeq[i1+1]-gaugeq[i1]) / (gauget[i1+1]-gauget[i1])
            qt = gaugeq[i1] + slope * (t-gauget[i1])
        except:
            qt = gaugeq[0]
            print "Warning: t out of range"
        pylab.plot([t], [qt], 'ro')

        

    elif pp_plot_type == '1d_empty':
        # no plot to create (user might make one in afteritem or
        # afteraxes)
        pass

    else:
        raise ValueError("Unrecognized plot_type: %s" % pp_plot_type)
        return None


    # call an aftergrid function if present:
    if pp_aftergrid:
        if isinstance(pp_aftergrid, str):
            # a string to be executed
            exec(pp_aftergrid)
        else:
            # assume it's a function
            try:
                # set values that may be needed in aftergrid:
                current_data.gridno = gridno
                current_data.plotitem = plotitem
                current_data.grid = grid
                current_data.var = var
                current_data.xlower = grid.dimensions[0].lower
                current_data.xupper = grid.dimensions[0].upper
                current_data.x = X_center # cell centers
		current_data.dx = grid.d[0]
                output = pp_aftergrid(current_data)
                if output: current_data = output
            except:
                print '*** Error in aftergrid ***'
                raise


    try:
        plotitem._current_pobj = pobj
    except:
        pass # if no plot was done



    return current_data
    

#==================================================================
def plotitem2(framesoln, plotitem, current_data, gridno):
#==================================================================
    """
    Make a 2d plot for a single plot item for the solution in framesoln.

    The current_data object holds data that should be passed into
    aftergrid or afteraxes if these functions are defined.  The functions
    may add to this object, so this function should return the possibly
    modified current_data for use in other plotitems or in afteraxes or
    afterframe.

    """

    import numpy as np
    from pyclaw.plotters import colormaps

    plotdata = plotitem._plotdata
    plotfigure = plotitem._plotfigure
    plotaxes = plotitem._plotaxes

    # The following plot parameters should be set and independent of 
    # which AMR level a grid is on:

    plot_params = """
             plot_type  afteritem  mapc2p  MappedGrid
             """.split()

    # Convert each plot parameter into a local variable starting with 'pp_'.

    for plot_param in plot_params:
        cmd = "pp_%s = getattr(plotitem, '%s',None)" \
             % (plot_param, plot_param)
        exec(cmd)

    if pp_mapc2p is None:
        # if this item does not have a mapping, check for a global mapping:
        pp_mapc2p = getattr(plotdata, 'mapc2p', None)

    grid = framesoln.grids[gridno]
    level = grid.level

    current_data.grid = grid
    current_data.q = grid.q
    current_data.aux = grid.aux
    current_data.level = level

    t = framesoln.t

    # The following plot parameters may be set, depending on what
    # plot_type was requested:

    plot_params = """
             plot_var  aftergrid  kwargs 
             gridlines_show  gridlines_color  grid_bgcolor
             gridedges_show  gridedges_color  add_colorbar
             pcolor_cmap  pcolor_cmin  pcolor_cmax
             imshow_cmap  imshow_cmin  imshow_cmax
             contour_levels  contour_nlevels  contour_min  contour_max
             contour_colors   contour_cmap  contour_show
             schlieren_cmap  schlieren_cmin schlieren_cmax
             quiver_coarsening  quiver_var_x  quiver_var_y  quiver_key_show
             quiver_key_scale  quiver_key_label_x  quiver_key_label_y  
             quiver_key_scale  quiver_key_units  quiver_key_kwargs
             """.split()

    # For each plot_param check if there is an amr_ list for this
    # parameter and if so select the value corresponding to the level of
    # this grid.  Otherwise, use plotitem attribute of this name.
    # The resulting variable starts with 'pp_'.

    for plot_param in plot_params:
        amr_plot_param = "amr_%s" % plot_param
        amr_list = getattr(plotitem, amr_plot_param, [])
        if len(amr_list) > 0:
            index = min(grid.level, len(amr_list)) - 1
            exec("pp_%s = amr_list[%i]" % (plot_param, index))
        else:
            exec("pp_%s = getattr(plotitem, '%s', None)" \
                 % (plot_param, plot_param))



    # turn grid background color into a colormap for use with pcolor cmd:
    pp_grid_bgcolormap = colormaps.make_colormap({0.: pp_grid_bgcolor, \
                                             1.: pp_grid_bgcolor})
    
    thisgridvar = get_gridvar(grid,pp_plot_var,2,current_data)

    xc_center = thisgridvar.xc_center   # cell centers (on mapped grid)
    yc_center = thisgridvar.yc_center
    xc_edge = thisgridvar.xc_edge       # cell edges (on mapped grid)
    yc_edge = thisgridvar.yc_edge
    var = thisgridvar.var             # variable to be plotted


    # Grid mapping:

    if pp_MappedGrid is None:
        pp_MappedGrid = (pp_mapc2p is not None)

    if (pp_MappedGrid & (pp_mapc2p is None)):
        print "*** Warning: MappedGrid == True but no mapc2p specified"
    elif pp_MappedGrid:
        X_center, Y_center = pp_mapc2p(xc_center, yc_center)
        X_edge, Y_edge = pp_mapc2p(xc_edge, yc_edge)
    else:
        X_center, Y_center = xc_center, yc_center
        X_edge, Y_edge = xc_edge, yc_edge


    # The plot commands using matplotlib:

    pylab.hold(True)


    if ma.isMaskedArray(var):
        # If var is a masked array: plotting should work ok unless all 
        # values are masked, in which case pcolor complains and there's
        # no need to try to plot.  Check for this case...
        var_all_masked = (ma.count(var) == 0)
    else:
        # not a masked array, so certainly not all masked:
        var_all_masked = False

    if pp_plot_type == '2d_pcolor':

        pcolor_cmd = "pobj = pylab.pcolor(X_edge, Y_edge, var, \
                        cmap=pp_pcolor_cmap"

        if pp_gridlines_show:
            pcolor_cmd += ", edgecolors=pp_gridlines_color"
        else: 
            pcolor_cmd += ", shading='flat'"

        pcolor_cmd += ", **pp_kwargs)"

        if not var_all_masked:
            exec(pcolor_cmd)


            if (pp_pcolor_cmin not in ['auto',None]) and \
                     (pp_pcolor_cmax not in ['auto',None]):
                pylab.clim(pp_pcolor_cmin, pp_pcolor_cmax) 
        else:
            #print '*** Not doing pcolor on totally masked array'
            pass

    elif pp_plot_type == '2d_imshow':

        if not var_all_masked:
            if pp_imshow_cmin in ['auto',None]:
                pp_imshow_cmin = np.min(var)
            if pp_imshow_cmax in ['auto',None]:
                pp_imshow_cmax = np.max(var)
            color_norm = Normalize(pp_imshow_cmin,pp_imshow_cmax,clip=True)

            xylimits = (X_edge[0,0],X_edge[-1,-1],Y_edge[0,0],Y_edge[-1,-1])
            pobj = pylab.imshow(pylab.flipud(var.T), extent=xylimits, \
                    cmap=pp_imshow_cmap, interpolation='nearest', \
                    norm=color_norm)

            if pp_gridlines_show:
                # This draws grid for labels shown.  Levels not shown will
                # not have lower levels blanked out however.  There doesn't
                # seem to be an easy way to do this. 
                pobj = pylab.plot(X_edge, Y_edge, color=pp_gridlines_color)
                pobj = pylab.plot(X_edge.T, Y_edge.T, color=pp_gridlines_color)

        else:
            #print '*** Not doing imshow on totally masked array'
            pass


    elif pp_plot_type == '2d_contour':
        levels_set = True
        if pp_contour_levels is None:
            levels_set = False
            if pp_contour_nlevels is None:
                print '*** Error in plotitem2:'
                print '    contour_levels or contour_nlevels must be set'
                raise
                return
            if (pp_contour_min is not None) and \
                    (pp_contour_max is not None):

                pp_contour_levels = pylab.linspace(pp_contour_min, \
                       pp_contour_max, pp_contour_nlevels)
                levels_set = True 


        if pp_gridlines_show:
            pobj = pylab.pcolor(X_edge, Y_edge, pylab.zeros(var.shape), \
                    cmap=pp_grid_bgcolormap, edgecolors=pp_gridlines_color)
        elif pp_grid_bgcolor is not 'w': 
            pobj = pylab.pcolor(X_edge, Y_edge, pylab.zeros(var.shape), \
                    cmap=pp_grid_bgcolormap, edgecolors='None')
        pylab.hold(True)


        # create the contour command:
        contourcmd = "pobj = pylab.contour(X_center, Y_center, var, "
        if levels_set:
            contourcmd += "pp_contour_levels"
        else:
            contourcmd += "pp_contour_nlevels"

        if pp_contour_cmap:
            if (pp_kwargs is None) or ('cmap' not in pp_kwargs):
                contourcmd += ", cmap = pp_contour_cmap"
        elif pp_contour_colors:
            if (pp_kwargs is None) or ('colors' not in pp_kwargs):
                contourcmd += ", colors = pp_contour_colors"

        contourcmd += ", **pp_kwargs)"

        if (pp_contour_show and not var_all_masked):
            # may suppress plotting at coarse levels
            exec(contourcmd)

    elif pp_plot_type == '2d_grid':
        # plot only the grids, no data:
        if pp_gridlines_show:
            pobj = pylab.pcolor(X_edge, Y_edge, pylab.zeros(var.shape), \
                    cmap=pp_grid_bgcolormap, edgecolors=pp_gridlines_color,\
                    shading='faceted')
        else: 
            pobj = pylab.pcolor(X_edge, Y_edge, pylab.zeros(var.shape), \
                    cmap=pp_grid_bgcolormap, shading='flat')


    elif pp_plot_type == '2d_schlieren':
        # plot 2-norm of gradient of variable var:
        (vx,vy) = pylab.gradient(var)
        vs = pylab.sqrt(vx**2 + vy**2)

        pcolor_cmd = "pobj = pylab.pcolor(X_edge, Y_edge, vs, \
                        cmap=pp_schlieren_cmap"

        if pp_gridlines_show:
            pcolor_cmd += ", edgecolors=pp_gridlines_color"
        else: 
            pcolor_cmd += ", edgecolors='None'"

        pcolor_cmd += ", **pp_kwargs)"

        if not var_all_masked:
            exec(pcolor_cmd)

            if (pp_schlieren_cmin not in ['auto',None]) and \
                     (pp_schlieren_cmax not in ['auto',None]):
                pylab.clim(pp_schlieren_cmin, pp_schlieren_cmax) 

    elif pp_plot_type == '2d_quiver':
        if pp_quiver_coarsening > 0:
            var_x = get_gridvar(grid,pp_quiver_var_x,2,current_data).var
            var_y = get_gridvar(grid,pp_quiver_var_y,2,current_data).var
            Q = pylab.quiver(X_center[::pp_quiver_coarsening,::pp_quiver_coarsening],
                             Y_center[::pp_quiver_coarsening,::pp_quiver_coarsening],
                             var_x[::pp_quiver_coarsening,::pp_quiver_coarsening],
                             var_y[::pp_quiver_coarsening,::pp_quiver_coarsening],
                             **pp_kwargs)
                             # units=pp_quiver_units,scale=pp_quiver_scale)

            # Show key
            if pp_quiver_key_show:
                if pp_quiver_key_scale is None:
                    key_scale = np.max(np.sqrt(var_x**2+var_y**2))*0.5
                else:
                    key_scale = pp_quiver_key_scale
                label = r"%s %s" % (str(np.ceil(key_scale)),pp_quiver_key_units)
                pylab.quiverkey(Q,pp_quiver_key_label_x,pp_quiver_key_label_y,
                                key_scale,label,**pp_quiver_key_kwargs)

    elif pp_plot_type == '2d_empty':
        # no plot to create (user might make one in afteritem or
        # afteraxes)
        pass

    else:
        raise ValueError("Unrecognized plot_type: %s" % pp_plot_type)
        return None

    # end of various plot types



    # plot grid patch edges if desired:

    if pp_gridedges_show:
        for i in [0, X_edge.shape[0]-1]:
            X1 = X_edge[i,:]
            Y1 = Y_edge[i,:]
            pylab.plot(X1, Y1, pp_gridedges_color)
        for i in [0, X_edge.shape[1]-1]:
            X1 = X_edge[:,i] 
            Y1 = Y_edge[:,i]
            pylab.plot(X1, Y1, pp_gridedges_color)


    if pp_aftergrid:
        try:
            if isinstance(pp_aftergrid, str):
                exec(pp_aftergrid)
            else:
                # assume it's a function
                current_data.gridno = gridno
                current_data.plotitem = plotitem
                current_data.grid = grid
                current_data.var = var
                current_data.xlower = grid.dimensions[0].lower
                current_data.xupper = grid.dimensions[0].upper
                current_data.ylower = grid.dimensions[0].lower
                current_data.yupper = grid.dimensions[0].upper
                current_data.x = X_center # cell centers
                current_data.y = Y_center # cell centers
                current_data.dx = grid.d[0]
                current_data.dy = grid.d[1]

                output = pp_aftergrid(current_data)
                if output: current_data = output
        except:
            print '*** Warning: could not execute aftergrid'
            raise


    try:
        plotitem._current_pobj = pobj
    except:
        pass # if no plot was done


    return current_data


#--------------------------------------
def get_gridvar(grid, plot_var, ndim, current_data):
#--------------------------------------
    """
    Return arrays for spatial variable(s) (on mapped grid if necessary)
    and variable to be plotted on a single grid in ndim space dimensions.
    """

    #grid.compute_physical_coordinates()
    #grid.compute_computational_coordinates()
    
    if ndim == 1:
    
        # +++ until bug in solution.py fixed.
        xc_center = grid.c_center[0]
        xc_edge = grid.c_edge[0]
        #xc_center = grid.p_center[0]
        #xc_edge = grid.p_edge[0]
        current_data.x = xc_center
        current_data.dx = grid.d[0]

    
        if isinstance(plot_var, int):
            var = grid.q[plot_var,:]
        else:
            try:
                #var = plot_var(grid.q, xc_center, grid.t)
                var = plot_var(current_data)
            except:
                print '*** Error applying function plot_var = ',plot_var
                traceback.print_exc()
                return  

        thisgridvar = Data()
        thisgridvar.var = var
        thisgridvar.xc_center = xc_center
        thisgridvar.xc_edge = xc_edge
        thisgridvar.mx = grid.dimensions[0].n

        # end of 1d case
    
    elif ndim == 2:
        # +++ until bug in solution.py fixed.
        #xc_center, yc_center = grid.c_center
        #xc_edge, yc_edge = grid.c_edge
        xc_center, yc_center = grid.p_center
        xc_edge, yc_edge = grid.p_edge
        current_data.x = xc_center
        current_data.y = yc_center
        current_data.dx = grid.d[0]
        current_data.dy = grid.d[1]

        if isinstance(plot_var, int):
            var = grid.q[plot_var,:,:]
        else:
            try:
                var = plot_var(current_data)
            except:
                print '*** Error applying function plot_var = ',plot_var
                traceback.print_exc()
                return 

        thisgridvar = Data()
        thisgridvar.var = var
        thisgridvar.xc_center = xc_center
        thisgridvar.yc_center = yc_center
        thisgridvar.xc_edge = xc_edge
        thisgridvar.yc_edge = yc_edge
        thisgridvar.mx = grid.dimensions[0].n
        thisgridvar.my = grid.dimensions[1].n

        # end of 2d case

    return thisgridvar


#------------------------------------------------------------------------
def printfig(fname='',frameno='', figno='', format='png', plotdir='.', \
             verbose=True):
#------------------------------------------------------------------------
    """
    Save the current plot to file fname or standard name from frame/fig.
.  
    If fname is nonempty it is used as the filename, with extension
    determined by format if it does not already have a valid extension.

    If fname=='' then save to file frame000NfigJ.ext  where N is the frame
    number frameno passed in, J is the figure number figno passed in,
    and the extension ext is determined by format.  
    If figno='' then the figJ part is omitted.
    """

    if fname == '':
        fname = 'frame' + str(frameno).rjust(4,'0') 
        if isinstance(figno,int):
            fname = fname + 'fig%s' % figno
    splitfname = os.path.splitext(fname)
    if splitfname[1] not in ('.png','.emf','.eps','.pdf'):
        fname = splitfname[0] + '.%s' % format
    if figno=='':
        figno = 1
    pylab.figure(figno)
    if plotdir != '.':
       fname = os.path.join(plotdir,fname)
    if verbose:  print '    Saving plot to file ', fname
    pylab.savefig(fname)



#======================================================================
def printframes(plotdata=None, verbose=True):
#======================================================================

    """
    Deprecated: use plotpages.plotclaw_driver instead to get gauges as well.
      - RJL, 1/1/10

    Produce a set of png files for all the figures specified by plotdata.
    Also produce a set of html files for viewing the figures and navigating
    between them.  These will all be in directorey plotdata.plotdir.

    The ClawPlotData object plotdata will be initialized by a call to
    function setplot unless plotdata.setplot=False.  

    If plotdata.setplot=True then it is assumed that the current directory
    contains a module setplot.py that defines this function.

    If plotdata.setplot is a string then it is assumed this is the name of
    a module to import that contains the function setplot.

    If plotdata.setplot is a function then this function will be used.
    """

    import glob
    from pyclaw.plotters.data import ClawPlotData



    if not sys.modules.has_key('matplotlib'):
        print '*** Error: matplotlib not found, no plots will be done'
        return plotdata
        
    if not isinstance(plotdata,ClawPlotData):
        print '*** Error, plotdata must be an object of type ClawPlotData'
        return plotdata

    plotdata._mode = 'printframes'

    plotdata = call_setplot(plotdata.setplot, plotdata)

    try:
        plotdata.rundir = os.path.abspath(plotdata.rundir)
        plotdata.outdir = os.path.abspath(plotdata.outdir)
        plotdata.plotdir = os.path.abspath(plotdata.plotdir)

        framenos = plotdata.print_framenos # frames to plot
        fignos = plotdata.print_fignos     # figures to plot at each frame
        fignames = {}                      # names to use in html files

        rundir = plotdata.rundir       # directory containing *.data files
        outdir = plotdata.outdir       # directory containing fort.* files
        plotdir = plotdata.plotdir     # where to put png and html files
        overwrite = plotdata.overwrite # ok to overwrite?
        msgfile = plotdata.msgfile     # where to write error messages
    except:
        print '*** Error in printframes: plotdata missing attribute'
        print '  *** plotdata = ',plotdata
        return plotdata

    if fignos == 'all':
        fignos = plotdata._fignos
        #for (figname,plotfigure) in plotdata.plotfigure_dict.iteritems():
        #    fignos.append(plotfigure.figno)


    # filter out the fignos that will be empty, i.e.  plotfigure._show=False.
    plotdata = set_show(plotdata)
    fignos_to_show = []
    for figname in plotdata._fignames:
        figno = plotdata.plotfigure_dict[figname].figno
        if (figno in fignos) and plotdata.plotfigure_dict[figname]._show:
            fignos_to_show.append(figno)
    fignos = fignos_to_show
        

    rootdir = os.getcwd()

    # annoying fix needed when EPD is used for plotting under cygwin:
    if rootdir[0:9] == 'C:\cygwin' and outdir[0:9] != 'C:\cygwin':
        outdir = 'C:\cygwin' + outdir
        plotdata.outdir = outdir
    if rootdir[0:9] == 'C:\cygwin' and rundir[0:9] != 'C:\cygwin':
        rundir = 'C:\cygwin' + rundir
        plotdata.rundir = rundir
    if rootdir[0:9] == 'C:\cygwin' and plotdir[0:9] != 'C:\cygwin':
        plotdir = 'C:\cygwin' + plotdir
        plotdata.plotdir = plotdir

    try:
        os.chdir(rundir)
    except:
        print '*** Error: cannot move to run directory ',rundir
        print 'rootdir = ',rootdir
        return plotdata


    if msgfile != '':
        sys.stdout = open(msgfile, 'w')
        sys.stderr = sys.stdout


    try:
        plotpages.cd_plotdir(plotdata)
    except:
        print "*** Error, aborting plotframes"
        return plotdata


    framefiles = glob.glob(os.path.join(plotdir,'frame*.png')) + \
                    glob.glob(os.path.join(plotdir,'frame*.html'))
    if overwrite:
        # remove any old versions:
        for file in framefiles:
            os.remove(file)
    else:
        if len(framefiles) > 1:
            print "*** Remove frame*.png and frame*.html and try again,"
            print "  or use overwrite=True in call to printframes"
            return plotdata

    
    # Create each of the figures
    #---------------------------

    try:
        os.chdir(outdir)
    except:
        print '*** Error printframes: cannot move to outdir = ',outdir
        return plotdata


    fortfile = {}
    pngfile = {}
    frametimes = {}

    import glob
    for file in glob.glob('fort.q*'):
        frameno = int(file[7:10])
        fortfile[frameno] = file
        for figno in fignos:
            pngfile[frameno,figno] = 'frame' + file[-4:] + 'fig%s.png' % figno
    
    if len(fortfile) == 0:
        print '*** No fort.q files found in directory ', os.getcwd()
        return plotdata
    
    # Discard frames that are not from latest run, based on
    # file modification time:
    framenos = only_most_recent(framenos, plotdata.outdir)

    numframes = len(framenos)

    print "Will plot %i frames numbered:" % numframes, framenos
    print 'Will make %i figure(s) for each frame, numbered: ' % len(fignos),\
          fignos

    fignames = {}
    for figname in plotdata._fignames:
        figno = plotdata.plotfigure_dict[figname].figno
        fignames[figno] = figname


    
    # Make png files by calling plotframe:
    # ------------------------------------

    for frameno in framenos:
        #plotframe(frameno, plotdata, verbose)
        frametimes[frameno] = plotdata.getframe(frameno, plotdata.outdir).t
        #print 'Frame %i at time t = %s' % (frameno, frametimes[frameno])

    plotdata.timeframes_framenos = framenos
    plotdata.timeframes_frametimes = frametimes
    plotdata.timeframes_fignos = fignos
    plotdata.timeframes_fignames = fignames

    os.chdir(plotdir)

    if plotdata.html:
        plotpages.timeframes2html(plotdata)

    if not plotdata.printfigs:
        print "Using previously printed figure files"
    else:
        print "Now making png files for all figures..."
        for frameno in framenos:
            plotframe(frameno, plotdata, verbose)
            #frametimes[frameno] = plotdata.framesoln_dict[frameno].t
            print 'Frame %i at time t = %s' % (frameno, frametimes[frameno])

    if plotdata.latex:
        plotpages.timeframes2latex(plotdata)
    

    # Movie:
    #-------
    
    if plotdata.gif_movie:
        print 'Making gif movies.  This may take some time....'
        for figno in fignos:
            try:
                os.system('convert -delay 20 frame*fig%s.png moviefig%s.gif' \
                   % (figno,figno))
                print '    Created moviefig%s.gif' % figno
            except:
                print '*** Error creating moviefig%s.gif' % figno
    
    os.chdir(rootdir)

    # print out pointers to html index page:
    path_to_html_index = os.path.join(os.path.abspath(plotdata.plotdir), \
                               plotdata.html_index_fname)
    plotpages.print_html_pointers(path_to_html_index)

    # reset stdout for future print statements
    sys.stdout = sys.__stdout__

    return plotdata
    # end of printframes

#======================================================================
def plotclaw2html(plotdata=None, verbose=True):
#======================================================================
    """
    Old plotting routine no longer exists.
    """
    print '*** Error: plotclaw2html name is deprecated, use printframes.'
    print '*** Plotting command syntax has also changed, you may need to'
    print '*** update your setplot function.'



#------------------------------------------------------------------
def var_limits(plotdata,vars,padding=0.1):
#------------------------------------------------------------------
    """
    Determine range of values encountered in data for all frames
    in the most recent computation in order to determine appropriate
    y-axis limits for plots.

    vars is a list of variables to return the limits for.
    Each element of the list can be a number, indicating the 
    corresponding element of the q vector, or a function that should
    be applied to q to obtain the value of interest.  In other words,
    any valid plot_var attribute of a ClawPlotItem object.

    If vars=='all', use vars=[0,1,...] over all components of q.
    *** not implemented yet ***

    Returns: varmin, varmax, varlim,  
    each is a dictionary with keys consisting by the elements of vars.

    For each var in vars, 
      varmin[var] is the minimum of var over all frames,
      varmax[var] is the maximum of var over all frames,
      varlim[var] is of the form [v1, v2] suitable to use as the 
            ylimits attritube of an object of class ClawPlotAxes with
                v1 = vmin - padding*(vmax-vmin)
                v2 = vmax + padding*(vmax-vmin)
            to give some space above and below the min and max.

    """

    varlim = {}
    vmin,vmax = var_minmax(plotdata,'all',vars)
    varmin = {}
    varmax = {}
    for var in vars:
        varmin[var] = vmin[var]['all']   # min over all frames
        varmax[var] = vmax[var]['all']   # max over all frames
        v1 = varmin[var] - padding*(varmax[var]-varmin[var])
        v2 = varmax[var] + padding*(varmax[var]-varmin[var])
        varlim[var] = [v1,v2]

    return varmin,varmax,varlim
        

#------------------------------------------------------------------
def var_minmax(plotdata,framenos,vars):
#------------------------------------------------------------------
    """
    Determine range of values encountered in data for all frames
    in the list framenos.  If framenos='all', search over all frames
    in the most recent computation.

    vars is a list of variables to return the min and max for.
    Each element of the list can be a number, indicating the 
    corresponding element of the q vector, or a function that should
    be applied to q to obtain the value of interest.  In other words,
    any valid plot_var attribute of a ClawPlotItem object.

    If vars=='all', use vars=[0,1,...] over all components of q.
    *** not implemented yet ***

    Returns (varmin, varmax) where varmin and varmax are each 
    dictionaries with keys given by the variables.  

    For each var in vars, varmin[var] and varmax[var] are dictionaries, with
    keys given by frame numbers in framenos.
    
    Also varmin[var]['all'] and varmax[var]['all'] are the min and max of
    variable var over all frames considered.

    For example if vars = [0,'machnumber'] then
       varmin[0][3] is the minimum of q[0] (the first component of the
           solution vector q) over all grids in frame 3.
       varmin['machnumber']['all'] is the minimum of machnumber 
           over all grids in all frames.

    """

    from pylab import inf
    framenos = only_most_recent(framenos, plotdata.outdir)
    if len(framenos) == 0:
        print '*** No frames found in var_minmax!'

    print 'Determining min and max of variables: ',vars
    print '   over frames: ',framenos

    varmin = {}
    varmax = {}
    if vars=='all':
        print "*** Error in var_minmax: vars == 'all' is not implemented yet"
        return (varmin,varmax)
    for var in vars:
        varmin[var] = {'all': inf}
        varmax[var] = {'all': -inf}
        for frameno in framenos:
            varmin[var][frameno] = inf
            varmax[var][frameno] = -inf

    for frameno in framenos:
        solution = plotdata.getframe(frameno, plotdata.outdir)
        ndim = solution.ndim

        for igrid in range(len(solution.grids)):
            grid = solution.grids[igrid]
            for var in vars:
                if isinstance(var,int):
                    if ndim == 1:
                        qvar = grid.q[var,:]
                    elif ndim == 2:
                        qvar = grid.q[var,:,:]
                    elif ndim == 3:
                        qvar = grid.q[var,:,:,:]
                else:
                    t = solution.t
                    #grid.compute_physical_coordinates()
                    if ndim == 1:
                        X_center = grid.p_center[0]
                        qvar = var(grid.q, X_center, t)
                    elif ndim == 2:
                        X_center, Y_center = grid.p_center
                        qvar = var(grid.q, X_center, \
                                   Y_center, t)
                    elif ndim == 3:
                        X_center, Y_center, Z_center = grid.p_center
                        qvar = var(grid.q, X_center, \
                                   Y_center, Z_center, t)
                varmin[var][frameno] = min(varmin[var][frameno], qvar.min())
                varmax[var][frameno] = max(varmax[var][frameno], qvar.max())
                varmin[var]['all'] = min(varmin[var]['all'], \
                                         varmin[var][frameno])
                varmax[var]['all'] = max(varmax[var]['all'], \
                                         varmax[var][frameno])
    return (varmin, varmax)
        


#------------------------------------------------------------------
def only_most_recent(framenos,outdir='.',verbose=True):
#------------------------------------------------------------------

    """
    Filter the list framenos of frame numbers and return a new
    list that only contains the ones from the most recent run.
    This is determined by comparing modification times of 
    fort.q files.

    Returns the filtered list.
    """

    import glob,time,os

    startdir = os.getcwd()
    if outdir != '.':
        try:
            os.chdir(outdir)
        except:
            print "*** Could not chdir to ", outdir
            return framenos

    fortfile = {}
    for file in glob.glob('fort.q*'):
        frameno = int(file[7:10])
        fortfile[frameno] = file

    #DK: In PetClaw, we don't output fort.q* files.  Instead count the
    #claw.pkl* files.
    if len(fortfile) == 0:
        for file in glob.glob('claw.pkl*'):
            frameno = int(file[9:12])
            fortfile[frameno] = file
    
    if len(fortfile) == 0:
        print '*** No fort.q or claw.pkl files found in directory ', os.getcwd()
        framenos = []
        return framenos
    
    # Figure out which files are from latest run:
    numframes = 0
    mtime = 0
    framekeys = fortfile.keys()
    framekeys.sort()
    for frameno in framekeys:
        mtimeprev = mtime
        mtime = os.path.getmtime(fortfile[frameno])
        # sometimes later fort files are closed a few seconds after
        # earlier ones, so include a possible delaytime:
        delaytime = 5  # seconds
        if mtime < mtimeprev-delaytime:
            break
        numframes = numframes + 1
    
    newframes = framekeys[:numframes]
    if (numframes < len(framekeys)) & verbose:
        print '*** Frames %s and above appear to be from an old run' \
                       % framekeys[numframes]
        print '***    and will be ignored.'
        time.sleep(2)

    #print 'framenos = ',framenos
    if framenos == 'all':
        framenos = newframes
    else:
        # compute intersection of framenos and newframes:
        framenos = list(set(framenos).intersection(set(newframes)))
    framenos.sort()
    os.chdir(startdir)
    return framenos

#------------------------------------------------------------------------
def call_setplot(setplot, plotdata, verbose=True):
#------------------------------------------------------------------
    """
    Try to apply setplot to plotdata and return the result.
    If setplot is False or None, return plotdata unchanged.
    If setplot is True, setplot function is in setplot.py.
    If setplot is a string, setplot function is in module named by string.
    Otherwise assume setplot is a function.
    """

    if not setplot:
        if verbose:
            print '*** Warning: no setplot specified'
        return plotdata

    if setplot is True:
        # indicates we should import setplot.py from current directory
        # and the setplot function is in this module.  
        try:
            sys.path.insert(0,os.getcwd())
            import setplot as SetPlot
            # make sure setplot from current directory is used:
            reload(SetPlot)
            setplot = SetPlot.setplot   # should be a function

            if verbose:
                print 'Imported setplot from ', os.getcwd()
        except:
            print """*** Error in call_setplot: 
                  Problem importing setplot.py in directory %s""" \
                  % os.getcwd()
            raise
            return plotdata

    elif isinstance(setplot,str):
        # assume setplot specifies module containing setplot function
        # strip off the .py if it is there:
        setplotmod = os.path.splitext(setplot)[0]
        try:
            sys.path.insert(0,os.getcwd())
            exec('import %s as SetPlot' % setplotmod)
            # make sure setplot from current directory is used:
            reload(SetPlot)
            setplot = SetPlot.setplot   # should be a function
        except:
            print """*** Error in call_setplot: 
                  Problem importing from %s.py in directory %s""" \
                  % (setplotmod, os.getcwd())
            raise
            return plotdata
    else:
        # assume setplot is the setplot function itself:
        pass

    # execute setplot:
    try:
        plotdata = setplot(plotdata)
        if verbose:
            print 'Executed setplot successfully'
    except: #Exception as error:
        print '*** Error in call_setplot: Problem executing function setplot'
        raise
        return plotdata

    if plotdata is None:
        print '*** Error in printframes:  plotdata is None after call to setplot'
        print '*** Did you forget the "return plotdata" statement?'
        raise
        return plotdata
    
    return plotdata

#------------------------------------------------------------------
def clawpack_header():
#------------------------------------------------------------------
    pylab.axes([.3, .8, .98, 1.])
    pylab.text(.1,.13,'Clawpack Plots',fontsize=30,color='brown')
    pylab.axis('off')

#------------------------------------------------------------------
def errors_2d_vs_1d(solution,reference,var_2d,var_1d,map_2d_to_1d):
#------------------------------------------------------------------
    """
    Input: 
      solution: object of class ClawSolution with 2d computed solution

      reference: object of class ClawSolution with 1d reference solution

      var_2d: variable to compare from solution
      if integer, compare q[var_2d,:,:]
      if function, apply to q to obtain the variable

      var_1d: variable to compare from reference
      if integer, compare to q[var_2d,:]
      if function, apply to q to obtain the variable

      map_2d_to_1d: function mapping 2d data to 1d for comparison::

          xs, qs = map_2d_to_1d(xcenter, ycenter, q)
    """

    from numpy import interp

    t = solution.t
    xs = {}
    qs = {}
    qint = {}

    errmax = 0.
    for gridno in range(len(solution.grids)):
        grid = solution.grids[gridno]

        #grid.compute_physical_coordinates()
        X_center, Y_center = grid.p_center
        X_edge, Y_edge = grid.p_center
    
        if isinstance(var_2d, int):
            q = grid.q[var_2d,:,:]
        else:
            try:
                q = var_2d(grid.q, X_center, Y_center, t)
            except:
                print '*** Error applying function plot_var = ',plot_var
                traceback.print_exc()
                return 
        
        xs1, qs1 = map_2d_to_1d(q, X_center, Y_center, t)

        if hasattr(reference,'grids'):
            if len(reference.grids) > 1:
                print '*** Warning in errors_2d_vs_1d: reference solution'
                print '*** has more than one grid -- only using grid[0]'
            refgrid = reference.grids[0]
        else:
            refgrid = reference   # assume this contains true solution or
                                  # something set separately rather than
                                  # a framesoln

        #refgrid.compute_physical_coordinates()
        xref = grid.p_center[0]
        if isinstance(var_1d, int):
            qref = refgrid.q[var_1d,:].T
        else:
            try:
                qref = var_1d(refgrid.q, xref, t)
            except:
                print '*** Error applying function var_1d'
                return 
            
        qint1 = interp(xs1, xref, qref)

        xs[gridno] = xs1
        qs[gridno] = qs1
        qint[gridno] = qint1
        errabs = abs(qs1 - qint1)
        errmax = max(errmax, errabs.max())
        #import pdb
        #pdb.set_trace()

    return errmax, xs, qs, qint
        


#------------------------------------------------------------------
def set_show(plotdata):
#------------------------------------------------------------------
    """
    Determine which figures and axes should be shown.
    plotaxes._show should be true only if plotaxes.show and at least one
    item is showing or if the axes have attribute type=='empty', in which
    case something may be plotting in an afteraxes command, for example.

    plotfigure._show should be true only if plotfigure.show and at least
    one axes is showing.
    """

    for figname in plotdata._fignames:
        plotfigure = plotdata.plotfigure_dict[figname]
        plotfigure._show = False
        if plotfigure.show:
            # Loop through all axes to make sure at least some item is showing
            for plotaxes in plotfigure.plotaxes_dict.itervalues():
                plotaxes._show = False
                if plotaxes.show:
                    # Loop through plotitems checking each item to see if it 
                    # should be shown
                    for plotitem in plotaxes.plotitem_dict.itervalues():
                        plotitem._show = plotitem.show
                        if plotitem.show:
                            plotaxes._show = True
                            plotfigure._show = True
                    # Check to see if the axes are supposed to be empty or 
                    # something may be in the afteraxes function
                    if not plotaxes._show:
                        if plotaxes.afteraxes is not None or plotaxes.type == 'empty':
                            plotaxes._show = True
                            plotfigure._show = True
                
    return plotdata

#------------------------------------------------------------------
def set_outdirs(plotdata):
#------------------------------------------------------------------
    """
    Make a list of all outdir's for all plotitem's in the order they
    are first used.
    """

    outdir_list = []
    for figname in plotdata._fignames:
        plotfigure = plotdata.plotfigure_dict[figname]
        if not plotfigure._show:
            continue  # skip to next figure
        for axesname in plotfigure._axesnames:
            plotaxes = plotfigure.plotaxes_dict[axesname]
            if not plotaxes._show:
                continue  # skip to next axes
            for itemname in plotaxes._itemnames:
                plotitem = plotaxes.plotitem_dict[itemname]
                if not plotitem._show:
                    continue  # skip to next item
                if plotitem.outdir is not None:
                    outdir = plotitem.outdir
                else:
                    outdir = plotdata.outdir
                if outdir not in outdir_list:
                    outdir_list.append(outdir)
                
    plotdata._outdirs = outdir_list
    return plotdata
