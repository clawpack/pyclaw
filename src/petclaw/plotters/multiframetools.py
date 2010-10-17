
# This is not yet ready for use!
    

#==========================================
def plot_multiframes(plotdata, verbose=False):
#==========================================

    """
    Plot all figures that involve multiple frames

    """

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


    # initialize current_data containing data that will be passed
    # to afterframe, afteraxes, aftergrid commands
    current_data = Data()
    current_data.user = Data()   # for user specified attributes
                                 # to avoid potential conflicts
    current_data.plotdata = plotdata


 
    if plotdata._mode == 'iplotclaw':
        print '    Plotting Frame %s at t = %s' % (frameno,t)
        requested_fignos = plotdata.iplotclaw_fignos
    else:
        requested_fignos = plotdata.print_fignos
    plotted_fignos = []

    plotdata = set_show(plotdata)   # set _show attributes for which figures
                                    # and axes should be shown.

    # loop over figures to look for multiframe plots
    # -------------------------------------------

    for figname in plotdata._fignames:
        plotfigure = plotdata.plotfigure_dict[figname]
        if not plotfigure._show:
            continue  # skip to next figure 

        figno = plotfigure.figno
        if requested_fignos != 'all':
            if figno not in requested_fignos:
                continue # skip to next figure

        plotted_fignos.append(figno)


        try:
            plotaxes_dict = plotfigure.plotaxes_dict
        except:
            print '*** Error in plotframe: plotdata missing plotaxes_dict'
            print '*** This should not happen'
            return  None

        if (len(plotaxes_dict) == 0) or (len(plotfigure._axesnames) == 0):
            print '*** Warning in plotframe: plotdata has empty plotaxes_dict'
            print '*** Apparently no axes to plot in figno ',figno

        # loop over axes looking for type multiframe 
        # ------------------------------------------

        thisfigure_multiframe = False
        for axesname in plotfigure._axesnames:
            plotaxes = plotaxes_dict[axesname]
            if (plotaxes._show) and (plotaxes.type == 'multiframe'):
                thisfigure_multiframe = True
        if not thisfigure_multiframe:
            continue
        
        if not plotfigure.kwargs.has_key('facecolor'):
            # use Clawpack's default bg color (tan)
            plotfigure.kwargs['facecolor'] = '#ffeebb'   

        # create figure and set handle:
        plotfigure._handle = pylab.figure(num=figno, **plotfigure.kwargs)

        pylab.ioff()
        if plotfigure.clf_each_frame:
            pylab.clf()


        # loop over axes to appear on this figure:
        # ----------------------------------------

        for axesname in plotfigure._axesnames:
            plotaxes = plotaxes_dict[axesname]
            if (not plotaxes._show) or (plotaxes.type != 'multiframe'):
                continue   # skip this axes if no items show

            # create the axes:
            axescmd = getattr(plotaxes,'axescmd','subplot(1,1,1)')
            axescmd = 'plotaxes._handle = pylab.%s' % axescmd
            exec(axescmd)
            pylab.hold(True)

            # loop over plotitems on these axes:
            # ----------------------------------

            for itemname in plotaxes._itemnames:
                
                plotitem = plotaxes.plotitem_dict[itemname]
                if plotitem._show == False:
                    print 'Skipping item ', itemname
                    continue
                if plotitem.plot_type not in ['1d_vs_t', '1d_xt']:
                    continue

                # Get solution for required frames:
                framenos = only_most_recent(plotitem.framenos)
                if len(framenos)==0:
                    print "*** No frames found in plot_multiframes for "
                    print "*** plotitem named ",itemname
                    continue
                else if plotitem.plot_type == '1d_xt':
                    plot_var = plotitem.plot_var
                    var_over_t = []
                    times = []
                    for frameno in framenos:
                        framesoln = plotdata.getframe(frameno, plotdata.outdir)
                        if len(framesoln.grids) > 1:
                            print "*** plot_type 1d_xt can't be used with AMR"
                            return None
                        grid = framesoln.grids[0]
                        thisgridvar = get_gridvar(grid,plot_var,1,current_data)
                        var_over_t.append(thisgridvar)
                        times.append(framesoln.t)
                    x = thisgridvar.xc_center   # cell centers
                    var = array(var_over_t)
                    times = array(times)
                    contour(x,times,var)



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

            if plotaxes.scaled:
                pylab.axis('scaled')


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

