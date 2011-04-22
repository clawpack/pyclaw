"""
Tools for plotting data from gauges, gauge locations, etc.
"""

import os,sys,shutil,glob
import string,re
import time
import traceback


from pyclaw.data import Data
from pyclaw.plotters import plotpages
from pyclaw.plotters.frametools import set_show

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



#==========================================
def plotgauge(gaugeno, plotdata, verbose=False):
#==========================================

    """
    Plot all requested plots for a single gauge from the computation.
    The plots are requested by setting attributes of plotdata
    to ClawPlotFigure objects with plot_type="each_gauge".

    """


    if verbose:  
        gaugesoln = plotdata.getgauge(gaugeno)
        print '    Plotting gauge %s  at x = %g, y = %g ... '  \
                 % (gaugeno, gaugesoln.x, gaugesoln.y)

    if plotdata.mode() == 'iplotclaw':
        pylab.ion()

        
    try:
        plotfigure_dict = plotdata.plotfigure_dict
    except:
        print '*** Error in plotgauge: plotdata missing plotfigure_dict'
        print '*** This should not happen'
        return None

    if len(plotfigure_dict) == 0:
        print '*** Warning in plotgauge: plotdata has empty plotfigure_dict'
        print '*** Apparently no figures to plot'




    # initialize current_data containing data that will be passed
    # to beforegauge, aftergauge, afteraxes commands
    current_data = Data()
    current_data.user = Data()   # for user specified attributes
                                 # to avoid potential conflicts
    current_data.plotdata = plotdata
    current_data.gaugeno = gaugeno

    # call beforegauge if present, which might define additional 
    # attributes in current_data or otherwise set up plotting for this
    # gauge.

    beforegauge =  getattr(plotdata, 'beforegauge', None)
    if beforegauge:
        if isinstance(beforegauge, str):
            # a string to be executed
            exec(beforegauge)
        else:
            # assume it's a function
            try:
                output = beforegauge(current_data)
                if output: current_data = output
            except:
                print '*** Error in beforegauge ***'
                raise



    # iterate over each single plot that makes up this gauge:
    # -------------------------------------------------------
 
    if plotdata._mode == 'iplotclaw':
        gaugesoln = plotdata.getgauge(gaugeno)
        #import pdb; pdb.set_trace()
        print '    Plotting Gauge %s  at x = %g, y = %g ... '  \
                 % (gaugeno, gaugesoln.x, gaugesoln.y)
        requested_fignos = plotdata.iplotclaw_fignos
    else:
        requested_fignos = plotdata.print_fignos
    plotted_fignos = []

    plotdata = set_show(plotdata)   # set _show attributes for which figures
                                    # and axes should be shown.

    # loop over figures to appear for this gauge: 
    # -------------------------------------------

    for figname in plotdata._fignames:
        plotfigure = plotdata.plotfigure_dict[figname]
        if (not plotfigure._show) or (plotfigure.type != 'each_gauge'):
            continue  # skip to next figure 

        figno = plotfigure.figno
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
        if plotfigure.clf_each_gauge:
            pylab.clf()

        try:
            plotaxes_dict = plotfigure.plotaxes_dict
        except:
            print '*** Error in plotgauge: plotdata missing plotaxes_dict'
            print '*** This should not happen'
            return  None

        if (len(plotaxes_dict) == 0) or (len(plotfigure._axesnames) == 0):
            print '*** Warning in plotgauge: plotdata has empty plotaxes_dict'
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



            # loop over items:
            # ----------------

            for itemname in plotaxes._itemnames:
                
                plotitem = plotaxes.plotitem_dict[itemname]
                try:
                    outdir = plotitem.outdir
                    if outdir is None:
                        outdir = plotdata.outdir
                    gaugesoln = plotdata.getgauge(gaugeno, outdir)
                except:
                    print '*** Cannot find gauge number ',gaugeno
                    print '*** looking in directory ', outdir
                    print '*** cwd = ',os.getcwd()
                    return None

                #import pdb; pdb.set_trace()
                current_data.gaugesoln = gaugesoln
                current_data.q = gaugesoln.q
                current_data.t = gaugesoln.t

                if plotitem._show:
                    try:
                        output = plotgauge1(gaugesoln,plotitem,\
                            current_data)
                        if output: current_data = output
                        if verbose:  
                                print '      Plotted  plotitem ', itemname
                    except:
                        print '*** Error in plotgauge: problem calling plotgauge1'
                        traceback.print_exc()
                        return None

            # end of loop over plotitems


        for itemname in plotaxes._itemnames:
            plotitem = plotaxes.plotitem_dict[itemname]
            if plotitem.afteritem:
                print "*** ClawPlotItem.afteritem is deprecated"
                print "*** use ClawPlotAxes.afteraxes "
                print "*** or  ClawPlotItem.aftergrid instead"


        pylab.title("%s at gauge %s" % (plotaxes.title,gaugeno))


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


    # call an aftergauge function if present:
    aftergauge =  getattr(plotdata, 'aftergauge', None)
    if aftergauge:
        if isinstance(aftergauge, str):
            # a string to be executed
            exec(aftergauge)
        else:
            # assume it's a function
            try:
                output = aftergauge(current_data)
                if output: current_data = output
            except:
                print '*** Error in aftergauge ***'
                raise


    if plotdata.mode() == 'iplotclaw':
        pylab.ion()
    for figno in plotted_fignos:
        pylab.figure(figno)
        pylab.draw()

    if verbose:
        print '    Done with plotgauge for gauge %i' % (gaugeno)

    
    # print the figure(s) to file(s) if requested:
    if (plotdata.mode() != 'iplotclaw') & plotdata.printfigs:
        # iterate over all figures that are to be printed:
        for figno in plotted_fignos:
            printfig(gaugeno=gaugeno, figno=figno, \
                    format=plotdata.print_format, plotdir=plotdata.plotdir,\
                    verbose=verbose)

    return current_data

    # end of plotgauge

    
#==================================================================
def plotgauge1(gaugesoln, plotitem, current_data):
#==================================================================
    """
    Make a 1d plot for a single plot item for the gauge solution in
    gaugesoln.

    The current_data object holds data that should be passed into
    aftergrid or afteraxes if these functions are defined.  The functions
    may add to this object, so this function should return the possibly
    modified current_data for use in other plotitems or in afteraxes or
    afterframe.

    """

    plotdata = plotitem._plotdata
    plotfigure = plotitem._plotfigure
    plotaxes = plotitem._plotaxes


    # the following plot parameters may be set, depending on what
    # plot_type was requested:

    plot_params = """
             plot_var  aftergrid  plotstyle color kwargs 
             plot_var2 fill_where map_2d_to_1d 
             """.split()

    # No amr_ parameters for gauge data.

    plot_var = plotitem.plot_var
    plot_type = plotitem.plot_type
    kwargs = plotitem.kwargs
    color = plotitem.color
    plotstyle = plotitem.plotstyle

    t = gaugesoln.t
    if type(plot_var) is int:
        #import pdb pdb.set_trace()
        var = gaugesoln.q[:,plot_var]
    else:
        try:
            var = plot_var(gaugesoln)
        except:
            raise Exception("Problem applying plot_var to gaugesoln")
    tmax = t.max()
    varmax = var.max()

    # The plot commands using matplotlib:

    pylab.hold(True)

    pylab.title("%s at Gauge %i" % (plotitem._plotaxes.title,\
                 gaugesoln.gaugeno))

    pylab.xlabel("time")

    if (plot_type in ['1d_plot']) and (plotstyle != ''):
        if color:
            kwargs['color'] = color

        plotcommand = "pobj=pylab.plot(t,var,'%s', **kwargs)"  \
                      % plotstyle
        exec(plotcommand)


    elif plot_type == '1d_empty':
        # no plot to create (user might make one in afteritem or
        # afteraxes)
        pass

    else:
        raise ValueError("Unrecognized plot_type: %s" % plot_type)
        return None

    return current_data
    
def read_setgauges(datadir):
    """
    Read the info from setgauges.data.
    """
    import os
    import numpy as np
    from pyclaw.data import Data
    from matplotlib.mlab import find

    setgauges = Data()

    # default values if no gauges found:
    setgauges.numgauges = 0
    setgauges.gaugenos = []
    setgauges.x = {}
    setgauges.y = {}
    setgauges.t1 = {}
    setgauges.t2 = {}

    fname = os.path.join(datadir, 'setgauges.data')
    if not os.path.isfile(fname):
        #print "*** Warning in read_setgauges: missing file ",fname
        return setgauges

    sgfile = open(fname,'r')
    lines = sgfile.readlines()
    sgfile.close()

    lineno = 0
    ignore_line = True
    while ignore_line:
        line = lines[lineno] + "#"
        #print "+++ lineno = %s, line = %s" % (lineno,line)
        if line.split()[0][0]=="#":
            lineno = lineno+1
        else:
            ignore_line = False

    try:
        numgauges = int(line.split()[0])
    except:
        print "*** error setting numgauges"
        return

    if numgauges==0:
        return setgauges
        
    #print '+++ ignoring %s lines, numgauges = %s' %(lineno, numgauges)
    try:
        sgno, x, y, t1, t2 = np.loadtxt(fname, unpack=True, skiprows=lineno+1, \
                                usecols=range(5))
        if numgauges==1:
            # loadtxt returns numbers rather than arrays in this case:
            sgno = [sgno]; x = [x]; y = [y]; t1 = [t1]; t2 = [t2]

    except:
        print "*** problem reading gauges from setgauges.data"
        return setgauges

    sgno = np.array(sgno, dtype=int)  # convert to int
    setgauges.gaugenos = sgno
    setgauges.numgauges = numgauges

    for n in sgno:
        nn = find(sgno==n)
        if len(nn) > 1:
            print "*** Warning: found more than one gauge numbered ",n

        if len(nn) == 0:
            print "*** Error: didn't find gauge number %s in %s" % (n,fname)
        else:
            nn = nn[0]
            setgauges.x[n] = x[nn]
            setgauges.y[n] = y[nn]
            setgauges.t1[n] = t1[nn]
            setgauges.t2[n] = t2[nn]

    return setgauges


def plot_gauge_locations(plotdata, gaugenos='all', \
                format_string='ko',add_labels=True, \
                markersize=5, fontsize=15, xoffset=0, yoffset=0):
    """
    Plot gauge locations on current axes.
    format_string specifies the symbol to be plotted.
    If add_labels==True then labels will also be added (gauge number).
    This routine determines locations from the file setgauges.data in
    directory plotdata.rundir.  It does not require reading in the fort.gauge file
    produced by running the code.
    """

    from pylab import figure, plot, clf, title, text

    datadir = plotdata.rundir  # this should contain setgauges.data

    try:
        setgauges = read_setgauges(datadir)
    except:
        return

    if setgauges.numgauges == 0:
        print "*** plot_gauge_locations: No gauges specified in setgauges.data"
        return

    if gaugenos=='all':
        gaugenos = setgauges.gaugenos


    for n in gaugenos:
        try:
            xn = setgauges.x[n]
            yn = setgauges.y[n]
            #print "Gauge %s:  x = %g, y = %g" % (n,xn,yn)
            plot([xn], [yn], format_string, markersize=markersize)
            if add_labels: 
                xn = xn + xoffset
                yn = yn + yoffset
                text(xn,yn,'  %s' % n, fontsize=fontsize)
        except:
            print "*** plot_gauge_locations: warning: did not find x,y data for gauge ",n


#------------------------------------------------------------------------
def printfig(fname='',gaugeno='', figno='', format='png', plotdir='.', \
             verbose=True):
#------------------------------------------------------------------------
    """
    Save the current plot to file fname or standard name from gauge/fig.
.  
    If fname is nonempty it is used as the filename, with extension
    determined by format if it does not already have a valid extension.

    If fname=='' then save to file gauge000NfigJ.ext  where N is the gauge
    number gaugeno passed in, J is the figure number figno passed in,
    and the extension ext is determined by format.  
    If figno='' then the figJ part is omitted.
    """

    if fname == '':
        fname = 'gauge' + str(gaugeno).rjust(4,'0') 
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
def printgauges(plotdata=None, verbose=True):
#======================================================================

    """
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
        
    # figure out what type each figure is:
    fignos_each_frame = []
    fignos_each_gauge = []
    fignos_each_run = []
    for figno in fignos:
        figname = plotdata._figname_from_num[figno]
        if plotdata.plotfigure_dict[figname].type == 'each_frame':
            fignos_each_frame.append(figno)
        if plotdata.plotfigure_dict[figname].type == 'each_gauge':
            fignos_each_gauge.append(figno)
        if plotdata.plotfigure_dict[figname].type == 'each_run':
            fignos_each_run.append(figno)
        

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
        for figno in fignos_each_frame:
            pngfile[frameno,figno] = 'frame' + file[-4:] + 'fig%s.png' % figno
    
    if len(fortfile) == 0:
        print '*** No fort.q files found in directory ', os.getcwd()
        return plotdata
    
    # Discard frames that are not from latest run, based on
    # file modification time:
    framenos = only_most_recent(framenos, plotdata.outdir)

    numframes = len(framenos)

    print "Will plot %i frames numbered:" % numframes, framenos
    print 'Will make %i figure(s) for each frame, numbered: ' \
          % len(fignos_each_frame), fignos_each_frame

    #fignames = {}
    #for figname in plotdata._fignames:
        #figno = plotdata.plotfigure_dict[figname].figno
        #fignames[figno] = figname
    # use new attribute:
    fignames = plotdata._figname_from_num


    for frameno in framenos:
        frametimes[frameno] = plotdata.getframe(frameno, plotdata.outdir).t

    plotdata.timeframes_framenos = framenos
    plotdata.timeframes_frametimes = frametimes
    plotdata.timeframes_fignos = fignos_each_frame
    plotdata.timeframes_fignames = fignames

    # Gauges:
    gaugenos = plotdata.print_gaugenos

    # Make html files for time frame figures:
    # ---------------------------------------

    os.chdir(plotdir)

    if plotdata.html:
        plotpages.timeframes2html(plotdata)
    
    # Make png files for all frames and gauges:
    # -----------------------------------------

    if not plotdata.printfigs:
        print "Using previously printed figure files"
    else:
        print "Now making png files for all figures..."
        for frameno in framenos:
            plotframe(frameno, plotdata, verbose)
            print 'Frame %i at time t = %s' % (frameno, frametimes[frameno])
        for gaugeno in gaugenos:
            plotgauge(gaugeno, plotdata, verbose)
            print 'Gauge %i ' % gaugeno


    if plotdata.latex:
        plotpages.timeframes2latex(plotdata)
    

    # Movie:
    #-------
    
    if plotdata.gif_movie:
        print 'Making gif movies.  This may take some time....'
        for figno in fignos_each_frame:
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

