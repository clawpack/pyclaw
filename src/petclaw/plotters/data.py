#!/usr/bin/env python
# encoding: utf-8
# ======================================================================
#  Package:     petclaw.plotters
#  File:        data.py
#  Created:     Aug 7, 2009
#  Author:      R.J. LeVeque
# ============================================================================
"""
Plotting Data Module

Contains the general class definition and the subclasses of the Clawpack 
data objects specific to plotting.
"""
import os
import copy
import re
import logging
from petclaw.data import Data


# ============================================================================
#  Subclass ClawPlotData containing data for plotting results
# ============================================================================
class ClawPlotData(Data):
    """ClawPlotData class
    
    Data subclass containing plot data.

    """

    # ========== Initialization routine ======================================
    def __init__(self, data_files=[], controller=None):
        """Initialize a PlotData object
        
        Accepts a list of data_files to be read into and instantiate into one
        ClawPlotData object.  An empty object can be created by not passing 
        anything in for the data_file_list
        """
        plot_attrs = ['rundir','plotdir','outdir','overwrite','plotter',
                           'msgfile','printfig_format','afterframe',
                           'beforeframe','mapc2p',
                           'html_framenos','html_fignos','html_fignames',
                           'ndim','clear_figs', 'setplot', 'eagle',
                           'plotitem_dict', 'html_movies', 'params']

        # Initialize the data object and read the data files
        super(ClawPlotData,self).__init__(data_files,plot_attrs)

        # default values of attributes:

        if controller:
            controller.plotdata = self
            # inherit some values from controller
            self.rundir = copy.copy(controller.rundir)
            self.outdir = copy.copy(controller.outdir)
        else:
            self.rundir = os.getcwd()     # uses *.data from rundir
            self.outdir = os.getcwd()     # where to find fort.* files


        self.plotdir = os.getcwd()      # directory for plots *.png, *.html
        self.overwrite = True           # ok to overwrite old plotdir?
        self.plotter = 'matplotlib'     # backend for plots
        self.msgfile = ''               # where to write error messages
        self.verbose = True             # verbose output?

        self.ion = False                # call ion() or ioff()?

        self.user = Data()              # for user to pass things into
                                        # afterframe, for example
					# Deprecated.

        self.printfigs = True 
        self.print_format = 'png'     
        self.print_framenos = 'all'  # which frames to plot
        self.print_gaugenos = 'all'  # which gauges to plot
        self.print_fignos = 'all'    # which figures to plot each frame

        self.iplotclaw_fignos = 'all'    # which figures to plot interactively

        self.latex = True                # make latex files for figures
        self.latex_fname = 'plots'       # name of latex file
        self.latex_title = 'Clawpack Results'       
        self.latex_framesperpage = 'all' # number of frames on each page
        self.latex_framesperline = 2     # number of frames on each line
        self.latex_figsperline = 'all'   # number of figures on each line
        self.latex_makepdf = False       # run pdflatex on latex file

        self.html = True                # make html files for figures
        self.html_index_fname = '_PlotIndex.html'   # name of html index file
        self.html_index_title = 'Plot Index'   # title at top of index page
        self.html_homelink = None       # link to here from top of _PlotIndex.html
        self.html_movie = True          # make html with java script for movie
        self.html_eagle = False         # use EagleClaw titles on html pages?

        self.gif_movie = False          # make animated gif movie of frames


    #    self.clear_figs = True          # give clf() command in each figure
                                        # before plotting each frame

        self.setplot = False            # Execute setplot.py in plot routine
    #    self.setplot_caller = None      # Set before calling setplot

        self.mapc2p = None              # function to map computational
	                                # points to physical


        self.beforeframe = None         # function called before all plots 
                                        # in each frame are done
        self.afterframe = None          # function called after all plots 
                                        # in each frame are done

        self.plotfigure_dict = {}  

        self.framesoln_dict = {}        # dictionary for holding framesoln
                                        # objects associated with plots

        self.gaugesoln_dict = {}        # dictionary for holding gaugesoln
                                        # objects associated with plots
                                        
        self.save_frames = True         # True ==> Keep a copy of any frame
                                        # read in.  False ==> Clear the frame
                                        # solution dictionary before adding
                                        # another solution

        self.save_figures = True        # True ==> Keep a copy of and figure
                                        # created.  False ==> Clear the 
                                        # figure dictionary before adding
                                        # another solution

        self.refresh_frames = False     # False ==> don't re-read framesoln if 
                                        # already in framesoln_dict

        self.refresh_gauges = False     # False ==> don't re-read gaugesoln if 
                                        # already in gaugesoln_dict


        self._next_FIG = 1000
        self._fignames = []
        self._fignos = []
        self._mode = 'unknown'
        self._figname_from_num = {}

        #if data_file_list is not None:
        if len(data_files) > 0:
            # values in data files may overwrite some default values
            # or set parameter values in params dictionary
            for data_file in data_files:
                self.read(data_file)

    def new_plotfigure(self, name=None, figno=None, type='each_frame'):
        """
        Create a new figure for Clawpack plots.  
        If type='each_frame' it is a figure that will be plotted 
	for each time frame.
        If type='multi_frame' it is a figure that will be plotted based on
	all the frames, such as x-t plots or time series. (Not yet implemented)
        """
        if (self._mode != 'iplotclaw') and (name in self._fignames):
            print '*** Warning, figure named %s has already been created' % name
        if (self._mode != 'iplotclaw') and (figno in self._fignos):
            print '*** Warning, figure number %s has already been created' % figno
        if figno is None:
            self._next_FIG += 1
            figno = self._next_FIG
        if name is None:
            name = "FIG%s" % figno
        if name in self._fignames:
            print "*** Error in new_plotfigure: Figure name already used... ",name
            raise Exception("Figure name already used")
        elif figno in self._fignos:
            print "*** Error in new_plotfigure: Figure number already used... ",figno
            raise Exception("Figure number already used")

        self._fignames.append(name)
        self._fignos.append(figno)
        plotfigure = ClawPlotFigure(name, figno, type, self)
        if not self.save_figures:
            self.plotfigure_dict.clear()
        self.plotfigure_dict[name] = plotfigure
        self._figname_from_num[figno] = name
        return plotfigure


    def getframe(self,frameno,outdir=None):
        """
        ClawPlotData.getframe:
        Return an object of class Solution containing the solution
        for frame number frameno.

        If self.refresh_frames == True then this frame is read from the fort
        files, otherwise it is read from the fort files only if the
        the dictionary self.framesoln_dict has no key frameno.  If it does, the
        frame has previously been read and the dictionary value is returned.
        """

        from petclaw import solution

        framesoln_dict = self.framesoln_dict

        if 0:
            if outdir:
                key = (frameno, outdir)
            else:
                key = frameno
                outdir = self.outdir

        if outdir is None:
            outdir = self.outdir
        outdir = os.path.abspath(outdir)
        key = (frameno, outdir)

        if self.refresh_frames or (not framesoln_dict.has_key(key)):
            thisdir = os.getcwd()
            try:
                os.chdir(outdir)
            except:
                print '*** Error in getframe: cannot move to outdir = ',\
                       outdir
                print '*** thisdir = ',thisdir
                raise
                return
            try:
                framesoln = solution.Solution(frameno)
            except:
                print '*** Error reading frame in ClawPlotData.getframe'
                os.chdir(thisdir)
                raise
                return
            os.chdir(thisdir)
            if not self.save_frames:
                framesoln_dict.clear()
            framesoln_dict[key] = framesoln
            if key != frameno:
                print '    Reading  Frame %s at t = %g  from outdir = %s' \
                    % (frameno,framesoln.t,outdir)
            else:
                print '    Reading  Frame %s at t = %g  ' \
                    % (frameno,framesoln.t)
        else:
            framesoln = self.framesoln_dict[key]

        return framesoln
        
    def gettime(self,frameno,outdir='./'):
        r"""Fetch time from solution corresponding to frame number in outdir
        
        This method only works for ascii formatted files
        """

        from petclaw.io.ascii import read_ascii_t
        t,meqn,ngrids,maux,ndim = read_ascii_t(frameno,path=outdir)
        return t

    def clearfigures(self):
        """
        Clear all plot parameters specifying figures, axes, items.
	Does not clear the frames of solution data already read in.
	  For that use clearframes.
        """

	self.plotfigure_dict.clear()
	self._fignames = []
	self._fignos = []
	self._next_FIG = 1000


    def clearframes(self, framenos='all'):
        """
        Clear one or more frames from self.framesoln_dict.
        Need to add outdir option!
        """

        if isinstance(framenos, int):
            framenos = [framenos]  # turn into a list

        if framenos=='all':
            self.framesoln_dict.clear()
            print 'Cleared all frames'
        else:
            for frameno in framenos:
                xxx = self.plotdata.framesoln_dict.pop(frameno,None)
                if xxx is None:
                   print 'No frame data to clear for frame ',frameno
                else:
                   print 'Cleared data for frame ',frameno


    def getgauge(self, gaugeno, outdir=None):
        """
        ClawPlotData.getgauge:
        Return an object of class GaugeSolution containing the solution
        for gauge number gaugeno.

        If self.refresh_gauges == True then this gauge is read from the
        fort.gauge file, otherwise it is read only if the
        the dictionary self.gaugesoln_dict has no key gaugeno.  If it does, the
        gauge has previously been read and the dictionary value is returned.
        """

        gaugesoln_dict = self.gaugesoln_dict

        if outdir is None:
            outdir = self.outdir
        outdir = os.path.abspath(outdir)
        key = (gaugeno, outdir)

        if self.refresh_gauges or (not gaugesoln_dict.has_key(key)):
            thisdir = os.getcwd()
            try:
                os.chdir(outdir)
            except:
                print '*** Error in getgauge: cannot move to outdir = ',\
                       outdir
                print '*** thisdir = ',thisdir
                raise
                return
            try:
                gauges = self.read_gauges(outdir)
            except:
                print '*** Error reading gauges in ClawPlotData.getgauge'
                print '*** outdir = ', outdir
                print '*** thisdir = ', thisdir
                os.chdir(thisdir)
                raise
                return
            os.chdir(thisdir)

            try:
                for (k,v) in gauges.iteritems():
                    gaugesoln_dict[(k, outdir)] = v
            except:
                raise("*** Problem setting gaugesoln_dict in getgauge")

            #print '    Read all gauge data from %s/fort.gauge' % outdir

        try:
            gaugesoln = gaugesoln_dict[key]
        except:
            print "*** Cannot find key = ",key
            print "***   in gaugesoln_dict = ",gaugesoln_dict
            raise("*** Problem getting gaugesoln in getgauge")
                

        return gaugesoln

    def read_gauges(self, outdir='.'):
        """
        Read the gauge output in file fort.gauge in the directory specified by
        outdir.
    
        Returns a dictionary *gauges* with an entry for each gauge number.
        Each entry is an object of class GaugeSolution
    
        """
    
        import os
        import numpy as np
        from matplotlib.mlab import find
        from petclaw.plotters import gaugetools
        from StringIO import StringIO
    
        fname = outdir + '/fort.gauge'
        if not os.path.isfile(fname):
            print "*** Gauge file not found: ",fname
            gauges = {}

        print '    Reading gauge data from ',fname
        try:
            gdata = np.loadtxt(fname)
        except:
            try:
                print "*** Warning: incomplete last line, computation may "
                print "*** still be in progress "
                gdata_lines = open(fname,'r').read()
                gdata_end = gdata_lines.rfind('\n',-200,-1)
                gdata_file = StringIO(gdata_lines[:gdata_end+1])
                gdata = np.loadtxt(gdata_file)
            except:
                print "*** Problem reading file ",fname
                #print "*** Possibly an incomplete last line if computation is still in progress"
                raise Exception("Problem reading fort.gauge")
                gauges = {}

        gaugeno = np.array(gdata[:,0], dtype=int)
        level = np.array(gdata[:,1], dtype=int)
        t = gdata[:,2]
        q = gdata[:,3:]  # all remaining columns are stored in q
    
        
        setgauges = gaugetools.read_setgauges(datadir=outdir)
    
        gauges = {}
        gaugenos = set(gaugeno)   # reduces to unique elements
        for n in gaugenos:
            n = int(n)
            gauges[n] = GaugeSolution()
            gauges[n].gaugeno = n
            nn = find(gaugeno==n)
            gauges[n].level = level[nn]
            gauges[n].t = t[nn]
            gauges[n].q = q[nn,:]
    
            # Locations:
            try:
                gauges[n].x = setgauges.x[n]
                gauges[n].y = setgauges.y[n]
                gauges[n].t1 = setgauges.t1[n]
                gauges[n].t2 = setgauges.t2[n]
            except:
                print "*** Could not extract gauge locations for gaugeno = ",n
    
        print '    Found gauge numbers: ',gauges.keys()
        return gauges
    


    def plotframe(self, frameno):
        from petclaw.plotters import frametools
        frametools.plotframe(frameno, self)
        
    def printframes(self, verbose=True):
        #from petclaw.plotters import frametools
        #frametools.printframes(self, verbose)
        print "*** printframes is deprecated.  Use plotpages.plotclaw_driver"
        print "*** for added capabilities."
        
    def fignos(self):
        """
        Return a list of the figure numbers actually used.
        Useful in afterframe function for example to loop over all
        figures and do something.
        """
        return self._fignos

    def mode(self):
        """
        Return self._mode, which is set internally to 
           'iplotclaw' if Iplotclaw is in use,
           'printframes' if printframes is being used
        Useful in afterframe function if you want to do different things
           for interactive or print modes.
        """
        return self._mode

    def iplotclaw(self):
        """
        Return True if interactive plotting with iplotclaw is being done.
        """
        return (self._mode == 'iplotclaw')


    def getfigure(self,figname):
        try:
            plotfigure = self.plotfigure_dict[figname]
        except:
            print '*** Error accessing plotfigure_dict[%s]' % figname
            return None
        return plotfigure

    def getaxes(self,axesname,figname=None):
        found = True
        if not figname:
            found = False
            for fig in self._fignames:
                plotfigure = self.getfigure(fig)
                if axesname in plotfigure._axesnames:
                    if found == True: # already found!
                        print '*** Ambiguous... must specify figname'
                        print '    try getaxes(axesname, figname)'
                        return None
                    figname = fig
                    found = True
        if not found:
            print '*** No axes found with name = ',axesname
            return None
        try:
            plotfigure = self.getfigure(figname)
            plotaxes = plotfigure.plotaxes_dict[axesname]
        except:
            print '*** Error accessing plotaxes[%s]' % axesname
            print '*** figname = %s' % figname
            return None
        return plotaxes

    def getitem(self,itemname,axesname=None,figname=None):
        found = True
        if not figname:
            # search over all figures looking for the item
            found = False
            for fign in self._fignames:
                plotfigure = self.getfigure(fign)
                if not axesname:
                    # search over all axes looking for the item
                    for axesn in plotfigure._axesnames:
                        plotaxes = self.getaxes(axesn,fign)
                        if itemname in plotaxes._itemnames:
                            if found == True: # already found!
                                print '*** Ambiguous... must specify figname and/or axesname'
                                print '    try getitem(itemname, axesname, figname)'
                                return None
                            axesname = axesn
                            figname = fign
                            found = True
                else:
                    # axesname was specified (but not figname)
                    plotaxes = self.getaxes(axesname,fign)
                    if itemname in plotaxes._itemnames:
                        if found == True: # already found!
                            print '*** Ambiguous... must specify figname and/or axesname'
                            print '    try getitem(itemname, axesname, figname)'
                            return None
                        figname = fign
                        found = True

        elif not axesname:
            # figname was specified but not axesname.
            # search over all axes looking for the item
            found = False
            plotfigure = self.getfigure(figname)
            for axesn in plotfigure._axesnames:
                plotaxes = self.getaxes(axesn,figname)
                if itemname in plotaxes._itemnames:
                    if found == True: # already found!
                        print '*** Ambiguous... must specify axesname'
                        print '    try getitem(itemname, axesname, figname)'
                        return None
                    axesname = axesn
                    found = True

        if not found:
            print '*** No item found with name = ',itemname
            return None
        try:
            plotaxes = self.getaxes(axesname,figname)
            plotitem = plotaxes.plotitem_dict[itemname]
        except:
            print '*** Error accessing plotitem[%s]' % itemname
            print '*** figname = ',figname
            print '*** axesname = ',axesname
            return None
        return plotitem


    def showitems(self):
        fignames = self._fignames
        print "\n\nCurrent plot figures, axes, and items:"
        print "---------------------------------------"
        for figname in fignames:
            plotfigure = self.getfigure(figname)
            s =  "  figname = %s, figno = %s" % (figname, plotfigure.figno)
            if not plotfigure._show: 
                s = s + "  [Not showing]"
            print s
            axesnames = plotfigure._axesnames
            for axesname in axesnames:
                plotaxes = self.getaxes(axesname,figname)
                s =  "     axesname = %s, axescmd = %s" \
                       % (axesname, plotaxes.axescmd)
                if not plotaxes._show: 
                    s = s + "  [Not showing]"
                print s
                for itemname in plotaxes._itemnames:
                    plotitem = self.getitem(itemname,axesname,figname)
                    plot_type = plotitem.plot_type
                    s =  "        itemname = %s,  plot_type = %s" \
                          % (itemname,plot_type)
                    if not plotitem._show: 
                        s = s + "  [Not showing]"
                    print s
            print " "


    def getq(self,frameno):
        solution = self.getframe(frameno)
        grids = solution.grids
        if len(grids) > 1:
            print '*** Warning: more than 1 grid, q on grid[0] is returned'
        q = grids[0].q
        return q


# ============================================================================
#  Subclass ClawPlotFigure containing data for plotting a figure
# ============================================================================
class ClawPlotFigure(Data):
    """
    
    Data subclass containing plot data needed to plot a single figure.
    This may consist of several ClawPlotAxes objects.

    """

    # ========================================================================
    #  Initialization routine
    # ========================================================================
    def __init__(self, name, figno, type, plotdata):
        """
        Initialize a ClawPlotFigure object
        """

        attributes = ['name','figno','_plotdata','clf','plotaxes_dict', \
                           '_axesnames','show','_show','kwargs','_handle',\
			   '_type']

        super(ClawPlotFigure, self).__init__(attributes = attributes)    

        self._plotdata = plotdata           # parent ClawPlotData object
        self.name = name
        self.figno = figno
        self.kwargs = {}
        self.clf_each_frame = True
        self.clf_each_gauge = True
        self._axesnames = []
        self.show = True
        self._show = True
        self.plotaxes_dict  = {}
        self.type = type   # = 'each_frame' or 'each_run' or 'each_gauge'
        self._next_AXES = 0

    def new_plotaxes(self, name=None, type='each_frame'):
        """
        Create a new axes that will be plotted in this figure.
        If type='each_frame' it is an axes that will be plotted 
	for each time frame.
        If type='multi_frame' it is an axes that will be plotted based on
	all the frames, such as x-t plots or time series. (Not yet implemented)
        If type='empty' it is created without doing any plots using the
        petclaw tools.  Presumably the user will create a plot within an
        afteraxes command, for example.
        """
        if name is None:
            self._next_AXES += 1
            name = "AXES%s" % self._next_AXES
        if name in self._axesnames:
            print '*** Warning, axes named %s has already been created' % name

        if name not in self._axesnames:
            self._axesnames.append(name)
        plotaxes = ClawPlotAxes(name, self)
        self.plotaxes_dict[name] = plotaxes
        plotaxes.type = type
        return plotaxes

    def gethandle(self):
        _handle = getattr(self,'_handle',None)
        return _handle


# ============================================================================
#  Subclass ClawPlotAxes containing data for plotting axes within a figure
# ============================================================================
class ClawPlotAxes(Data):
    """
    
    Data subclass containing plot data needed to plot a single axes.
    This may consist of several ClawPlotItem objects.

    """

    # ========================================================================
    #  Initialization routine
    # ========================================================================
    def __init__(self, name, plotfigure):
        """
        Initialize a ClawPlotAxes object
        """

        attributes = ['name','type','figno','plotdata','plotfigure','title',\
                           'axescmd','xlimits','ylimits','plotitem_dict', 'user',\
                           'afteraxes','_itemnames','show','_show','_handle', \
                           '_plotfigure','_plotdata', 'scaled']
        super(ClawPlotAxes, self).__init__(attributes = attributes)    

        self._plotfigure = plotfigure                   # figure this item is on
        self._plotdata = plotfigure._plotdata           # parent ClawPlotData object

        self.name = name
        self.title = name
        self.title_with_t = True        # creates title of form 'title at time t = ...'
        self.axescmd = 'subplot(1,1,1)'
        self.user = Data()          # for user to pass things into
                                        # afteraxes, for example
					# Deprecated.
        self.afteraxes = None
        self.xlimits = None
        self.ylimits = None
        self.scaled = False              # true so x- and y-axis scaled same
        self.plotitem_dict  = {}
        self.type = 'each_frame'
        self._itemnames = []
        self.show = True
        self._show = True
        self._handle = None
        self._next_ITEM = 0
        self.figno = self._plotfigure.figno

    def new_plotitem(self, name=None, plot_type=None):
        # Create a new entry in self.plotitem_dict

        if name is None:
            self._next_ITEM += 1
            name = "ITEM%s" % self._next_ITEM
        
        if name not in self._itemnames:
            self._itemnames.append(name)

        plotitem = ClawPlotItem(name, plot_type, plotaxes=self)

        self.plotitem_dict[name] = plotitem
        
        return plotitem

    def get_plotdata(self):
        plotdata = getattr(self,'_plotdata',None)
        return self._plotdata

    def get_plotfigure(self):
        plotfigure = getattr(self,'_plotfigure',None)
        return self._plotfigure

    def gethandle(self):
        _handle = getattr(self,'_handle',None)
        return self._handle

# ============================================================================
#  Subclass ClawPlotItem containing data for plotting a single object
# ============================================================================
class ClawPlotItem(Data):
    """
    
    Data subclass containing plot data needed to plot a single object.
    This may be a single curve, set of points, contour plot, etc.

    """

    # ========================================================================
    #  Initialization routine
    # ========================================================================
    def __init__(self, name, plot_type, plotaxes):
        """
        Initialize a ClawPlotItem object
        """
        attributes = ['ndim','outdir','refresh_frames',\
                           'plot_var','plot_var_title', \
                           'MappedGrid', 'mapc2p', \
                           'figno', 'handle', 'params', \
                           'aftergrid','afteritem','framesoln_dict', \
                           '_pobjs']  

        super(ClawPlotItem, self).__init__(attributes = attributes)    


        self._plotaxes = plotaxes                       # axes this item is on
        self._plotfigure = plotaxes._plotfigure         # figure this item is on
        self._plotdata = plotaxes._plotfigure._plotdata           # parent ClawPlotData object

        try:
            ndim = int(plot_type[0])   # first character of plot_type should be ndim
        except:
            print '*** Error: could not determine ndim from plot_type = ',plot_type

        self.ndim = ndim
        self.name = name
        self.figno = plotaxes.figno

        self.outdir = None              # indicates data comes from 
                                        #   self._plotdata.outdir

        self.plot_type = plot_type
        self.plot_var = 0
        self.plot_show = True

        self.MappedGrid = None          # False to plot on comput. grid even
                                        # if _plotdata.mapc2p is not None.

        self.mapc2p = None              # function to map computational
	                                # points to physical (over-rides
	                                # plotdata.mapc2p if set for item


        self.aftergrid = None           # function called after each grid is
                                        # plotted within each single plotitem.
        self.afteritem = None           # function called after the item is
                                        # plotted for each frame

        self.user = Data()              # for user to pass things into
                                        # aftergrid, for example
                                        # Deprecated.

        self.show = True                # False => suppress showing this item
        self._show = True               # Internal

        self._current_pobj = None


        if ndim == 1:
            self.add_attribute('plotstyle','-')
            self.add_attribute('color',None)
            self.add_attribute('kwargs',{})

            if plot_type == '1d_fill_between':
                zero_function = lambda current_data: 0.
                self.add_attribute('plot_var2',zero_function)
                self.add_attribute('fill_where',None)

            if plot_type == '1d_from_2d_data':
                self.add_attribute('map2d_to_1d',None)

        elif ndim == 2:

            # default values specifying this single plot:
            self.add_attribute('plot_type',plot_type)
            self.add_attribute('gridlines_show',0)
            self.add_attribute('gridlines_color','k')
            self.add_attribute('grid_bgcolor','w')
            self.add_attribute('gridedges_show',0)
            self.add_attribute('gridedges_color','k')
            self.add_attribute('kwargs',{})

            if plot_type == '2d_pcolor':
                # from pylab import cm
                # self.add_attribute('pcolor_cmap',cm.RdYlBu,True)
                from petclaw.plotters import colormaps
                self.add_attribute('pcolor_cmap',colormaps.yellow_red_blue)
                self.add_attribute('pcolor_cmin',None)
                self.add_attribute('pcolor_cmax',None)
                self.add_attribute('add_colorbar',True)

            elif plot_type == '2d_imshow':
                # from pylab import cm
                # self.add_attribute('pcolor_cmap',cm.RdYlBu,True)
                from petclaw.plotters import colormaps
                self.add_attribute('imshow_cmap',colormaps.yellow_red_blue)
                self.add_attribute('imshow_cmin',None)
                self.add_attribute('imshow_cmax',None)
                self.add_attribute('add_colorbar',True)


            elif plot_type == '2d_contour':
                self.add_attribute('contour_nlevels',20)
                self.add_attribute('contour_levels',None)
                self.add_attribute('contour_min',None)
                self.add_attribute('contour_max',None)
                self.add_attribute('contour_show',1)
                self.add_attribute('contour_color','k')
                self.add_attribute('contour_cmap',None)
                self.add_attribute('add_colorbar',False)

            elif plot_type == '2d_schlieren':
                from petclaw.plotters import colormaps
                self.add_attribute('schlieren_cmap',colormaps.schlieren_grays)
                self.add_attribute('schlieren_cmin',None)
                self.add_attribute('schlieren_cmax',None)
                self.add_attribute('add_colorbar',False)

            elif plot_type == '2d_grid':
                self.add_attribute('max_density',None)
                self.gridlines_show = True
                
            elif plot_type == '2d_quiver':
                self.add_attribute('quiver_var_x',None)
                self.add_attribute('quiver_var_y',None)
                self.add_attribute('quiver_coarsening',1)
                self.add_attribute('quiver_key_show',False)
                self.add_attribute('quiver_key_label_x',0.15)
                self.add_attribute('quiver_key_label_y',0.95)
                self.add_attribute('quiver_key_units','')
                self.add_attribute('quiver_key_scale',None)
                self.add_attribute('quiver_key_kwargs',{})

            else:
                 print '*** Warning 2d plot type %s not recognized' % plot_type

        elif ndim == 3:
            print '*** Warning- ClawPlotItem not yet set up for ndim = 3'
    
        else:
            print '*** Warning- Unrecognized plot_type in ClawPlotItem'

        self.params = {}  # dictionary to hold optional parameters
        

    def getframe(self,frameno):
        """
        ClawPlotItem.getframe:
        Return an object of class Solution containing the solution
        for frame number frameno.

        If self.refresh_frames == True then this frame is read from the fort
        files, otherwise it is read from the fort files only if the
        the dictionary self.framesoln_dict has key frameno.  If it does, the
        frame has previously been read and the dictionary value is returned.
        """

        plotdata = self._plotdata
        outdir = self.outdir
        framesoln = plotdata.getframe(frameno, outdir)

        return framesoln


    def getgauge(self,gauge):
        """
        ClawPlotItem.getgauge:
        Return an object of class GaugeSolution containing the solution
        for gauge number gaugeno.

        If self.refresh_gauges == True then this gauge is read from the
        fort.gauge file, otherwise it is read only if the
        the dictionary self.gaugesoln_dict has no key gaugeno.  If it does, the
        gauge has previously been read and the dictionary value is returned.
        """

        plotdata = self._plotdata
        outdir = self.outdir
        gaugesoln = plotdata.getgauge(gauge, outdir)

        return gaugesoln

#-----------------------------------------------------------------------
# New classes and functions for dealing with data in setrun function.

class ClawInputData(Data):
    """
    Object that will be written out to claw.data.
    """
    def __init__(self, ndim):
        super(ClawInputData,self).__init__()
        self.add_attribute('ndim',ndim)
        
        # Set default values:
        if ndim == 1:
            self.add_attribute('mx',100)
            self.add_attribute('nout',5)
            self.add_attribute('outstyle',1)
            self.add_attribute('tfinal',1.0)
            self.add_attribute('dt_initial',1.e-5)
            self.add_attribute('dt_max',1.e99)
            self.add_attribute('cfl_desired',0.9)
            self.add_attribute('cfl_max',1.0)
            self.add_attribute('max_steps',5000)
            self.add_attribute('dt_variable',1)
            self.add_attribute('order',2)
            self.add_attribute('order_trans',0)
            self.add_attribute('verbosity',0)
            self.add_attribute('src_split',0)
            self.add_attribute('mcapa',0)
            self.add_attribute('maux',0)
            self.add_attribute('meqn',1)
            self.add_attribute('mwaves',1)
            self.add_attribute('mthlim',[4])
            self.add_attribute('t0',0.)
            self.add_attribute('xlower',0.)
            self.add_attribute('xupper',1.)
            self.add_attribute('mbc',2)
            self.add_attribute('mthbc_xlower',1)
            self.add_attribute('mthbc_xupper',1)
            self.add_attribute('restart',0)
            self.add_attribute('N_restart',0)


        elif ndim == 2:
            self.add_attribute('mx',100)
            self.add_attribute('my',100)
            self.add_attribute('nout',5)
            self.add_attribute('outstyle',1)
            self.add_attribute('tfinal',1.0)
            self.add_attribute('dt_initial',1.e-5)
            self.add_attribute('dt_max',1.e99)
            self.add_attribute('cfl_desired',0.9)
            self.add_attribute('cfl_max',1.0)
            self.add_attribute('max_steps',5000)
            self.add_attribute('dt_variable',1)
            self.add_attribute('order',2)
            self.add_attribute('order_trans',2)
            self.add_attribute('verbosity',0)
            self.add_attribute('src_split',0)
            self.add_attribute('mcapa',0)
            self.add_attribute('maux',0)
            self.add_attribute('meqn',1)
            self.add_attribute('mwaves',1)
            self.add_attribute('mthlim',[4])
            self.add_attribute('t0',0.)
            self.add_attribute('xlower',0.)
            self.add_attribute('xupper',1.)
            self.add_attribute('ylower',0.)
            self.add_attribute('yupper',1.)
            self.add_attribute('mbc',2)
            self.add_attribute('mthbc_xlower',1)
            self.add_attribute('mthbc_xupper',1)
            self.add_attribute('mthbc_ylower',1)
            self.add_attribute('mthbc_yupper',1)
            self.add_attribute('restart',0)
            self.add_attribute('N_restart',0)

        else:
            print '*** Error: only ndim=1 or 2 supported so far ***'
            raise()

    def write(self):
        print 'Creating data file claw.data for use with xclaw'
        make_clawdatafile(self)   



class AmrclawInputData(Data):
    """
    Object that will be written out to amr2ez.data.
    """
    def __init__(self, ndim):
        super(AmrclawInputData,self).__init__()
        self.add_attribute('ndim',ndim)
        
        # Set default values:
        if ndim == 1:
            self.add_attribute('mx',100)
            self.add_attribute('nout',5)
            self.add_attribute('outstyle',1)
            self.add_attribute('tfinal',1.0)
            self.add_attribute('dt_initial',1.e-5)
            self.add_attribute('dt_max',1.e99)
            self.add_attribute('cfl_desired',0.9)
            self.add_attribute('cfl_max',1.0)
            self.add_attribute('max_steps',5000)
            self.add_attribute('dt_variable',1)
            self.add_attribute('order',2)
            self.add_attribute('order_trans',0)
            self.add_attribute('verbosity',0)
            self.add_attribute('src_split',0)
            self.add_attribute('mcapa',0)
            self.add_attribute('maux',0)
            self.add_attribute('meqn',1)
            self.add_attribute('mwaves',1)
            self.add_attribute('mthlim',[4])
            self.add_attribute('t0',0.)
            self.add_attribute('xlower',0.)
            self.add_attribute('xupper',1.)
            self.add_attribute('mbc',2)
            self.add_attribute('mthbc_xlower',1)
            self.add_attribute('mthbc_xupper',1)
            self.add_attribute('restart',0)
            self.add_attribute('N_restart',0)

            # attributes need only since AMR is done using 2d amrclaw:
            self.add_attribute('my',1)
            self.add_attribute('ylower',0.)
            self.add_attribute('yupper',1.)
            self.add_attribute('mthbc_ylower',1)
            self.add_attribute('mthbc_yupper',1)
            self.add_attribute('inraty',[1,1,1,1,1,1])

        elif ndim == 2:
            self.add_attribute('mx',100)
            self.add_attribute('my',100)
            self.add_attribute('nout',5)
            self.add_attribute('outstyle',1)
            self.add_attribute('tfinal',1.0)
            self.add_attribute('dt_initial',1.e-5)
            self.add_attribute('dt_max',1.e99)
            self.add_attribute('cfl_desired',0.9)
            self.add_attribute('cfl_max',1.0)
            self.add_attribute('max_steps',5000)
            self.add_attribute('dt_variable',1)
            self.add_attribute('order',2)
            self.add_attribute('order_trans',2)
            self.add_attribute('verbosity',0)
            self.add_attribute('src_split',0)
            self.add_attribute('mcapa',0)
            self.add_attribute('maux',0)
            self.add_attribute('meqn',1)
            self.add_attribute('mwaves',1)
            self.add_attribute('mthlim',[4])
            self.add_attribute('t0',0.)
            self.add_attribute('xlower',0.)
            self.add_attribute('xupper',1.)
            self.add_attribute('ylower',0.)
            self.add_attribute('yupper',1.)
            self.add_attribute('mbc',2)
            self.add_attribute('mthbc_xlower',1)
            self.add_attribute('mthbc_xupper',1)
            self.add_attribute('mthbc_ylower',1)
            self.add_attribute('mthbc_yupper',1)
            self.add_attribute('restart',0)
            self.add_attribute('N_restart',0)
            self.add_attribute('inraty',[1])

        if ndim <= 2:
            # AMR parameters:
            self.add_attribute('mxnest',-1)
            self.add_attribute('inratx',[1])
            self.add_attribute('inratt',[1])
            self.add_attribute('auxtype',[])
            self.add_attribute('restart',False)
            self.add_attribute('checkpt_iousr',1000)
            self.add_attribute('tol',-1.0)
            self.add_attribute('tolsp',0.05)
            self.add_attribute('kcheck',2)
            self.add_attribute('ibuff',3)
            self.add_attribute('cutoff',0.7)
            self.add_attribute('PRINT',False)
            self.add_attribute('NCAR',False)
            self.add_attribute('fortq',True)
            self.add_attribute('dprint',False)
            self.add_attribute('eprint',False)
            self.add_attribute('edebug',False)
            self.add_attribute('gprint',False)
            self.add_attribute('nprint',False)
            self.add_attribute('pprint',False)
            self.add_attribute('rprint',False)
            self.add_attribute('sprint',False)
            self.add_attribute('tprint',False)
            self.add_attribute('uprint',False)
        else:
            print '*** Error: only ndim=1 or 2 supported so far ***'
            raise()

    def write(self):
        print 'Creating data file amr2ez.data for use with xamr'
        make_amrclawdatafile(self)   



def open_datafile(name, datasource='setrun.py'):
    """
    Open a data file and write a warning header.
    Warning header starts with '#' character.  These lines are skipped if
    data file is opened using the library routine opendatafile.

    INPUT:
        name - name of data file
    OUTPUT:
        file - file object
    """
    
    import string

    source = string.ljust(datasource,25)
    file = open(name, 'w')
    file.write('########################################################\n')
    file.write('### DO NOT EDIT THIS FILE:  GENERATED AUTOMATICALLY ####\n')
    file.write('### To modify data, edit  %s ####\n' % source)
    file.write('###    and then "make .data"                        ####\n')
    file.write('########################################################\n\n')

    return file


def data_write(file, dataobj, name=None, descr=''):

    """
    Write out value to data file, in the form
       value =: name  descr
    Remove brackets and commas from lists, and replace booleans by T/F.

    INPUTS
       name, normally a string defining the variable
             if name==None, write a blank line.
       descr, a short description to appear on the line
    """

    import string
    if name is None:
        file.write('\n')
    else:
        try:
            value = getattr(dataobj, name)
        except:
            print "Variable missing: ",name
            print "  from dataobj = ", dataobj
            raise
        # Convert value to an appropriate string repr
        if isinstance(value,tuple) | isinstance(value,list):
            # Remove [], (), and ','
            string_value = repr(value)[1:-1]
            string_value = string_value.replace(',','')
        elif isinstance(value,bool):
            if value:
                string_value = 'T'
            else:
                string_value = 'F'
        else:
            string_value = repr(value)
        padded_value = string.ljust(string_value, 25)
        padded_name = string.ljust(name, 12)
        file.write('%s =: %s %s\n' % (padded_value, padded_name, descr))
        

def make_clawdatafile(clawdata):
    """
    Take the data specified in clawdata and write it to claw.data in the
    form required by the Fortran code lib/main.f95.
    """


    # open file and write a warning header:
    file = open_datafile('claw.data')

    ndim = clawdata.ndim
    data_write(file, clawdata, 'ndim', '(number of dimensions)')
    data_write(file, clawdata, 'mx', '(cells in x direction)')
    if ndim > 1:
        data_write(file, clawdata, 'my', '(cells in y direction)')
    if ndim == 3:
        data_write(file, clawdata, 'mz', '(cells in z direction)')
    data_write(file, clawdata, None)  # writes blank line

    data_write(file, clawdata, 'nout', '(number of output times)')
    data_write(file, clawdata, 'outstyle', '(style of specifying output times)')
    if clawdata.outstyle == 1:
        data_write(file, clawdata, 'tfinal', '(final time)')
    elif clawdata.outstyle == 2:
        data_write(file, clawdata, 'tout', '(output times)')
    elif clawdata.outstyle == 3:
        data_write(file, clawdata, 'iout', '(output every iout steps)')
    else:
        print '*** Error: unrecognized outstyle'
        raise
        return

    data_write(file, clawdata, None)
    data_write(file, clawdata, 'dt_initial', '(initial time step dt)')
    data_write(file, clawdata, 'dt_max', '(max allowable dt)')
    data_write(file, clawdata, 'cfl_max', '(max allowable Courant number)')
    data_write(file, clawdata, 'cfl_desired', '(desired Courant number)')
    data_write(file, clawdata, 'max_steps', '(max time steps per call to claw)')
    data_write(file, clawdata, None)
    data_write(file, clawdata, 'dt_variable', '(1 for variable dt, 0 for fixed)')
    data_write(file, clawdata, 'order', '(1 or 2)')
    if ndim == 1:
        data_write(file, clawdata, 'order_trans', '(not used in 1d)')
    else:
        data_write(file, clawdata, 'order_trans', '(transverse order)')
    data_write(file, clawdata, 'verbosity', '(verbosity of output)')
    data_write(file, clawdata, 'src_split', '(source term splitting)')
    data_write(file, clawdata, 'mcapa', '(aux index for capacity fcn)')
    data_write(file, clawdata, 'maux', '(number of aux variables)')
    data_write(file, clawdata, None)
    
    data_write(file, clawdata, 'meqn', '(number of equations)')
    data_write(file, clawdata, 'mwaves', '(number of waves)')
    data_write(file, clawdata, 'mthlim', '(limiter choice for each wave)')
    data_write(file, clawdata, None)
    
    data_write(file, clawdata, 't0', '(initial time)')
    data_write(file, clawdata, 'xlower', '(xlower)')
    data_write(file, clawdata, 'xupper', '(xupper)')
    if ndim > 1:
        data_write(file, clawdata, 'ylower', '(ylower)')
        data_write(file, clawdata, 'yupper', '(yupper)')
    if ndim == 3:
        data_write(file, clawdata, 'zlower', '(zlower)')
        data_write(file, clawdata, 'zupper', '(zupper)')
    data_write(file, clawdata, None)
    
    data_write(file, clawdata, 'mbc', '(number of ghost cells)')
    data_write(file, clawdata, 'mthbc_xlower', '(type of BC at xlower)')
    data_write(file, clawdata, 'mthbc_xupper', '(type of BC at xupper)')
    if ndim > 1:
        data_write(file, clawdata, 'mthbc_ylower', '(type of BC at ylower)')
        data_write(file, clawdata, 'mthbc_yupper', '(type of BC at yupper)')
    if ndim == 3:
        data_write(file, clawdata, 'mthbc_zlower', '(type of BC at zlower)')
        data_write(file, clawdata, 'mthbc_zupper', '(type of BC at zupper)')
    
    file.close()

def make_amrclawdatafile(clawdata):
    """
    Take the data specified in clawdata and write it to claw.data in the
    form required by the Fortran code lib/main.f95.
    """


    # open file and write a warning header:
    file = open_datafile('amr2ez.data')

    ndim = clawdata.ndim
    #data_write(file, clawdata, 'ndim', '(number of dimensions)')
    data_write(file, clawdata, 'mx', '(cells in x direction)')
    data_write(file, clawdata, 'my', '(cells in y direction)')
    if ndim == 3:
        data_write(file, clawdata, 'mz', '(cells in z direction)')
    data_write(file, clawdata, 'mxnest', '(max number of grid levels)')
    data_write(file, clawdata, 'inratx', '(refinement ratios)')
    data_write(file, clawdata, 'inraty', '(refinement ratios)')
    if ndim == 3:
        data_write(file, clawdata, 'inratz', '(refinement ratios)')
    data_write(file, clawdata, 'inratt', '(refinement ratios)')
    data_write(file, clawdata, None)  # writes blank line

    data_write(file, clawdata, 'nout', '(number of output times)')
    data_write(file, clawdata, 'outstyle', '(style of specifying output times)')
    if clawdata.outstyle == 1:
        data_write(file, clawdata, 'tfinal', '(final time)')
    elif clawdata.outstyle == 2:
        data_write(file, clawdata, 'tout', '(output times)')
    elif clawdata.outstyle == 3:
        data_write(file, clawdata, 'iout', '(output every iout steps)')
    else:
        print '*** Error: unrecognized outstyle'
        raise
        return

    data_write(file, clawdata, None)
    data_write(file, clawdata, 'dt_initial', '(initial time step dt)')
    data_write(file, clawdata, 'dt_max', '(max allowable dt)')
    data_write(file, clawdata, 'cfl_max', '(max allowable Courant number)')
    data_write(file, clawdata, 'cfl_desired', '(desired Courant number)')
    data_write(file, clawdata, 'max_steps', '(max time steps per call to claw)')
    data_write(file, clawdata, None)
    data_write(file, clawdata, 'dt_variable', '(1 for variable dt, 0 for fixed)')
    data_write(file, clawdata, 'order', '(1 or 2)')
    if ndim == 1:
        data_write(file, clawdata, 'order_trans', '(not used in 1d)')
    else:
        data_write(file, clawdata, 'order_trans', '(transverse order)')
    data_write(file, clawdata, 'verbosity', '(verbosity of output)')
    data_write(file, clawdata, 'src_split', '(source term splitting)')
    data_write(file, clawdata, 'mcapa', '(aux index for capacity fcn)')
    data_write(file, clawdata, 'maux', '(number of aux variables)')
    if len(clawdata.auxtype) != clawdata.maux:
        file.close()
        raise AttributeError, "require len(clawdata.auxtype) == clawdata.maux"
    for i in range(clawdata.maux):
        file.write("'%s'\n" % clawdata.auxtype[i])
    data_write(file, clawdata, None)
    
    data_write(file, clawdata, 'meqn', '(number of equations)')
    data_write(file, clawdata, 'mwaves', '(number of waves)')
    data_write(file, clawdata, 'mthlim', '(limiter choice for each wave)')
    data_write(file, clawdata, None)
    
    data_write(file, clawdata, 't0', '(initial time)')
    data_write(file, clawdata, 'xlower', '(xlower)')
    data_write(file, clawdata, 'xupper', '(xupper)')
    data_write(file, clawdata, 'ylower', '(ylower)')
    data_write(file, clawdata, 'yupper', '(yupper)')
    if ndim == 3:
        data_write(file, clawdata, 'zlower', '(zlower)')
        data_write(file, clawdata, 'zupper', '(zupper)')
    data_write(file, clawdata, None)
    
    data_write(file, clawdata, 'mbc', '(number of ghost cells)')
    data_write(file, clawdata, 'mthbc_xlower', '(type of BC at xlower)')
    data_write(file, clawdata, 'mthbc_xupper', '(type of BC at xupper)')
    data_write(file, clawdata, 'mthbc_ylower', '(type of BC at ylower)')
    data_write(file, clawdata, 'mthbc_yupper', '(type of BC at yupper)')
    if ndim == 3:
        data_write(file, clawdata, 'mthbc_zlower', '(type of BC at zlower)')
        data_write(file, clawdata, 'mthbc_zupper', '(type of BC at zupper)')
    data_write(file, clawdata, None)

    data_write(file, clawdata, 'restart', '(1 to restart from a past run)')
    data_write(file, clawdata, 'checkpt_iousr', '(how often to checkpoint)')
    data_write(file, clawdata, None)

    data_write(file, clawdata, 'tol', '(tolerance for Richardson extrap)')
    data_write(file, clawdata, 'tolsp', '(tolerance used in flag2refine)')
    data_write(file, clawdata, 'kcheck', '(how often to regrid)')
    data_write(file, clawdata, 'ibuff', '(buffer zone around flagged pts)')
    data_write(file, clawdata, 'cutoff', '(efficiency cutoff for grid gen.)')
    data_write(file, clawdata, None)

    data_write(file, clawdata, 'PRINT', '(print to fort.amr)')
    data_write(file, clawdata, 'NCAR', '(obsolete!)')
    data_write(file, clawdata, 'fortq', '(Output to fort.q* files)')
    data_write(file, clawdata, None)

    data_write(file, clawdata, 'dprint', '(print domain flags)')
    data_write(file, clawdata, 'eprint', '(print err est flags)')
    data_write(file, clawdata, 'edebug', '(even more err est flags)')
    data_write(file, clawdata, 'gprint', '(grid bisection/clustering)')
    data_write(file, clawdata, 'nprint', '(proper nesting output)')
    data_write(file, clawdata, 'pprint', '(proj. of tagged points)')
    data_write(file, clawdata, 'rprint', '(print regridding summary)')
    data_write(file, clawdata, 'sprint', '(space/memory output)')
    data_write(file, clawdata, 'tprint', '(time step reporting each level)')
    data_write(file, clawdata, 'uprint', '(update/upbnd reporting)')
    
    file.close()

def make_userdatafile(userdata):

    """
    Create the data file using the parameters in userdata.
    The parameters will be written to this file in the same order they were
    specified using userdata.add_attribute.
    Presumably the user will read these in using a Fortran routine, such as
    setprob.f95, and the order is important.
    """

    # open file and write a warning header:
    file = open_datafile(userdata._UserData__fname)

    # write all the parameters:
    for param in userdata.attributes:
        data_write(file, userdata, param, \
                   userdata._UserData__descr[param])

    file.close()

class GaugeSolution(Data):
    """
    Holds gaugeno, t, q, x, y, t1, t2 for a single gauge.
    """

    def __init__(self):

        data_files = []
        gauge_attrs = 'gaugeno level t q x y t1 t2'.split()

        # Initialize the data object and read the data files
        super(GaugeSolution,self).__init__(data_files,gauge_attrs)

        # default values of attributes:
        # none.

