"""
Module Iplotclaw for interactive plotting of Clawpack results.

For more instructions see the Usage notes in class Iplotclaw below. 

For options during looping type:
  >>> ip = Iplotclaw()
  >>> ip.plotloop()
  PLOTCLAW> help

"""

import cmd, os, sys, string


if not sys.modules.has_key('matplotlib'):
    try:
        import matplotlib
        # Override system defaults before importing pylab
        # If you set the backend here, you must start ipython w/o pylab
        # before importing this package.
        #matplotlib.use('Agg')  # Use an image backend
        #matplotlib.use('TkAgg')
        matplotlib.rc('text', usetex=False)
        matplotlib.interactive(True)

    except:
        print "problem importing matplotlib"
        sys.exit(1)

try:
    import pylab 
except:
    print "problem importing pylab"
    sys.exit(1)

from pyclaw.plotters import data, frametools, gaugetools

#rundata = data.ClawData('claw.data')

#------------------------
class Iplotclaw(cmd.Cmd):
#------------------------
    """
    Class for interactively stepping through Clawpack plots.
    Uses frametools.plotframe to actually plot each frame.

    Usage:
    ------
    >>> from pyclaw.plotters.Iplotclaw import Iplotclaw 
    >>> ip = Iplotclaw()             # new instantiation
    >>> ip.plotloop()                # to start looping
    PLOTCLAW > help                  # for list of available commands
    PLOTCLAW > q                     # to quit looping and return to python
    >>> ip.plotloop()                # to restart looping at previous frame,
                                     #    with data, outdir, etc. preserved.
    Arguments:
    ----------
    The defaults are:
       Iplotclaw(setplot='setplot.py', outdir=None) 
    The setplot argument gives the file providing the setplot function.
    The outdir argument gives the directory providing the fort.* files.
      If outdir==None, the file .output created by 'make .output' is
      examined to determine where the most recent output might be found.
      If .output does not exist, outdir defaults to '.', the current directory.
    Other arguments of Iplotclaw rarely need to be changed:
       completekey='tab', stdin=None, stdout=None
    """
    from pyclaw.plotters import frametools, data

    # initialization:
    # ---------------

    prompt = 'PLOTCLAW > '
    lastcmd = 'n'             # initialize so <return> advances frame'

    def __init__(self, setplot='setplot.py', outdir=None, \
                 completekey='tab', stdin=None, stdout=None):
        """Instantiate a line-oriented interpreter framework.


        The optional argument 'completekey' is the readline name of a
        completion key; it defaults to the Tab key. If completekey is
        not None and the readline module is available, command completion
        is done automatically. The optional arguments stdin and stdout
        specify alternate input and output file objects; if not specified,
        sys.stdin and sys.stdout are used.

        """
        import sys
        if stdin is not None:
            self.stdin = stdin
        else:
            self.stdin = sys.stdin
        if stdout is not None:
            self.stdout = stdout
        else:
            self.stdout = sys.stdout
        self.cmdqueue = []
        self.completekey = completekey
 
        self.setplot = setplot
        plotdata = data.ClawPlotData()
        plotdata.setplot = self.setplot
        plotdata._mode = 'iplotclaw'

        if outdir is None:
            try:
                # possibly set by 'make .output' and stored in .output:
                dot_output = open('.output','r')    
                outdir = dot_output.readline().strip()
                dot_output.close()
            except:
                outdir = '.'
        plotdata.outdir = outdir
        # Note also that outdir might be reset by setplot!

        try:
            plotdata = frametools.call_setplot(self.setplot,plotdata)
        except:
            print '*** Problem executing setplot in Iplotclaw'
            #print '    plotdata = ', plotdata
            print '    setplot = ', self.setplot
            print '*** Either this file does not exist or '
            print '    there is a problem executing the function setplot in this file.'
            print '*** PLOT PARAMETERS MAY NOT BE SET! ***'
            raise
            #return
        self.plotdata = plotdata
        self.prevplotdata = plotdata
        self.restart = False
        self.prevframeno = 0


    def preloop(self):

        print '\nInteractive plotting for Clawpack output... '
        print '\nPlotting data from outdir = ', self.plotdata.outdir
        print 'Type ? at PLOTCLAW prompt for list of commands'

        startframeno = raw_input('\n    Start at which frame [default=%i] ? '\
                                % self.prevframeno)

        makeplot = True 

        if startframeno == '':
            self.frameno = self.prevframeno
            self.plotdata = self.prevplotdata
            if self.restart:
                replot = raw_input('    Replot data for frame %s [no] ? ' \
                                % self.frameno)

                if replot not in ['y','yes','Y']:
                    makeplot = False

                if makeplot:
                    reload = raw_input('    Reload data for frame %s [no] ? ' \
                                    % self.frameno)
                    if reload in ['y','yes','Y']:
                        self.plotdata.refresh_frames = True
                    else:
                        self.plotdata.refresh_frames = False
        else:
            self.plotdata = self.prevplotdata
            try:
                self.frameno = int(startframeno)
            except:
                print '\n    *** Error: frameno must be an integer, resetting to 0'
                self.frameno = 0
        

        if makeplot:
            self.current_data = frametools.plotframe(self.frameno, self.plotdata)
            self.plotdata.refresh_frames=False
            #pylab.draw()

        self.lastcmd = 'n'
        self.restart = True 

    # end of initialization 
    # ---------------------

    def postloop(self):
        #self.framedata = self.plotdata.getframe(self.frameno)
        #grid0 = self.framedata.grids[0]
        #self.q = grid0.q
        #self.mx = grid0.mx
        #if self.ndim > 1:
        #    self.my = grid0.my
        self.prevframeno = self.frameno


    # Commands that can be typed at the PLOTCLAW> prompt:
    # ---------------------------------------------------

    # help command:
    # -------------
    def help_help(self):
        print 'print this list of valid commands\n'

    # next frame:
    # -----------
    def do_n(self, rest):
        from pyclaw.plotters import frametools, data
        #print '    frameno = ',self.frameno
        self.frameno = self.frameno+1
	self.current_data = frametools.plotframe(self.frameno, self.plotdata)
        pylab.draw()
    def help_n(self):
        print 'n: advance to next frame\n'

    # previous frame:
    # ---------------
    def do_p(self, rest):
        #print '    frameno = ',self.frameno
        self.frameno = max(self.frameno-1, 0)
        self.current_data = frametools.plotframe(self.frameno, self.plotdata)
    def help_p(self):
        print 'p: go back to previous frame\n'

    # jump to arbitrary frame:
    # ------------------------
    def do_j(self, rest):
        try:
            newframeno = int(rest)
        except:
            newframeno = raw_input('\n    Jump to which frame? ')
        if newframeno == 'n': 
            self.do_n(rest)
            self.lastcmd = 'n'
        elif newframeno == 'p': 
            self.do_p(rest)
            self.lastcmd = 'p'
        else:
            try:
                newframeno = int(newframeno)
                self.frameno = newframeno
                self.current_data = frametools.plotframe(self.frameno, self.plotdata)
            except:
                print '\n    *** Error: frameno must be an integer, n, or p'
                #print '\n Requested frameno = %s  %s' %(newframeno,type(newframeno))
    def help_j(self):
        print 'j N: jump to frame N\n'
        print 'j:   jump to some other frame (will prompt for N)\n'

    # redraw frame:
    # -------------
    def do_r(self, rest):
        self.current_data = frametools.plotframe(self.frameno, self.plotdata)
    def help_r(self):
        print 'r: redraw the current frame,  rr: reload and redraw\n'

    def do_rr(self, rest):
        #self.plotdata.refresh_frames=True
        outdir = os.path.abspath(self.plotdata.outdir)
        key = (self.frameno, outdir)
        xxx = self.plotdata.framesoln_dict.pop(key,None)
        if xxx is None:
           print 'No frame data to clear for frame ',self.frameno
        else:
           print 'Cleared data for frame ',self.frameno
        print 'Reading data from outdir = ',self.plotdata.outdir
        self.current_data = frametools.plotframe(self.frameno, self.plotdata)
        self.plotdata.refresh_frames=False
    def help_rr(self):
        print 'r: redraw the current frame,  rr: reload and redraw\n'

    # call setplot again
    # --------------------
    def do_resetplot(self, rest):
        if rest:
            self.setplot = rest
            print '*** Resetting setplot to: ',rest
            self.plotdata.setplot = self.setplot
        print 'Executing setplot from ',self.setplot
        try:
            plotdata = frametools.call_setplot(self.setplot,self.plotdata)
        except:
            print '*** Problem re-executing setplot'
            raise

    def help_resetplot(self):
        print 'resetplot: re-execute the function setplot'
        print '           The easiest way to change plotting parameters'
        print '           is to modify setplot.py and then do resetplot.'
        print ' '
        print 'resetplot <new>: switch to a different setplot function'
        print '           as specified by <new>, which is a function or'
        print '           a string specifying the module containing setplot.'

    # show plot parameters:
    # ---------------------
    def do_show(self, rest):
        self.plotdata.showitems()

    def help_show(self):
        print 'show: show the current plot items'

    # clearframes
    # ---------
    def do_clearframes(self, rest):
        if rest=='':
            self.plotdata.framesoln_dict.clear()
            print 'Cleared all frames'
        else:
            try:
                for framestr in rest.split():
                    frameno = int(framestr)
                    outdir = os.path.abspath(self.plotdata.outdir)
                    key = (frameno, outdir)
                    xxx = self.plotdata.framesoln_dict.pop(key,None)
                    if xxx is None:
                       print 'No frame data to clear for frame ',frameno
                    else:
                       print 'Cleared data for frame ',frameno
            except:
                print 'Error in clearframes: unrecognized input'

    def help_clearframes(self):
        print 'clearframes: delete frame data from cache to replot'
        print '    use if you have rerun the code and want to plot the'
        print '    latest results'
        print '          clearframes framenos  clears one or more frames'
        print '          clearframes           clears all frames'

    # cleargauges
    # ---------
    def do_cleargauges(self, rest):
        if rest=='':
            self.plotdata.gaugesoln_dict.clear()
            print 'Cleared all gauges'
        else:
            print 'Not implemented: try cleargauges'

    def help_cleargauges(self):
        print 'cleargauges: delete gauge data from cache to replot'
        print '    use if you have rerun the code and want to plot the'
        print '    latest results'

    # save
    # ---------
    def do_save(self, rest):
        rest = rest.split()
        if len(rest)==2:
            try:
                figno = int(rest[0])
            except:
                print "*** Expected figure number, got: ",rest[0]
            try:
                fname = rest[1]
                pylab.figure(figno)
                pylab.savefig(fname)
                print "Saved figure number %s to file %s" % (figno,fname)
            except:
                print "*** Problem executing savefig"
        else:
            print "*** save requires two arguments: figno, fname"
            print "*** got: ",rest

    def help_save(self):
        print 'save figno fname: save figure figno to file fname using savefig.'


    # print working directory:
    # ------------------------
    def do_pwd(self, rest):
        print '  now in directory: ',os.getcwd()
        print '  data from outdir: ',self.plotdata.outdir
    def help_pwd(self):
        print 'pwd: print current working directory and outdir'
        print '     fort.* files in outdir provide frame data\n'



    # print figure to a file:
    # -----------------------
    def do_print(self, rest):
        #from pyclaw.plotters import frametools
        fname = rest
        try:
            for figno in self.plotdata._fignos:
                if len(fname)>0:
                    # doesn't work properly!
                    frametools.printfig(fname, self.frameno, figno)
                else:
                    frametools.printfig(frameno=self.frameno, figno=figno)
        except:
            print '    *** Error saving figure to file'

    def help_print(self):
        print 'print: print all figures for this frame to files of the form'
        print '      frame000NfigJ.png'
        print 'To print a single figure or with different style, try e.g.'
        print '     PLOTCLAW > q'
        print '     figure(2)'
        print '     savefig("myname.jpg")\n'
        

    # use vi e.g. to edit setplot.py:
    # -------------------------------
    def do_vi(self, rest):
        exitcode = os.system('vi %s' % rest)
        if exitcode != 0:
            print '*** System vi command failed.  Try "help edit"'

    def help_vi(self):
        print 'Edit file using vi, for example to change the plot parameters:'
        print '    PLOTCLAW> vi setplot.py '
        print '    PLOTCLAW> resetplot '
        print 'See also "help edit" for use of other editors.\n'
        

    # edit a file using editor specified by environment variable EDITOR:
    # -----------------------------------------------------------------
    def do_edit(self, rest):
        try:
            editor = os.environ['EDITOR']
            eval("os.system('%s %s')" % (editor,rest))
        except:
            print '*** Environment variable EDITOR not set... '
            print '*** Type "help edit" for more info'

    def help_edit(self):
        print 'Edit file, for example to change the plot parameters:'
        print '    PLOTCLAW> edit setplot.py '
        print '    PLOTCLAW> resetplot '
        print 'Specify the editor by setting environment variable EDITOR'
        print '  before starting Python shell.'
        print 'If you want to use vi, see also "help vi".\n'

        
    # plotgauge commands:
    # --------------
    def do_plotgauge(self, rest):
        gaugesoln_dict = self.plotdata.gaugesoln_dict
        if rest in ['','all']:
            if len(gaugesoln_dict) == 0:
                outdir = os.path.abspath(self.plotdata.outdir)
                try:
                    gauges = self.plotdata.read_gauges(outdir)
                except:
                    print '*** Error reading gauges in Iplotclaw'
                    gauges = {}

                #print '+++ gauges.keys = ',gauges.keys()
                #print '+++ gaugesoln_dict = ',gaugesoln_dict
                try:
                    for (k,v) in gauges.iteritems():
                        gaugesoln_dict[(k, outdir)] = v
                except:
                    raise Exception("*** Problem setting gaugesoln_dict in Iplotclaw")

            if len(gaugesoln_dict) > 0:
                keys = gaugesoln_dict.keys()
                #print '+++ keys = ',keys
                gaugenos = set([keys[k][0] for k in range(len(keys))])
                #print '+++ gaugenos = ',gaugenos
                gaugenos = list(gaugenos)
                gaugenos.sort()
                n = 0
                ans = ""
                while (ans != "q") and (n<len(gaugenos)):
                    gaugeno = gaugenos[n]
                    try:
                        gaugetools.plotgauge(gaugeno, self.plotdata)
                    except:
                        print "*** Problem executing gaugetools.plotgauge with gaugeno = ", gaugeno
                    if n < len(gaugenos)-1:
                        ans = raw_input("      Hit return for next gauge or q to quit ")
                    n = n+1

        else:

            try:
                gaugeno = int(rest)
            except:
                print "Expected gauge number or 'all'"
                gaugeno = None

            try:
                gaugetools.plotgauge(gaugeno, self.plotdata)
            except:
                print "*** Problem executing gaugetools.plotgauge"


    def help_plotgauge(self):
        print 'plotgauge n  : plot figure for gauge number n, if found'
        print 'plotgauge all: loop through plots of all gauges'
        print 'plotgauge    : loop through plots of all gauges'
        
        

    # quit commands:
    # --------------
    def do_quit(self, rest):
        print 'quitting...'
        return True
    def help_quit(self):
        print 'q or quit: terminates the command loop\n'
        
    def do_q(self, rest):
        print 'quitting...'
        return True
    def help_q(self):
        print 'q or quit: terminates the command loop\n'
        
    def do_k(self, rest):
        print 'quitting...'
        return True
    def help_k(self):
        print 'k: terminates the command loop\n'
        
    def do_EOF(self, rest):
        print "quitting..."
        return True
    def help_EOF(self):
        print "Terminates the command loop\n"
        
    # alias plotloop = cmdloop:
    # -------------------------
    def plotloop(self):
        self.cmdloop()


    # Convenience functions for examining solution or making additional
    # plots from the ipython command line:
    # -----------------------------------------------------------------

    def get_frame(self, frameno=None):
        """
        Return the framesoln for Frame for frameno.  
        If frameno is not specified, use the most recently plotted frameno.
        """
        if frameno is None:
            frameno = self.frameno

        return self.plotdata.getframe(frameno)

    def get_t(self, frameno=None):
        """
        Return the time for Frame for frameno.  
        If frameno is not specified, use the most recently plotted frameno.
        """
        if frameno is None:
            frameno = self.frameno

        return self.plotdata.getframe(frameno).t


    def get_grids(self, frameno=None):
        """
        Return the list of grids for frameno.  
        If frameno is not specified, use the most recently plotted frameno.
        """
        if frameno is None:
            frameno = self.frameno

        return self.plotdata.getframe(frameno).grids

    def get_grid(self, frameno=None):
        """
        Return the final grid for frameno.  
        If frameno is not specified, use the most recently plotted frameno.  
        If AMR is not used and there is only one grid, then return this one 
        (rather than a list with one element, as get_grids would return).  
        If AMR is used, then the final grid plotted is returned, 
        similar to claw/matlab/plotclaw behavior where only the final grid 
        is easily available after the plots are made.
        """
        if frameno is None:
            frameno = self.frameno
        
        return self.plotdata.getframe(frameno).grids[-1]

    def otherfigures(self):
        """
        Create any other figures specified in plotdata.otherfigure_dict.
        """
        
        plotdata = self.plotdata
        if len(plotdata.otherfigure_dict)==0:
            print "No other figures specified."
        else:
            for name in plotdata.otherfigure_dict.iterkeys():
                otherfigure = plotdata.otherfigure_dict[name]
                fname = otherfigure.fname
                makefig = otherfigure.makefig
                if makefig:
                    if type(makefig)==str:
                        try:
                            exec(makefig)
                        except:
                            print "*** Problem executing makefig "
                            print "    for otherfigure ",name
                            import pdb; pdb.set_trace()
                    else:
                        try:
                            makefig(plotdata)
                        except:
                            print "*** Problem executing makefig function"
                            print "    for otherfigure ",name
                else:
                    print "No makefig function specified for ",name



# end of Iplotclaw.
#----------------------------------------------------------------------

