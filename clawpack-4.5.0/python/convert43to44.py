
import os,sys,glob,re,shutil

#=================================================================
def convert():
    if os.path.isfile('claw1ez.data'):
        ndim = 1
        rundata = make_rundata(1)
        make_setplot1(rundata)
        make_Makefile1()
    elif os.path.isfile('claw2ez.data'):
        ndim = 2
        rundata = make_rundata(2)
        make_setplot2(rundata)
        make_Makefile2()
    elif os.path.isfile('amr2ez.data'):
        ndim = 2
        #rundata = make_rundata(2, amr=True)
        print "2d amr not yet implemented"
    elif os.path.isfile('claw3ez.data'):
        ndim = 3
        #rundata = make_rundata(3)
        print "3d not yet implemented"
    elif os.path.isfile('amr3ez.data'):
        ndim = 3
        #rundata = make_rundata(2, amr=True)
        print "3d amr not yet implemented"
    else:
        print "Could not find a clawpack data file"

    make_README(ndim)
    fix_setprob(ndim)

    
#=================================================================
def make_rundata(ndim):
    fname = 'claw%sez.data' % ndim
    fname43 = 'claw%sez.data.claw43' % ndim
    if not os.path.isfile(fname43):
        try:
	    shutil.move(fname, fname43)
	    print "=== Moved %s to %s"  % (fname, fname43)
	except:
	    print "*** Could not find ", fname
            raise
	    return
    clawdata_file = open(fname43,'r')
    lines = clawdata_file.readlines()
    class rundata(object):
        mx          = int(next(lines))  # returns number from next nonempty line
        if ndim>1:
            my          = int(next(lines))
        if ndim>2:
            mz          = int(next(lines))
        nout        = int(next(lines))
        outstyle    = int(next(lines))
	tfinal      = float(next(lines))
        if outstyle!=1:
            outstyle = 1
            tfinal = 1.0
            print "*** Warning: need to check outstyle in setrun.py"
        dt_initial  = float(next(lines))
        dt_max      = float(next(lines))
        cfl_max     = float(next(lines))
        cfl_desired = float(next(lines))
        max_steps   = int(next(lines))
        dt_variable = int(next(lines))
        order       = int(next(lines))
        order_trans = int(next(lines))
        verbosity   = int(next(lines))
        src_split   = int(next(lines))
        mcapa       = int(next(lines))
        maux        = int(next(lines))
        meqn        = int(next(lines))
        mwaves      = int(next(lines))
        mthlim_1    = int(next(lines))
        mthlim      = mwaves*[mthlim_1]   # set all elements to first seen!
        t0          = float(next(lines))
        xlower      = float(next(lines))
        xupper      = float(next(lines))
        if ndim>1:
            ylower      = float(next(lines))
            yupper      = float(next(lines))
        if ndim>2:
            zlower      = float(next(lines))
            zupper      = float(next(lines))
        mbc          = int(next(lines))
        mthbc_xlower = int(next(lines))
        mthbc_xupper = int(next(lines))
        if ndim>1:
            mthbc_ylower = int(next(lines))
            mthbc_yupper = int(next(lines))
        if ndim>2:
            mthbc_zlower = int(next(lines))
            mthbc_zupper = int(next(lines))

    make_setrun(rundata, ndim)
    return rundata

    # end of make_rundata

    
#=================================================================
def next(lines):
    """
    Return a string with the first number from the first nonzero line in
    lines and remove this line (and blanks) from list.
    """
    line = '\n'
    while (line.strip() == '') or (line.strip()[0] == '#'):
        line = lines.pop(0)   # returns first line of lines, and removes it
    pattern = re.compile(r"^[ ]*(?P<string>[\S]+?)[, ]")
    result = pattern.search(line)
    try:
        numstring = result.group('string')
        numstring = numstring.replace('d','e')  # for floating notation
    except:
        print '*** Failed on line: ',line
	numstring = '0'
    return numstring


#========================================================================

def make_setrun(d,ndim):
    """
    Create a file setrun.py using the run-time parameters in d.
    """
    setrun = open('setrun.py', 'w')
    setrun.write('""" ')
    setrun.write("""
Module to set up run time parameters for Clawpack.

The values set in the function setrun are then written out to data files
that will be read in by the Fortran code.
    """)
    setrun.write('\n""" ')
    setrun.write("""

import os
from pyclaw import data 


#------------------------------
def setrun(claw_pkg='classic'):
#------------------------------
    """)
    setrun.write('\n    """ ')
    setrun.write("""
    Define the parameters used for running Clawpack.

    INPUT:
        claw_pkg expected to be "classic" for this setrun.

    OUTPUT:
        rundata - object of class ClawRunData 
    """)
    setrun.write('\n    """ ')
    setrun.write("""
    
    assert claw_pkg.lower() == 'classic',  "Expected claw_pkg = 'classic'"

    ndim = %s
    rundata = data.ClawRunData(claw_pkg, ndim)

    #------------------------------------------------------------------
    # Problem-specific parameters to be written to setprob.data:
    #------------------------------------------------------------------

    # MUST BE ADDED BY USER IF DESIRED !!
    #probdata = rundata.new_UserData(name='probdata',fname='setprob.data')
    # For example:
    #probdata.add_param('u',     1.0,  'advection velocity')
    
    #------------------------------------------------------------------
    # Standard Clawpack parameters to be written to claw.data:
    #------------------------------------------------------------------

    clawdata = rundata.clawdata  # initialized when rundata instantiated

    # ---------------
    # Spatial domain:
    # ---------------

    # Number of space dimensions:
    clawdata.ndim = ndim
    
    # Lower and upper edge of computational domain:
    clawdata.xlower = %s
    clawdata.xupper = %s
    """ % (ndim, d.xlower, d.xupper))

    if ndim>1:
        setrun.write("""
    clawdata.ylower = %s
    clawdata.yupper = %s
        """ % (d.ylower, d.yupper))
    if ndim>2:
        setrun.write("""
    clawdata.zlower = %s
    clawdata.zupper = %s
        """ % (d.zlower, d.zupper))

    setrun.write("""

    # Number of grid cells:
    clawdata.mx = %s
    """ % d.mx)

    if ndim>1:
        setrun.write("""
    clawdata.my = %s
        """ % d.my)
    if ndim>2:
        setrun.write("""
    clawdata.mz = %s
        """ % d.mz)

    setrun.write("""

    # ---------------
    # Size of system:
    # ---------------

    # Number of equations in the system:
    clawdata.meqn = %s

    # Number of auxiliary variables in the aux array (initialized in setaux)
    clawdata.maux = %s
    
    # Index of aux array corresponding to capacity function, if there is one:
    clawdata.mcapa = %s
    
    
    
    # -------------
    # Initial time:
    # -------------

    clawdata.t0 = %s
    
    """ % (d.meqn, d.maux, d.mcapa, d.t0))



    setrun.write("""
    # -------------
    # Output times:
    #--------------

    # Specify at what times the results should be written to fort.q files.
    # Note that the time integration stops after the final output time.
    # The solution at initial time t0 is always written in addition.

    clawdata.outstyle = %s

    if clawdata.outstyle==1:
        # Output nout frames at equally spaced times up to tfinal:
        clawdata.nout = %s
        clawdata.tfinal = %s

    elif clawdata.outstyle == 2:
        # Specify a list of output times.  
        clawdata.tout =  [0.5, 1.0]   # used if outstyle == 2
        clawdata.nout = len(clawdata.tout)

    elif clawdata.outstyle == 3:
        # Output every iout timesteps with a total of ntot time steps:
        iout = 1
        ntot = 5
        clawdata.iout = [iout, ntot]
    """ % (d.outstyle, d.nout, d.tfinal))
    setrun.write("""


    # ---------------------------------------------------
    # Verbosity of messages to screen during integration:  
    # ---------------------------------------------------

    # The current t, dt, and cfl will be printed every time step
    # at AMR levels <= verbosity.  Set verbosity = 0 for no printing.
    #   (E.g. verbosity == 2 means print only on levels 1 and 2.)
    clawdata.verbosity = %s
    
    

    # --------------
    # Time stepping:
    # --------------

    # if dt_variable==1: variable time steps used based on cfl_desired,
    # if dt_variable==0: fixed time steps dt = dt_initial will always be used.
    clawdata.dt_variable = %s
    
    # Initial time step for variable dt.  
    # If dt_variable==0 then dt=dt_initial for all steps:
    clawdata.dt_initial = %s
    
    # Max time step to be allowed if variable dt used:
    clawdata.dt_max = %s
    
    # Desired Courant number if variable dt used, and max to allow without 
    # retaking step with a smaller dt:
    clawdata.cfl_desired = %s
    clawdata.cfl_max = %s
    
    # Maximum number of time steps to allow between output times:
    clawdata.max_steps = %s

    """ % (d.verbosity, d.dt_variable, d.dt_initial, d.dt_max, \
           d.cfl_desired, d.cfl_max, d.max_steps))

    setrun.write("""
    

    # ------------------
    # Method to be used:
    # ------------------

    # Order of accuracy:  1 => Godunov,  2 => Lax-Wendroff plus limiters
    clawdata.order = %s
    
    # Transverse order for 2d or 3d (not used in 1d):
    clawdata.order_trans = %s
    
    # Number of waves in the Riemann solution:
    clawdata.mwaves = %s
    
    # List of limiters to use for each wave family:  
    # Required:  len(mthlim) == mwaves
    clawdata.mthlim = %s
    
    # Source terms splitting:
    #   src_split == 0  => no source term (src routine never called)
    #   src_split == 1  => Godunov (1st order) splitting used, 
    #   src_split == 2  => Strang (2nd order) splitting used,  not recommended.
    clawdata.src_split = %s
    
    """ % (d.order, d.order_trans, d.mwaves, d.mthlim, d.src_split))
    
    setrun.write("""
    # --------------------
    # Boundary conditions:
    # --------------------

    # Number of ghost cells (usually 2)
    clawdata.mbc = %s
    
    # Choice of BCs at xlower and xupper:
    #   0 => user specified (must modify bcN.f to use this option)
    #   1 => extrapolation (non-reflecting outflow)
    #   2 => periodic (must specify this at both boundaries)
    #   3 => solid wall for systems where q(2) is normal velocity
    
    clawdata.mthbc_xlower = %s
    clawdata.mthbc_xupper = %s
    """ % (d.mbc, d.mthbc_xlower, d.mthbc_xupper))

    if ndim>1:
        setrun.write("""
    clawdata.mthbc_ylower = %s
    clawdata.mthbc_yupper = %s
    """ % (d.mthbc_ylower, d.mthbc_yupper))
    if ndim>2:
        setrun.write("""
    clawdata.mthbc_zlower = %s
    clawdata.mthbc_zupper = %s
    """ % (d.mthbc_zlower, d.mthbc_zupper))

    setrun.write("""
    return rundata
    # end of function setrun
    # ----------------------


if __name__ == '__main__':
    # Set up run-time parameters and write all data files.
    import sys
    if len(sys.argv) == 2:
	rundata = setrun(sys.argv[1])
    else:
	rundata = setrun()

    rundata.write()
    """)


    setrun.close()
    print "=== Created setrun.py"
    # end of make_setrun

#========================================================================

def make_setplot1(d):
    """
    Create a file setplot.py using the plotting parameters in d.
    """
    setplot = open('setplot.py', 'w')
    setplot.write('\n""" ')
    setplot.write("""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    """)
    setplot.write('\n""" ')
    setplot.write("""

#--------------------------
def setplot(plotdata):
#--------------------------
    """)
    setplot.write('\n    """ ')
    setplot.write("""
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    """)
    setplot.write('\n    """ ')
    setplot.write("""

    plotdata.clearfigures()  # clear any old figures,axes,items data

    """)

    # create a figure for each component of q:

    for iq in range(d.meqn):
        setplot.write("""

    # Figure for q[%s]
    plotfigure = plotdata.new_plotfigure(name='q[%s]', figno=%s)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q[%s]'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = %s
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    """ % (iq,iq,iq,iq,iq))


    setplot.write("""

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    """)
    setplot.close()
    print "=== Created setplot.py"
    # end of make_setplot1

#========================================================================

def make_setplot2(d):
    """
    Create a file setplot.py using the plotting parameters in d.
    """
    setplot = open('setplot.py', 'w')
    setplot.write('\n""" ')
    setplot.write("""
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    """)
    setplot.write('\n""" ')
    setplot.write("""

#--------------------------
def setplot(plotdata):
#--------------------------
    """)
    setplot.write('\n    """ ')
    setplot.write("""
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    """)
    setplot.write('\n    """ ')
    setplot.write("""


    from pyclaw.plotters import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data
    """)

    # create a figure for each component of q:

    for iq in range(d.meqn):
        setplot.write("""

    # Figure for q[%s]
    plotfigure = plotdata.new_plotfigure(name='q[%s]', figno=%s)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'q[%s]'
    plotaxes.afteraxes = "pylab.axis('scaled')" 

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = %s
    plotitem.pcolor_cmap = colormaps.yellow_red_blue
    plotitem.add_colorbar = True
    plotitem.show = True       # show on plot?
    """ % (iq,iq,iq,iq,iq))


    setplot.write("""

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

    """)
    setplot.close()
    print "=== Created setplot.py"

    # end of make_setplot2

# =================================================================
def make_Makefile1():
    if not os.path.isfile('Makefile.claw43'):
        try:
	    shutil.move('Makefile','Makefile.claw43')
	    print "=== Moved Makefile to Makefile.claw43"
	except:
	    print "*** Could not find Makefile"
            raise
	    return
    oldmake = open('Makefile.claw43','r')
    mkfile = oldmake.read()

    newmake = open('Makefile','w')
    newmake.write("""
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/util/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = Classic                  # Clawpack package to use
CLAW_EXE = xclaw                    # Executable to create
CLAW_setrun_file = setrun.py        # File containing function to make data
CLAW_OUTDIR = _output               # Directory for output
CLAW_setplot_file = setplot.py      # File containing function to set plots
CLAW_PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# Add any desired compiler flags such as -g here:
FFLAGS =


# ---------------------------------
# List of sources for this program:
# ---------------------------------

    """)

    pattern = re.compile(r"OBJECTS =(?P<objs>.*?)LIBOBJECTS", re.DOTALL)
    result = pattern.search(mkfile)
    if result:
        objs = result.group('objs')  # local object files
    else:
        print "*** No local files?"

    newmake.write("\nCLAW_SOURCES =")
    objs = objs.replace('.o','.f',100)
    newmake.write(objs)

    newmake.write("""
# Clawpack library to be used:
CLAW_LIB = $(CLAW)/clawpack/1d/lib
    """)

    pattern = re.compile(r"LIBOBJECTS =(?P<libobjs>.*?)[\s]*SOURCES", re.DOTALL)
    result = pattern.search(mkfile)
    if result:
        libobjs = result.group('libobjs')  # library object files
    else:
        print "*** No library files?"

    libobjs = libobjs.replace('.o','.f',100)
    libobjs = libobjs.replace(r'$(CLAW)/clawpack/1d/lib','$(CLAW_LIB)',100)
    if 'out1' not in libobjs:
        # Sometimes Makefile has a switch between out1 and out1_hdf
        # Here just add out1.f:
        libobjs = libobjs + "\\\n  $(CLAW_LIB)/out1.f"
    libobjs = libobjs + "\\\n  $(CLAW_LIB)/opendatafile.f"
    newmake.write("\nCLAW_LIBSOURCES =")
    newmake.write(libobjs)

    newmake.write("""

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


### DO NOT remove this line - make depends on it ###
    """)

    newmake.close()
    print "=== Created Makefile"

    ## end of make_Makefile1



# =================================================================
def make_Makefile2():
    if not os.path.isfile('Makefile.claw43'):
        try:
	    shutil.move('Makefile','Makefile.claw43')
	    print "=== Moved Makefile to Makefile.claw43"
	except:
	    print "*** Could not find Makefile"
            raise
	    return
    oldmake = open('Makefile.claw43','r')
    mkfile = oldmake.read()

    newmake = open('Makefile','w')
    newmake.write("""
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/util/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = Classic                  # Clawpack package to use
CLAW_EXE = xclaw                    # Executable to create
CLAW_setrun_file = setrun.py        # File containing function to make data
CLAW_OUTDIR = _output               # Directory for output
CLAW_setplot_file = setplot.py      # File containing function to set plots
CLAW_PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# Add any desired compiler flags such as -g here:
FFLAGS =


# ---------------------------------
# List of sources for this program:
# ---------------------------------

    """)

    pattern = re.compile(r"OBJECTS =(?P<objs>.*?)LIBOBJECTS", re.DOTALL)
    result = pattern.search(mkfile)
    if result:
        objs = result.group('objs')  # local object files
    else:
        print "*** No local files?"

    newmake.write("\nCLAW_SOURCES =")
    objs = objs.replace('.o','.f',100)
    newmake.write(objs)

    newmake.write("""
# Clawpack library to be used:
CLAW_LIB = $(CLAW)/clawpack/2d/lib
    """)

    pattern = re.compile(r"LIBOBJECTS =(?P<libobjs>.*?)[\s]*SOURCES", re.DOTALL)
    result = pattern.search(mkfile)
    if result:
        libobjs = result.group('libobjs')  # library object files
    else:
        print "*** No library files?"

    libobjs = libobjs.replace('.o','.f',100)
    libobjs = libobjs.replace(r'$(CLAW)/clawpack/2d/lib','$(CLAW_LIB)',100)
    if 'out2' not in libobjs:
        # Sometimes Makefile has a switch between out2 and out2_hdf
        # Here just add out2.f:
        libobjs = libobjs + "\\\n  $(CLAW_LIB)/out2.f"
    if 'restart2' not in libobjs:
        # Sometimes Makefile has a switch between restart2 and restart2_hdf
        # Here just add restart2.f:
        libobjs = libobjs + "\\\n  $(CLAW_LIB)/restart2.f"
    libobjs = libobjs + "\\\n  $(CLAW_LIB)/opendatafile.f"
    newmake.write("\nCLAW_LIBSOURCES =")
    newmake.write(libobjs)

    newmake.write("""

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


### DO NOT remove this line - make depends on it ###
    """)

    newmake.close()
    print "=== Created Makefile"

    ## end of make_Makefile1




# =================================================================

def make_README(ndim):



    readme = open('README.txt','w')
    readme.write("""

begin_html [use: jsMath] [use: doc/doc.css]
<!--   For a more readable version of this file, execute
                  unix>  make htmls
       in this directory and then point your browser to README.html
     --------------------------------------------------------------  -->

<h2>
CLAWPACK Sample Code
</h2>

Add a description here!
    """)

    regexp = re.compile(r'.*book/chap(?P<ex>.*)')
    result = regexp.search(os.getcwd())
    if result:
        readme.write("""
Example [book/chap%s]
to accompany the book <br> &nbsp;&nbsp;
  [www.clawpack.org/book.html Finite Volume Methods for Hyperbolic Problems]
  by R. J. LeVeque.

Converted to [www.clawpack.org Clawpack 4.4] form in 2009.
        """ % result.group('ex'))


    readme.write("""
<h4>
Plots of results
</h4>
After running this code and creating plots via "make .plots", you should be
able to view the plots in [link: _plots/_PlotIndex.html].


<h4>
Fortran files
</h4>

<dl>
<dt>[code: Makefile]
<dd> Determines which version of fortran files
are used when compiling the code with make and specifies where output and
plots should be directed.  Type "make .help" at the Unix prompt for options.

<dt>[code: driver.f]
<dd>
The driver routine allocates storage and then calls the main Clawpack
routine.

<dt>[code: setprob.f]
<dd>
A routine by this name is called by the library routine
[clawcode: clawpack/%sd/lib/claw%sez.f]
and is generally used to set any values needed for the specific problem
being solved.
    """ % (ndim,ndim))

    if ndim==1:
        readme.write("""
    
<dt>[code: rp1.f]
<dd>
This is the Riemann solver, which takes the $q$ values stored in the
arrays <tt>ql</tt> and <tt>qr</tt> and returns the waves in the array
<tt>wave</tt> and speeds in the array <tt>s</tt> that result in solving the
Riemann problem at each cell interface, and the fluctuations <tt>amdq</tt>
and <tt>apdq</tt>.  See [claw:doc/rp1.html] for more information about 1d
Riemann solvers.
         """)
    if ndim==2:
        readme.write("""
    
<dt>[code: rpn2.f]
<dd> Normal Riemann solver (normal to cell interface).

<dt>[code: rpt2.f]
<dd> Transverse Riemann solver.

         """)

    readme.write("""
<dt>[code: qinit.f]
<dd>
This subroutine sets the initial data at time $t=0$.

</dl>

<h4>
Python files
</h4>
<dl>

<dt>[code: setrun.py]
<dd> This file contains a function that
specifies what run-time parameters will be used.

Some parameters that you might want to modify are described in the
[www.clawpack.org/doc.html documentation].

<dt>[code: setplot.py]
<dd> This file contains a function that
specifies what plots will be done and
sets various plotting parameters.

</dl>


<h4>
Data files
</h4>
<font color='red'>Warning:</font> These files are generally changed
when setting up a run, usually in [code: setrun.py].

<dl>
<dt>[code: claw.data]
<dd> This file contains a number of
parameter values that are used by CLAWPACK.
The values in this file are read by the library routine
[clawcode: clawpack/%sd/lib/claw%sez.f].


<dt> [code: setprob.data]
<dd> This file may contain various
parameters used in setting the initial conditions or otherwise setting up
the problem.


</dl>

<h4>Library routines</h3>

In addition to the Fortran routines in this library, several library
routines from [claw:clawpack/%sd/lib] are used.  See the [code: Makefile]
to determine which ones are used.

end_html

    """ % (ndim,ndim,ndim))
    readme.close()
    print "=== Created README.txt"

    ## end of make_README


def fix_setprob(ndim):
    """
    First pass at modifying setprob to use the opendatafile subroutine in
    place of open.
    """
    import re

    if (os.path.isfile('setprob.f') and \
          (not os.path.isfile('setprob.f.claw43'))):
	shutil.move('setprob.f','setprob.f.claw43')
	print "=== Moved setprob.f to setprob.f.claw43"
        lines = open('setprob.f.claw43','r').readlines()
        setprob = open('setprob.f','w')
	for line in lines:
	    if line.find('implicit') > -1:
	        setprob.write(line)
		setprob.write("      character*12 fname\n")
	    elif line.find('open(') > -1:
	        regexp = re.compile(r"open.*unit\s*?=\s*?(?P<iunit>[0-9]+)\s*,")
		result = regexp.search(line)
		if result:
		    iunit = result.group('iunit')
		else:
		    print '*** Oops, expected to find unit number in setprob.f'
		    print '*** setprob.f is corrupted, revert from setprob.f.claw43'
		    setprob.close()
                    raise
		    return
		setprob.write("c\n      iunit = %s" % iunit)
		setprob.write("""
      fname = 'setprob.data'
c     # open the unit with new routine from Clawpack 4.4 to skip over
c     # comment lines starting with #:
      call opendatafile(iunit, fname)
                \n""")

            else:
	        setprob.write(line)
	setprob.close()
  
	print "=== Modified setprob.f"

    ## end of fix_setprob
   


# =================================================================
    
if __name__ == "__main__":
    try:
        convert()
        print "Done converting."
    except:
        print "Problem converting"
