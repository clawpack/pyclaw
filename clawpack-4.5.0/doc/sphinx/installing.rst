.. _installing:

**************************************
Installation instructions
**************************************

Prerequisites
-------------

**Operating systems.**
Clawpack should work fine on Unix/Linux or Mac OS X systems.  Much
of it will work under Windows using Cygwin, but this is not officially
supported.


**Fortran.**
The main Clawpack routines are written in Fortran (a mixture of
Fortran 77 and Fortran 90/95) and so compiling and running the code
requires a Fortran compiler, such as `gfortran
<http://gcc.gnu.org/wiki/GFortran>`_.

Makefiles are used in libraries and directories and you will need some
version of *make*.

**Python.**
Starting with Version 4.4, we use Python for visualization of results
(see :ref:`plotting`) and also for user input (see :ref:`setrun`).
Older Matlab plotting scripts are still available but are no longer
being developed and the examples now included in Clawpack include
`setplot.py` files to facilitate use of the Python plotting tools
(see :ref:`setplot`).

You will need Python Version 2.5 or above (but **not** 3.0 or above,
which is not backwards compatible).  You will also need *NumPy* and
*matplotlib* for plotting.  See :ref:`python` for information on
installing the required modules and to get started using Python if
you are not familiar with it.

.. _downloading:

Downloading Clawpack
--------------------

For instructions on using the version of Clawpack in the Subversion repository instead of
the tar file described below, see the `Clawpack wiki
<http://kingkong.amath.washington.edu/trac/clawpack>`_ 

First download the tar file from the Clawpack download page:

  *  `<http://kingkong.amath.washington.edu/clawpack/clawdownload>`_

This file will be of the form clawpack-N.tar.gz  where N is the 
version number.

Move this tar file to the directory where you want to install claw and then::

  $  tar -zxvf clawpack-N.tar.gz
  $  cd clawpack-N


.. _setenv:

Setting environment variables
-----------------------------

In this claw directory modify the *setenv.py* file if necessary and then::

  $  python setenv.py 

This will provide files to set environment variables appropriately.
In particular, the variable `CLAW` should be set to point to this directory.  

Now execute ::

  $ source setenv.bash

if you are using the bash shell, or ::

  $ source setenv.csh

if you use `csh`.  If you don't know what shell you are using, try both and see which one
doesn't give errors, you won't hurt anything.

If you don't know about Unix shells, see these `class notes 
<http://kingkong.amath.washington.edu/uwamath583/sphinx/notes/html/shells.html>`_, for an
introduction and other links.


Consider putting the commands  contained in the appropriate file
`setenv.bash` or `setenv.csh` in your .cshrc or .bashrc
file (which is executed automatically in each new shell you create).   

In particular, the commands found in these files set the following
`environment variables
<http://kingkong.amath.washington.edu/uwamath583/sphinx/notes/html/vars.html>`_

 * `CLAW` is set to the path to the main directory of the Clawpack files.  
 * `PYTHONPATH` is a list of paths that should include $CLAW/python. 
   If this variable is already set in the shell from which you execute `setenv.py`
   then it should provide an extension of the original path to include this.
 * `FC` is set to `gfortran` as the default compiler to use for Fortran.  You may 
   want to change this.

.. _first_test:

Testing your installation and running an example
------------------------------------------------

There are a number of test cases bundled with Clawpack in the directories
`$CLAW/apps` and `$CLAW/book`.  Here and below it is assumed that the
environment variable `CLAW` has been set properly as described above.

As a first test, go to the directory
`$CLAW/apps/advection/1d/example1 <claw/apps/advection/1d/example1>`_.
You can try the following test in this directory, or you may want to first
make a copy of it (see the instructions in :ref:`copyex`).

The Makefiles are set up to do dependency checking so that in many
application directories you can simply type::

  $ make .plots

and the Fortran code will be compiled, data files created, the code
run, and the results plotted automatically, resulting in a set of webpages
showing the results.

However, if this is your first attempt to run a code, it is useful to go
through these steps one at a time, both to understand the steps and so that
any problems with your installation can be properly identified.

You might want to start by examining the Makefile.  This sets a number of
variables, which at some point you might need to modify for other examples,
see :ref:`makefiles` for more about this.  At the bottom of the Makefile is
an `include` statement that points to a common Makefile that is used by most
applications, and where all the details of the make process can be found.

To compile the code, type::

  $ make .exe    

If this gives an error, see :ref:`trouble_makeexe`.

This should compile the example code (after first compiling the required
library routines) and produce an executable named `xclaw` in this directory.

Before running the code, it is necessary to also create a set of data files
that are read in by the Fortran code.  This can be done via::
  
  $ make .data

If this gives an error, see :ref:`trouble_makedata`.

This uses the Python code in `setrun.py` to create data files that have the
form `*.data`.  For the 1d advection example, two files are created,
`claw.data` and `setprob.data`.  The file `claw.data` 
contains standard run-time
parameters of Clawpack (such as the number of grid cells `mx`, indications
of what method to use, what boundary conditions to impose, etc.).  
The file `setprob.data` typically contains parameters specific to a
particular application, in this case the advection velocity `u`.

In Clawpack 4.3 and earlier versions, the user would modify the `claw.data`
and `setprob.data` files directly.  Starting with Clawpack 4.4, the
recommended approach is to only modify the Python function `setrun` defined
in the file `setrun.py`, and use "make .data" to create the `*.data` files.
See :ref:`setrun` for more details.

Once the executable and the data files all exist, we can run the code.  The
recommended way to do this is to type::

  $ make .output

If this gives an error, see :ref:`trouble_makeoutput`.

One could run the code by typing "./xclaw", but using the make option has
several advantages.  For one thing,
this checks dependencies to make sure the executable and data files are up
to date, so you could have typed "make .output" without the first two steps
above.

Also, before running the code a subdirectory `_output` is created
and the output of the code (often a large number of files) is directed to
this subdirectory.  This is convenient if you want to do several runs with
different parameter values and keep the results organized.  After the code
has run you can rename the subdirectory, or you can modify the variable
`OUTDIR` in the Makefile to direct results to a different directory.  See
:ref:`makefiles` for more details.  Copies of all the data files are also
placed in the output directory for future reference.

If the code runs successfully, you should see output like the following::

  Reading data file, first 5 lines are comments: claw.data   
   running...
    
  Reading data file, first 5 lines are comments: setprob.data
  CLAW1EZ: Frame    0 output plot files done at time t =  0.0000D+00
  
  CLAW1... Step   1   Courant number = 5.000  dt =  0.1000D+00  t =  0.1000D+00
  CLAW1 rejecting step... Courant number too large
  CLAW1... Step   1   Courant number = 0.900  dt =  0.1800D-01  t =  0.1800D-01
  CLAW1... Step   2   Courant number = 0.900  dt =  0.1800D-01  t =  0.3600D-01
  CLAW1... Step   3   Courant number = 0.900  dt =  0.1800D-01  t =  0.5400D-01
  CLAW1... Step   4   Courant number = 0.900  dt =  0.1800D-01  t =  0.7200D-01
  CLAW1... Step   5   Courant number = 0.900  dt =  0.1800D-01  t =  0.9000D-01
  CLAW1... Step   6   Courant number = 0.500  dt =  0.1000D-01  t =  0.1000D+00
  CLAW1EZ: Frame    1 output plot files done at time t =  0.1000D+00
  
  --- etc --- etc ---
  
  CLAW1EZ: Frame    9 output plot files done at time t =  0.9000D+00
  
  CLAW1... Step   1   Courant number = 0.900  dt =  0.1800D-01  t =  0.9180D+00
  CLAW1... Step   2   Courant number = 0.900  dt =  0.1800D-01  t =  0.9360D+00
  CLAW1... Step   3   Courant number = 0.900  dt =  0.1800D-01  t =  0.9540D+00
  CLAW1... Step   4   Courant number = 0.900  dt =  0.1800D-01  t =  0.9720D+00
  CLAW1... Step   5   Courant number = 0.900  dt =  0.1800D-01  t =  0.9900D+00
  CLAW1... Step   6   Courant number = 0.500  dt =  0.1000D-01  t =  0.1000D+01
  CLAW1EZ: Frame   10 output plot files done at time t =  0.1000D+01
  
If you don't like seeing output from every time step, you can suppress this by setting
`verbosity = 0` in the file `setrun.py`.  You might try doing that and then typing::

  $ make .output

It should recreate the data files and rerun the code, with less output along the way.

If the code runs properly, the subdirectory `_output` should contain the following files::

    claw.data   fort.q0003  fort.q0008  fort.t0002  fort.t0007
    fort.info   fort.q0004  fort.q0009  fort.t0003  fort.t0008
    fort.q0000  fort.q0005  fort.q0010  fort.t0004  fort.t0009
    fort.q0001  fort.q0006  fort.t0000  fort.t0005  fort.t0010
    fort.q0002  fort.q0007  fort.t0001  fort.t0006  setprob.data

The `fort.info` file contains information about the run just completed.  The files
with names of the form `fort.t000N` and `fort.q000N` contain the computed results for
Frame `N`.  See :ref:`fortfiles` for more information about the contents of these files.

Normally you will not want to examine these files directly, but instead will use a
plotting tool to plot the results.


**Plotting the results.**  
Once the code has run and the files listed above have been created, there are several
options for plotting the results.  

To try the Python tools, type::

  $ make .plots

If this gives an error, see :ref:`trouble_makeplots`.

If this works, it will create a subdirectory named `_plots` that contains a number of
image files (the `*.png` files) and a set of html files that can be used to view the
results from a web browser.  See :ref:`plotting_makeplots` for more details.

An alternative is to view the plots from an interactive Python session, as described in
the section :ref:`plotting_Iplotclaw`.

If you wish to use Matlab instead, see :ref:`matlabplots`.

Other visualization packages could also be used to display the results, but you will need
to figure out how to read in the data.  See :ref:`fortfiles` for information about the
format of the files produced by Clawpack.


**Creating html versions of source files.***

To best view the results, and the source code and README files,
type::

  $ make .htmls

and view the resulting README.html file with a web browser.  

.. _startserver:

Starting a Python web server
-----------------------------

This part is not required, but 
to best view README.html and other Clawpack generated html files,
it is convenient to start a local webserver via::

  $ cd $CLAW
  $ python python/startserver.py

Note that this will take over the window, so do this in a new window, or
else do::

  $ xterm -e python python/startserver.py &

to execute it in a new xterm (if available).
The setenv commands described above will define an alias so that this last
command can be simplified to::

  $ clawserver

The main $CLAW directory will then be available at http://localhost:50005
and jsMath should work properly to display latex on the webpages (once you've
downloaded the required fonts, see
`<http://www.math.union.edu/locate/jsMath/users/fonts.html>`_).  
