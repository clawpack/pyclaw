.. _installation:

===============
Getting started 
===============
This page describes how to set up PyClaw, the serial code.  For the parallel
code, see :ref:`parallel`.

Dependencies
==================
PyClaw and several of its dependencies depend on the availability of a Fortran 95
compiler.  PyClaw is known to work with gfortran on OS X and Linux and the IBM
XLF compiler on the cross-compiled Blue Gene environment. The binaries files for
gfortran on OS X, Linux and Windows can be found at 
`GCC Wiki GFortranBinaries <http://gcc.gnu.org/wiki/GFortranBinaries>`_. 
We hope to support other Fortran compilers such as Intel, please email the list
if you are interested in helping to port PyClaw to your favorite compiler!

PyClaw relies on the usual Python libraries for scientific computing:

  * `numpy <http://numpy.scipy.org/>`_. Numpy is used both for handling
    arrays in Python and for interfacing between Python and Fortran
    (via f2py).  The current recommended version is 1.6.0.

  * `matplotlib <http://matplotlib.sourceforge.net/>`_.  Matplotlib is
    used for plotting results.  The current recommended version is 1.0.1.

**Obtaining numpy and matplotlib:**

Some Python distributions come already with numpy 1.5.x or 1.6.x and Matplotlib 
1.0.1 (see for instance `EPDChangelog <http://www.enthought.com/EPDChangelog.html>`_). 
However, in case you need to install it, you can use two different approaches:

    * Use `pip <http://pypi.python.org/pypi/pip>`_: ::

        $ pip install numpy==RELEASE-NUMBER
        $ pip install matplotlib==RELEASE-NUMBER
    

    * Use `easy_install <http://packages.python.org/distribute/easy_install.html>`_ ::
        
        $ easy_install "numpy==RELEASE-NUMBER"
        $ easy_install "matplotlib==RELEASE-NUMBER"

Both methods install numpy in the system. If you prefer to install numpy 
locally, i.e. only for your user account, append the option ``--user`` after 
"RELEASE-NUMBER".
 

To test the numpy functionality open a terminal and run python, i.e. ::
   
    $ python

Then type ::

    >>> import numpy
    >>> numpy.test()

You should get something like

    * Ran 2983 tests in 10.194s
    * OK (KNOWNFAIL=4, SKIP=1) <nose.result.TextTestResult run=2983 errors=0 failures=0>

To test Matplotlib follow the instructions at 
`<http://matplotlib.sourceforge.net/devel/coding_guide.html#testing>`_


Installation
==================
PyClaw requires the installation of four Clawpack projects:

*PyClaw*
    The actual package containing the PyClaw source code, tests, and examples
    
*Clawutil*
    A package containing important utilities for working with Clawpack projects
    
*Riemann*
    A package containing a collection of Riemann solvers for PyClaw and 
    Clawpack.
    
*VisClaw*
    A set of visualization tools built on top of Matplotlib
    
The best way to get PyClaw, Clawutil, Riemann and VisClaw right now is to clone
the Git repositories.  We recommend that you create a directory to contain all three 
packages.  If you wanted to call this directory ``clawpack`` say this would be 
accompilished by ::

    $ mkdir clawpack
    $ cd clawpack

The best way to get PyClaw, Clawutil, Riemann and VisClaw right now is to clone
the Git repositories into the directory you just created::

    $ git clone git@github.com:clawpack/clawutil.git
    $ git clone git@github.com:clawpack/pyclaw.git
    $ git clone git@github.com:clawpack/riemann.git
    $ git clone git@github.com:clawpack/visclaw.git

Eventually we will have an official release and there will be a tarball available.


Setting up the environment
============================
Now we need to setup the environment so we can easily refer to the projects that
we just cloned.  The following variables are used in PyClaw:

  * ``CLAW`` must point to the base directory you created and cloned the 
    repositories into, above we called this ``clawpack``.
  * ``CLAWUTIL`` must point to the path where you installed Clawutil (usually 
    ``$CLAW/clawutil``) 
  * ``PYCLAW`` must point to the path where you installed PyClaw (usually
    ``$CLAW/pyclaw``) 
  * ``RIEMANN`` must point to the path where you installed Riemann (usually 
    ``$CLAW/riemann``) 
  * Your ``PYTHONPATH`` must include PyClaw, Clawutil, Riemann, and VisClaw.

The easiest way to do this is to use a script provided in Clawutil that 
produces the appropriate environment variables for PyClaw (and the other
Clawpack projects).  To run the script, go into your base directory you 
created above and run ::

    $ python clawutil/src/python/setenv.py
    
This script should produce two files that contain the shell script for setting
the above variables.  By default these files are called ``setenv.bash`` and 
``setenv.csh``.  These can be used by running ::

    $ source setenv.bash
    
or ::
    
    $ source setenv.csh
    
depending on your shell (this can be checked by typing ``printenv SHELL`` at
your command line).  The shell code in these files can be copied to your
.bashrc, .cshrc, or .profile file to be run automatically when you open a 
terminal.


Testing your installation
============================
If you don't have it already, we recommend that you install nose ::

    $ easy_install nose

Now simply execute ::

    $ cd $PYCLAW
    $ nosetests -a petsc=False

If everything is set up correctly, this will compile the Fortran source,
run several tests, and inform you that the tests passed.  Note that the
tests *must* be run from the main PyClaw directory.

.. note::

    The flag `-a petsc=False` tells nose not to run the tests that require PETSc.
    If you have installed PETSc and petsc4py, you can run all tests by omitting this
    flag.

Running and plotting an example
================================
Next ::

    $ cd $PYCLAW/apps/advection_1d
    $ make
    $ python advection.py iplot=1

This will run the code and then place you in an interactive plotting shell.
To view the simulation output frames in sequence, simply press 'enter'
repeatedly.  To exit the shell, type 'q'.  For help, type '?' or see
this `Clawpack interactive python plotting help page <http://depts.washington.edu/clawpack/users/plotting.html>`_.

Next steps
================================
Now you're ready to set up your own PyClaw simulation.  Try the :ref:`pyclaw_tutorial`!
