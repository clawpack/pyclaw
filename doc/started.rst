.. _installation:

============================
Getting started with PyClaw
============================
This page describes how to set up PyClaw, the serial code.  For the parallel
code, see :ref:`petclaw_start`.

Dependencies
==================
PyClaw and several of its dependencies depend on the availability of a Fortran 95
compiler.  PyClaw is known to work with gfortran on OS X and Linux and the IBM
XLF compiler on the cross-compiled Blue Gene environment.  We hope to support
other Fortran compilers such as Intel, please email the list if you are
interested in helping to port PyClaw to your favorite compiler!

PyClaw relies on the usual Python libraries for scientific computing:

  * `numpy <http://numpy.scipy.org/>`_. Numpy is used both for handling
    arrays in Python and for interfacing between Python and Fortran
    (via f2py).  The current recommended version is 1.6.0.

  * `matplotlib <http://matplotlib.sourceforge.net/>`_.  Matplotlib is
    used for plotting results.  The current recommended version is 1.0.1.

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
PyClaw requires installation of two Clawpack projects: PyClaw itself;
Riemann, a collection of Riemann solvers; and Visclaw, a set of visualization tools
built on top of Matplotlib.  We recommend that you create
a directory to contain all three packages and set your `CLAW` environment 
variable to point to that.  This is accomplished by ::

    $ mkdir clawpack
    $ cd clawpack
    $ export CLAW=`pwd`

in bash.  In csh/tcsh, replace the last command above with ::

    $ setenv CLAW `pwd`

You will probably want to add the command to your `.cshrc` or `.bash_profile`.

The best way to get PyClaw, Riemann, and Visclaw right now is to clone the Git repositories ::

    $ cd $CLAW
    $ git clone  git@github.com:clawpack/pyclaw.git
    $ git clone  git@github.com:clawpack/riemann.git
    $ git clone  git@github.com:clawpack/visclaw.git

Eventually we will have an official release and there will be a tarball available.


Setting up the environment
============================
You will need the following environment variables set:

  * `PYCLAW` must point to the path where you installed PyClaw (usually `$CLAW/pyclaw`)
  * `RIEMANN` must point to the path where you installed Riemann (usually `$CLAW/riemann`)
  * Your `PYTHONPATH` must include PyClaw, Riemann, and VisClaw.

In bash this is accomplished via ::

    $ export RIEMANN=$CLAW/riemann
    $ export PYCLAW=$CLAW/pyclaw
    $ export VISCLAW=$CLAW/visclaw
    $ export PYTHONPATH=$PYCLAW/src:$RIEMANN/:$VISCLAW/src

In csh/tcsh, use ::

    $ setenv RIEMANN $CLAW/riemann
    $ setenv PYCLAW $CLAW/pyclaw
    $ setenv VISCLAW $CLAW/visclaw
    $ setenv PYTHONPATH $PYCLAW/src:$RIEMANN/:$VISCLAW/src

    
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

    $ cd $PYCLAW/apps/advection/1d/constant
    $ make
    $ python advection.py iplot=1

This will run the code and then place you in an interactive plotting shell.
To view the simulation output frames in sequence, simply press 'enter'
repeatedly.  To exit the shell, type 'q'.  For help, type '?' or see
this `Clawpack interactive python plotting help page <http://depts.washington.edu/clawpack/users/plotting.html>`_.

Next steps
================================
Now you're ready to set up your own PyClaw simulation.  Try the :ref:`pyclaw_tutorial`!
