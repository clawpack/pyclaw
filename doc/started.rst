.. _installation:

============================
Getting started with PyClaw
============================

Dependencies
==================
PyClaw relies on the usual Python libraries for scientific computing:

  * `numpy <http://numpy.scipy.org/>`_. Numpy is used both for handling
    arrays in Python and for interfacing between Python and Fortran
    (via f2py).  The current recommended version is 1.6.0.

  * `matplotlib <http://matplotlib.sourceforge.net/>`_.  Matplotlib is
    used for plotting results.  The current recommended version is 1.0.1.

These can be installed via easy_install ::

    $ easy_install numpy
    $ easy_install matplotlib

or pip ::

    $ pip install numpy
    $ pip install matplotlib

or by downloading the 
`Enthought Python Distribution <http://www.enthought.com/products/epd.php>`_.


Installation
==================
PyClaw requires installation of two Clawpack projects: PyClaw itself and
Riemann, a collection of Riemann solvers.  We recommend that you create
a directory to contain both and set your `CLAW` environment variable to point to that.
This is accomplished by ::

    $ mkdir clawpack
    $ cd clawpack
    $ export CLAW=`pwd`

in bash.  In csh/tcsh, replace the last command above with ::

    $ setenv CLAW `pwd`

You will probably want to add the command to your `.cshrc` or `.bash_profile`.

The best way to get PyClaw and Riemann right now is to clone the Git repositories ::

    $ cd $CLAW
    $ git clone  git@github.com:clawpack/pyclaw.git
    $ git clone  git@github.com:clawpack/riemann.git

Eventually we will have an official release and there will be a tarball available.


Setting up the environment
============================
You will need the following environment variables set:

  * `PYCLAW` must point to the path where you installed PyClaw (usually `$CLAW/pyclaw`)
  * `RIEMANN` must point to the path where you installed Riemann (usually `$CLAW/riemann`)

In bash this is accomplished via ::

    $ export RIEMANN=$CLAW/riemann
    $ export PYCLAW=$CLAW/pyclaw

In csh/tcsh, use ::

    $ setenv RIEMANN $CLAW/riemann
    $ setenv PYCLAW $CLAW/pyclaw

Testing your installation
============================
If you don't have it already, we recommend that you install nose ::

    $ easy_install nose

Now simply execute ::

    $ cd $PYCLAW
    $ nosetests

If everything is set up correctly, this will compile the Fortran source,
run several tests, and inform you that the tests passed.  Note that the
tests *must* be run from the main PyClaw directory.

.. note::

    At the moment, nosetests will run both the PyClaw and the PetClaw test suites,
    so the PetClaw tests will fail if you haven't installed PETSc.

Running and plotting an example
================================
Next ::

    $ cd $PYCLAW/apps/advection/1d/constant
    $ make
    $ python advection.py

This will run the code and then place you in an interactive plotting shell.
To view the simulation output frames in sequence, simply press 'enter'
repeatedly.  To exit the shell, type 'q'.  For help, type '?' or see
this `Clawpack interactive python plotting help page <http://kingkong.amath.washington.edu/clawpack/users/plotting.html#interactive-plotting-with-iplotclaw>`_.

Next steps
================================
Now you're ready to set up your own PyClaw simulation.  Try the :ref:`pyclaw_tutorial`!
