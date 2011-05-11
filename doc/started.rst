.. _installation:

============================
Getting started with PyClaw
============================

Installation
==================
PyClaw relies on the usual Python libraries for scientific computing:

  * numpy

  * matplotlib

These can be installed via easy_install or pip, or by downloading the Enthought
Python Distribution.
PyClaw requires installation of two Clawpack projects: PyClaw itself and
Riemann, a collection of Riemann solvers.  We recommend that you create
a directory to contain both and set your `CLAW` environment variable to point to that.

The best way to get PyClaw and Riemann right now is to clone the Git repositories ::

    $ cd $CLAW
    $ git clone  git@github.com:clawpack/pyclaw.git
    $ git clone  git@github.com:clawpack/riemann.git

Setting up the environment
============================
You will need the following environment variables set:

  * `PYCLAW` must point to the path where you installed PyClaw (usually `$CLAW/pyclaw`)
  * `RIEMANN` must point to the path where you installed Riemann (usually `$CLAW/riemann`)

Testing your installation
============================
If you don't have it already, we recommend that you install nose ::

    $ easy_install nose

Now simply execute ::

    $ cd $PYCLAW
    $ nosetests

If everything is set up correctly, this will compile the Fortran source,
run several examples, and inform you that the tests passed.  Note that the
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
To understand how to set up
a new problem, please read the `PyClaw tutorial <http://kingkong.amath.washington.edu/clawpack/users/pyclaw/tutorial.html>`_.
The `PyClaw reference documentation <http://kingkong.amath.washington.edu/clawpack/users/pyclaw/index.html>`_ may also be helpful.
