============================
Getting started
============================

Installation
==================
To run PetClaw you'll need to install PETSc, Python, 
numpy, and petsc4py, as well as PyClaw.
For some outdated installation instructions, see the 
`old bitbucket wiki <http://bitbucket.org/knepley/wiki/Home>`_.

Setting up the environment
============================
You will need the following environment variables set:

  * `CLAW` must contain the path where you installed clawpack4petclaw
  * `PETCLAW` must contain the path where you installed petclaw

Testing your installation
============================
If you don't have it already, install nose ::

    $ easy_install nose

Now simply execute ::

    $ cd $PYCLAW
    $ nosetests

If everything is set up correctly, this will compile the Fortran source,
run a couple of examples, and inform you that the tests passed.

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
PetClaw is based on the PyClaw package.  To understand how to set up
a new problem, please read the `PyClaw tutorial <http://kingkong.amath.washington.edu/clawpack/users/pyclaw/tutorial.html>`_.
The `PyClaw reference documentation <http://kingkong.amath.washington.edu/clawpack/users/pyclaw/index.html>`_ may also be helpful.
