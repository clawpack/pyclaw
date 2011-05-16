============================
Getting started with PetClaw
============================

Installation
==================
The PetClaw package is included in PyClaw.
To run PetClaw you'll need to install 

    * PyClaw and its dependencies

    * PETSc

    * petsc4py

    * mpi4py

For PyClaw installation instructions, see :ref:`installation`.

Setting up the environment
============================
You will need the following environment variables set:

  * `PYCLAW` must point to the path where you installed PyClaw
  * `RIEMANN` must point to the path where you installed Riemann

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
