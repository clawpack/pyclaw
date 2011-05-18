============================
Getting started with PetClaw
============================

Dependencies
==================
First make sure you have a working install of PyClaw.
For PyClaw installation instructions, see :ref:`installation`.

To run PetClaw you'll also need to install 

    * `PETSc <http://www.mcs.anl.gov/petsc/petsc-as/>`_  a toolkit for
      parallel scientific computing.  Installation instructions can be found at
      `<http://www.mcs.anl.gov/petsc/petsc-as/documentation/installation.html>`_.
      The current recommended version is 3.1.

    * `petsc4py <http://code.google.com/p/petsc4py/>`_: Python bindings for PETSc.
      The current recommended version is 1.1.2.

    * `mpi4py <http://mpi4py.scipy.org/docs/usrman/index.html>`_: Python bindings
      for MPI.  This can generally be installed using pip or easy_install as
      detailed here: `<http://mpi4py.scipy.org/docs/usrman/install.html#using-pip-or-easy-install>`_.


Installation
==================
The PetClaw package is included in PyClaw, and no additional installation is required.

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
    $ python advection.py use_PETSc=True

This will run the code and then place you in an interactive plotting shell.
To view the simulation output frames in sequence, simply press 'enter'
repeatedly.  To exit the shell, type 'q'.  For help, type '?' or see
this `Clawpack interactive python plotting help page <http://kingkong.amath.washington.edu/clawpack/users/plotting.html#interactive-plotting-with-iplotclaw>`_.

Next steps
================================
PetClaw is based on the PyClaw package.  To understand how to set up
a new problem, please read the `PyClaw tutorial <http://kingkong.amath.washington.edu/clawpack/users/pyclaw/tutorial.html>`_.
The `PyClaw reference documentation <http://kingkong.amath.washington.edu/clawpack/users/pyclaw/index.html>`_ may also be helpful.
