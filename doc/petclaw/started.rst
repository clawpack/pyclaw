.. _petclaw_start:

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
      The current recommended version is 3.1. We also support PETSc-dev which will be released
      soon as PETSc 3.2. 

    * `petsc4py <http://code.google.com/p/petsc4py/>`_: Python bindings for PETSc.
      The current recommended version is 1.1.2.

For detailed instruction on how to install the PetClaw dependencies on a Mac running 
OS X 10.6.x see :ref:`installationDepsPetClawMacOSX`. 


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

If everything is set up correctly, this will run all the regression tests
(which include pure python code and python/Fortran code) and inform you that
the tests passed.

.. note::

    At the moment, one of the tests requires PETSc-dev.  If you have installed
    a release version of PETSc, this one test will fail.  This test is related
    to an experimental feature (implicit time stepping) and will not affect
    regular functionality of PetClaw.

The testing phase can also be performed only on a sub-set of regression tests
(e.g. pure python code or python and fortran code, classic clawpack or
sharpclaw solver, explicit or implicit time stepping, etc.). This can be
accomplished by passing some attributes to nose. The attributes are already
defined in the regression tests suite and they are:

    * solver_type: classic or sharpclaw
    * kernel_language: python or fortran
    * petsc: True or False
    * time_stepping_mode: explicit or implicit
    * time_stepping_method: ForwardEuler or SSP33 (for the moment)
    * speed: fast or slow

The attribute 'time_stepping_method' is only used in combination with
'solver_type = sharpclaw' because the classic clawpack implements the
Lax-Wendroff scheme.

The attributes can be used in the following ways:

    * Logic AND: run only the regression tests that have the listed attributes ::
    
        $ nosetests -a attribute-1 = value-1,attribute-2 = value-2,attribute-3 = value-3

    * Logic OR: run the regression tests that have at least one of the listed attributes :: 
    
        $ nosetests -a attribute-1 = value-1 -a attribute-2 = value-2 -a attribute-3 = value-3



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
a new problem, please read the :ref:`pyclaw_tutorial`.
The :ref:`pyclaw_reference` and :ref:`riemann_reference` may also be helpful.
