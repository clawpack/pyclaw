.. _petclaw_start:

============================
Getting started with PetClaw
============================

Dependencies
==================
First make sure you have a working install of PyClaw.
For PyClaw installation instructions, see :ref:`installation`.

To run PetClaw you'll also need: 

    * `PETSc <http://www.mcs.anl.gov/petsc/petsc-as/>`_  a toolkit for
      parallel scientific computing. The current recommended version is 3.2. 

    * `petsc4py <http://code.google.com/p/petsc4py/>`_: Python bindings for PETSc.
      The current recommended version is the latest petsc4py-dev.

**Obtaining PETSc:**

PETSc 3.2 can be obtained using three approaches. Here we use `Mercurial <http://mercurial.selenic.com/>`_ to get it. Look at `<http://www.mcs.anl.gov/petsc/petsc-as/download/index.html>`_ for more information. 

Do: ::

    $ cd path/to/the/dir/where/you/want/download/Petsc-dev
    $ hg clone http://petsc.cs.iit.edu/petsc/petsc-dev
    $ cd petsc-dev/config
    $ hg clone http://petsc.cs.iit.edu/petsc/BuildSystem BuildSystem

For sh, bash, or zsh shells add the following lines to your shell start-up file: ::
    
    $ export PETSC_DIR=path/to/the/dir/where/you/downloaded/Petsc-dev/petsc-dev
    $ export PETSC_ARCH=your/architecture

whereas for csh/tcsh shells add the following lines to your shell start-up file: ::

    $ setenv PETSC_DIR path/to/the/dir/where/you/downloaded/Petsc-dev/petsc-dev
    $ setenv PETSC_ARCH your/architecture

For more information about PETSC_DIR and PETSC_ARCH, i.e. the variables that 
control the configuration and build process of PETSc, please look at 
`<http://www.mcs.anl.gov/petsc/petsc-as/documentation/installation.html>`_.

Then, if you want PETSc-dev configure for 32-bit use the following command: ::

    $ ./config/configure.py --with-cc=gcc --with-cxx=g++ --with-python=1 --download-mpich=1 --with-shared-libraries=1

whereas, if you want PETSc-dev 64-bit do: ::

    $ ./config/configure.py --with-cc=gcc --with-cxx=g++ --with-python=1 --download-mpich=1 --with-shared-libraries=1 --with-64-bit-indices=1

Note that one of the option is --download-mpich=1. This means that mpich is downloaded. If you do not need/want mpich, remove this option. Note that you need MPI when using PETSc. Therefore, if the option –download-mpich=1 is removed you should have MPI installed on your system or in your user account.

Once the configuration phase is completed, build PETSc libraries with ::

    $ make PETSC_DIR=path/to/the/dir/where/you/have/Petsc-dev PETSC_ARCH=your/architecture all

Check if the libraries are working by running ::

    $ make PETSC_DIR=path/to/the/dir/where/you/have/Petsc-dev PETSC_ARCH=your/architecture test

**Obtaining petsc4py:**

petsc4py is a python binding for PETSc. Since in the previous step PETSc-dev has been installed, we also need to install petsc4py-dev. To install this binding correctly make sure that the PETSC_DIR and PETSC_ARCH are part of your shell start-up file.

Obtain petsc4py-dev with mercurial: ::
    
    $ cd path/to/the/dir/where/you/want/download/petsc4py
    $ hg clone https://petsc4py.googlecode.com/hg/petsc4py -r latest-changeset

The prefered method for the petsc4py iinstallation is `pip <http://pypi.python.org/pypi/pip>`_ ::
    
    $ cd petsc4py-dev
    $ pip install . --user

Indeed, pip removes the old petsc4py installation, downloads the new version of 
`cython <http://cython.org/>`_ (if needed) and installs petsc4py.

To check petsc4py-dev installation do: ::
    
    $ cd petsc4py/test
    $ python runtests.py

All the tests cases should pass, i.e. OK should be printed at the screen.

**NOTE:** An alternative way to install petsc4py is simply using the python 
script setup.py inside petsc4py, i.e. ::
    
    $ cd petsc4py-dev
    $ python setup.py build 
    $ python setup.py install --user


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
the tests passed.  If any fail, please post the output and details of your 
platform on the petclaw-dev Google group.


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

