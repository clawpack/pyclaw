.. _installPetClawDeps:

====================================
Installation of PetClaw dependencies 
====================================
This section explains how to install all the PetClaw dependencies.
The following software will be installed:

    * `PETSc <http://www.mcs.anl.gov/petsc/petsc-as/>`_ 3.2     
    * `petsc4py <http://code.google.com/p/petsc4py/>`_-dev

.. note::
   
   PetClaw no longer requires mpi4py.

.. Installation of python 
.. ======================
.. Enthought python 2.7 academic distribution is considered in this notes. However, other python distributions also work fine. 

.. At the time of writing this document the most recent stable release is EPD-7.0-2. Download 32- or 64-bit Enthought python distribution from: `<http://www.enthought.com/products/edudownload.php>`_.

.. Then,

..    * Double-click on the .dmg file to mount the image   
..    * Double-click on the .mpkg file to run the installer

.. The installer will install EPD Python 2.7 in your system (typically in /Library/Frameworks/EPD64.framework/) and it will set the path for EPD-7.0-2.


.. Installation of the gcc (GNU compiler collection)
.. =================================================
.. We suggest you to install Xcode. Therefore, sign up at `<https://connect.apple.com/>`_, then download for FREE Xcode from `<http://developer.apple.com/xcode/>`_. 

.. Although Apple XCode Tools includes gcc 4.X, it is not a complete implementation and lacks `gfortran <http://gcc.gnu.org/wiki/GFortran>`_. However, various implementations of gfortran have been compiled and are available at `GCC Wiki GFortranBinaries <http://gcc.gnu.org/wiki/GFortranBinaries>`_. There are two other possible venues to install gfortran, `Tools - R (and Fortran) for Mac OS X <http://r.research.att.com/tools/>`_, or `High Performance Computing for Mac OS X <http://hpc.sourceforge.net/>`_.  Please take a look at this nice web page `<http://www.webmo.net/support/fortran_osx.html>`_ for more information about fortran and fortran compilers. Here the "GCC Wiki" installation is used. Thus, 
	* Visit the `GCC Wiki GFortranBinaries <http://gcc.gnu.org/wiki/GFortranBinaries>`_ page and scroll to the MacOS section
	* Download the gfortran-macosx-{**architecture you need**} (e.g. Intel processors, Tiger and later; Intel 64-bit processors, Snow Leopard; PowerPC Leopard)
	* Double-click on the .dmg file to mount the image   
    	* Double-click on the .pkg file to run the installer



Installation of PETSc
=====================
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


Installation of petsc4py
========================
`petsc4py <http://code.google.com/p/petsc4py/>`_ is a python binding for PETSc. Since in the previous step PETSc-dev has been installed, we also need to install petsc4py-dev. To install this binding correctly make sure that the PETSC_DIR and PETSC_ARCH are part of your shell start-up file.

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
