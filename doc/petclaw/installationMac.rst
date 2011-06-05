.. _installationDepsPetClawMacOSX:

=======================================================
Installation of PetClaw dependencies on Mac OS X 10.6.x
=======================================================
This section explains how to install all the PetClaw dependencies on Mac OS X 10.6.x (Snow Leopard).
The following software will be installed:

    * `Python <http://www.python.org/>`_.The current recommended python version is 2.7. 
      Python 3.0 will be supported in the near future.
    * `Numpy <http://numpy.scipy.org/>`_ 1.5 or 1.5.1 or 1.6 
    * `PETSc <http://www.mcs.anl.gov/petsc/petsc-as/>`_-dev. (soon PETSc 3.2)     
    * `petsc4py <http://code.google.com/p/petsc4py/>`_-dev

.. note::
   
   PetClaw no longer requires mpi4py.

Installation of python 
======================
Enthought python 2.7 academic distribution is considered in this notes. However, other python distributions also work fine. 

At the time of writing this document the most recent stable release is EPD-7.0-2. Download 32- or 64-bit Enthought python distribution from: `<http://download.enthought.com/academic-epd-7.0/>`_.

Then,

    * Double-click on the .dmg file to mount the image   
    * Double-click on the .mpkg file to run the installer

The installer will install EPD Python 2.7 in your system (typically in /Library/Frameworks/EPD64.framework/) and it will set the path for EPD-7.0-2.


Installation of numpy 
=====================
Enthought Python 2.7 comes already with numpy 1.5.1 (see `EPDChangelog <http://www.enthought.com/EPDChangelog.html>`_) which is one of the supported and recommended version for PetClaw. If you have a different Python distribution check if you have numpy 1.5 or 1.5.1 or 1.6.In case you need to install it, you can use three approaches:

    * Download and install manually the source file: ::
    
        $ cd path/to/the/dir/where/you/want/download/numpy
        $ svn co http://svn.scipy.org/svn/numpy/tags/RELEASE-NUMBER numpy
        $ cd numpy
        $ sudo python setup.py install

    * Use `pip <http://pypi.python.org/pypi/pip>`_: ::

        $ pip install numpy==RELEASE-NUMBER
    

    * Use `easy_install <http://packages.python.org/distribute/easy_install.html>`_ ::
        
        $ easy_install "numpy==RELEASE-NUMBER"

Both method install numpy in the system, i.e. inside the python distribution installed in the previous step. If you prefer to install numpy locally, i.e. only for your user account, append the option ``--user`` after "RELEASE-NUMBER".
 

To test the numpy functionality open a terminal and run python, i.e. ::
   
    $ python

Then type ::

    >>> import numpy
    >>> numpy.test()

You should get something like

    * Ran 2983 tests in 10.194s
    * OK (KNOWNFAIL=4, SKIP=1) <nose.result.TextTestResult run=2983 errors=0 failures=0>


Installation of the gcc (GNU compiler collection)
=================================================
We suggest you to install Xcode. Therefore, sign up at `<https://connect.apple.com/>`_, then download for FREE Xcode from `<http://developer.apple.com/xcode/>`_. 

Although Apple XCode Tools includes gcc 4.X, it is not a complete implementation and lacks `gfortran <http://gcc.gnu.org/wiki/GFortran>`_. However, various implementations of gfortran have been compiled and are available at `GCC Wiki GFortranBinaries <http://gcc.gnu.org/wiki/GFortranBinaries>`_. There are two other possible venues to install gfortran, `Tools - R (and Fortran) for Mac OS X <http://r.research.att.com/tools/>`_, or `High Performance Computing for Mac OS X <http://hpc.sourceforge.net/>`_.  Please take a look at this nice web page `<http://www.webmo.net/support/fortran_osx.html>`_ for more information about fortran and fortran compilers. Here the "GCC Wiki" installation is used. Thus, 
	* Visit the `GCC Wiki GFortranBinaries <http://gcc.gnu.org/wiki/GFortranBinaries>`_ page and scroll to the MacOS section
	* Download the gfortran-macosx-{**architecture you need**} (e.g. Intel processors, Tiger and later; Intel 64-bit processors, Snow Leopard; PowerPC Leopard)
	* Double-click on the .dmg file to mount the image   
    	* Double-click on the .pkg file to run the installer



Installation of PETSc
=====================
In order to use PetClaw in combination with the implicit time stepping schemes implemented in PETSc the development version of PETSc is needed. PETSc-dev can be obtained using three approaches. Here we use `Mercurial <http://mercurial.selenic.com/>`_ to get it. Look at `<http://www.mcs.anl.gov/petsc/petsc-as/developers/index.html>`_ for more information.

Do: ::

    $ cd path/to/the/dir/where/you/want/download/Petsc-dev
    $ hg clone http://petsc.cs.iit.edu/petsc/petsc-dev
    $ cd petsc-dev/config
    $ hg clone http://petsc.cs.iit.edu/petsc/BuildSystem BuildSystem

For sh, bash, or zsh shells add the following lines to your shell start-up file: ::
    
    $ export PETSC_DIR=path/to/the/dir/where/you/downloaded/Petsc-dev/petsc-dev
    $ export PETSC_ARCH=arch-darwin-c-debug

whereas for csh/tcsh shells add the following lines to your shell start-up file: ::

    $ setenv PETSC_DIR path/to/the/dir/where/you/downloaded/Petsc-dev/petsc-dev
    $ setenv PETSC_ARCH arch-darwin-c-debug

Then, if you want PETSc-dev configure for 32-bit use the following command: ::

    $ ./config/configure.py --with-cc=gcc --with-cxx=g++ --with-python=1 --download-mpich=1 --with-shared-libraries=1

whereas, if you want PETSc-dev 64-bit do: ::

    $ ./config/configure.py --with-cc=gcc --with-cxx=g++ --with-python=1 --download-mpich=1 --with-shared-libraries=1 --with-64-bit-indices=1

Note that one of the option is --download-mpich=1. This means that mpich is downloaded. If you do not need/want mpich, remove this option. Note that you need MPI when using PETSc. Therefore, if the option –download-mpich=1 is removed you should have MPI installed on your system or in your user account.

Once the configuration phase is completed, build PETSc libraries with ::

    $ make PETSC_DIR=path/to/the/dir/where/you/have/Petsc-dev PETSC_ARCH=arch-darwin-c-debug all

Check if the libraries are working by running ::

    $ make PETSC_DIR=path/to/the/dir/where/you/have/Petsc-dev PETSC_ARCH=arch-darwin-c-debug test


Installation of petsc4py
========================
`petsc4py <http://code.google.com/p/petsc4py/>`_ is a python binding for PETSc. Since in the previous step PETSc-dev has been installed, we also need to install petsc4py-dev. To install this binding correctly make sure that the PETSC_DIR and PETSC_ARCH are part of your shell start-up file.

Obtain petsc4py-dev with mercurial: ::
    
    $ cd path/to/the/dir/where/you/want/download/petsc4py
    $ hg clone https://petsc4py.googlecode.com/hg/ petsc4py -r latest-changeset

Install it: ::
    
    $ cd petsc4py-dev
    $ python setup.py build --petsc_arch=arch-darwin-c-debug
    $ python setup.py install --user

To check petsc4py-dev installation do: ::
    
    $ cd petsc4py/test
    $ python runtests.py

All the tests cases should pass, i.e. OK should be printed at the screen.



**NOTE:** An alternative easier way to install petsc4py is using again `pip <http://pypi.python.org/pypi/pip>`_., i.e. ::
    
    $ cd petsc4py-dev
    $ pip install . --user
