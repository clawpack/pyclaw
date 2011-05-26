.. _installationDepsPetClawMacOSX:

=======================================================
Installation of PetClaw dependencies on Mac OS X 10.6.x
=======================================================
This section explains how to install all the PetClaw dependencies on Mac OS X 10.6.x (Snow Leopard).
The following softwares will be installed:

    * `Python <http://www.python.org/>`_.The current recommended python version is 2.7. 
      Python 3.0 will be supported in the near future.
    * `Numpy <http://numpy.scipy.org/>`_ 1.5 or 1.5.1 or 1.6 
    * `PETSc <http://www.mcs.anl.gov/petsc/petsc-as/>`_-dev. (soon PETSc 3.2)     
    * `petsc4py <http://code.google.com/p/petsc4py/>`_-dev
    * `mpi4py <http://mpi4py.scipy.org/docs/usrman/index.html>`_. 
      The current recommended version is 1.2.2.


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
If you install the Mac OS X developer tools available at `<http://developer.apple.com/mac/>`_, the gcc GNU compiler suite will be installed as part of the package. To download the Mac OS X Developer tools an Apple developer connection login is needed. This can be obtained for free at the same address.


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

Then configure PETSc-dev with the following command ::

    $ ./config/configure.py --with-cc=gcc --with-cxx=g++ --with-python=1 --download-mpich=1 --with-shared-libraries=1

Note that one of the option is --download-mpich=1. This means that mpich is downloaded. If you do not need/want mpich, remove this option.

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


Installation of mpi4py
======================
`mpi4py <http://mpi4py.scipy.org/docs/usrman/index.html>`_ is a python bindings for MPI. Therefore, make sure that the MPI distributuion used by PETSc and petsc4py in your system is the same one that is used by mpi4py. During the PETSc configuration the option -–download-mpich=1 has been used. The binaries for the mpich used by PETSc can be found in the path $PETSC_DIR/$PETSC_ARCH/bin. This path should also be added to the environment variable PATH in the shell start-up file, i.e.: 

    * for sh, bash, or zsh shells add the following line to your shell start-up file ::
        
        $ export PATH=$PETSC_DIR/$PETSC_ARCH/bin:$PATH

    * for csh/tcsh shells add the following line to your shell start-up file ::

        $ setenv PATH "$PETSC_DIR/$PETSC_ARCH/bin:$PATH"

Do the aforementioned step before installing mpi4py to guarantee that mpi4py is using the same binaries of mpich. Overlooking this point might cause errors in importing petsc4py.PETSc mpi4py.MPI modules.

Next, add the following line to your shell start-up file::

    $ export ARCHFLAGS="-arch x86_64"

or ::
    
    $ setenv ARCHFLAGS "-arch x86_64"


The current recommended version is 1.2.2. Download it from `<http://code.google.com/p/mpi4py/downloads/list>`_. Afterwards go to the directory where you have got mpi4py-1.2.2.tar.gz and do: ::
    
    $ tar -xzvf (or -xvf) mpi4py-1.2.2.tar.gz
    $ cd mpi4py-1.1.2

Install it: ::

    $ python setup.py install --user

To check mpi4py installation do: ::
    
    $ mpiexec -n 4 python test/runalltest.py
    $ mpiexec -n 4 python demo/helloworld.py

All the tests cases should pass, i.e. you should get
    * OK 
and 
    * Hello, World! I am process 0 of 4 on kl-11638.local. 
    * Hello, World! I am process 1 of 4 on kl-11638.local.
    * Hello, World! I am process 2 of 4 on kl-11638.local.
    * Hello, World! I am process 3 of 4 on kl-11638.local.

for runalltest.py and helloworld.py, respectively.





