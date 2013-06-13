.. _parallel:

============================
Running in parallel
============================
PyClaw can be run in parallel on your desktop or on large supercomputers using the
PETSc library.
Running your PyClaw script in parallel is usually very easy; it mainly consists of
replacing::

    from clawpack import pyclaw

with::
    
    import clawpack.petclaw as pyclaw


Also, most of the provided scripts in `pyclaw/apps` are set up to run in parallel
simply by passing the command-line option `use_petsc=True`.


Installing PETSc
==================
First make sure you have a working install of PyClaw.
For PyClaw installation instructions, see :ref:`installation`.

To run in parallel you'll need: 

    * `PETSc <http://www.mcs.anl.gov/petsc/petsc-as/>`_  a toolkit for
      parallel scientific computing. The current recommended version is 3.3 with latest patch. 

    * `petsc4py <http://code.google.com/p/petsc4py/>`_: Python bindings for PETSc.
      The current recommended version is 3.3.

**Obtaining PETSc:**

PETSc 3.3 can be obtained using three approaches. Here we use `Mercurial <http://mercurial.selenic.com/>`_ to get it. Look at `<http://www.mcs.anl.gov/petsc/petsc-as/download/index.html>`_ for more information. 

Do: ::

    $ cd path/to/the/dir/where/you/want/download/Petsc-3.3
    $ hg clone https://bitbucket.org/petsc/petsc-3.3 petsc-3.3
    $ hg clone https://bitbucket.org/petsc/buildsystem-3.3 petsc-3.3/config/BuildSystem

For sh, bash, or zsh shells add the following lines to your shell start-up file: ::
    
    $ export PETSC_DIR=path/to/the/dir/where/you/downloaded/Petsc-3.3/petsc-3.3
    $ export PETSC_ARCH=your/architecture

whereas for csh/tcsh shells add the following lines to your shell start-up file: ::

    $ setenv PETSC_DIR path/to/the/dir/where/you/downloaded/Petsc-3.3/petsc-3.3
    $ setenv PETSC_ARCH your/architecture

For more information about PETSC_DIR and PETSC_ARCH, i.e. the variables that 
control the configuration and build process of PETSc, please look at 
`<http://www.mcs.anl.gov/petsc/petsc-as/documentation/installation.html>`_.

Then, if you want PETSc-3.3 configure for 32-bit use the following command: ::

    $ ./config/configure.py --with-cc=gcc --with-cxx=g++ --with-python=1 --download-mpich=1 --with-shared-libraries=1

whereas, if you want PETSc-3.3 64-bit do: ::

    $ ./config/configure.py --with-cc=gcc --with-cxx=g++ --with-python=1 --download-mpich=1 --with-shared-libraries=1 --with-64-bit-indices=1

Note that one of the option is --download-mpich=1. This means that mpich is downloaded. If you do not need/want mpich, remove this option. Note that you need MPI when using PETSc. Therefore, if the option –download-mpich=1 is removed you should have MPI installed on your system or in your user account.

Once the configuration phase is completed, build PETSc libraries with ::

    $ make PETSC_DIR=path/to/the/dir/where/you/have/Petsc-3.3/petsc-3.3 PETSC_ARCH=your/architecture all

Check if the libraries are working by running ::

    $ make PETSC_DIR=path/to/the/dir/where/you/have/Petsc-3.3/petsc-3.3 PETSC_ARCH=your/architecture test

**Obtaining petsc4py:**

petsc4py is a python binding for PETSc. We recommend installing petsc4py 3.3 because it is compatible with PETSc 3.3 and 3.2. To install this binding correctly make sure that the PETSC_DIR and PETSC_ARCH are part of your shell start-up file.

Obtain petsc4py-3.3 with mercurial: ::
    
    $ cd path/to/the/dir/where/you/want/download/petsc4py
    $ hg clone https://petsc4py.googlecode.com/hg/petsc4py-3.3 -r 3.3
 
The prefered method for the petsc4py iinstallation is `pip <http://pypi.python.org/pypi/pip>`_ ::
    
    $ cd petsc4py-3.3
    $ pip install . --user

Indeed, pip removes the old petsc4py installation, downloads the new version of 
`cython <http://cython.org/>`_ (if needed) and installs petsc4py.

To check petsc4py-3.3 installation do: ::
    
    $ cd petsc4py-3.3/test
    $ python runtests.py

All the tests cases should pass, i.e. OK should be printed at the screen.

**NOTE:** To run a python code that uses petsc4py in parallel you will need to use mpiexec or mpirun commands. It is important to remember to use the mpiexec or mpirun executables that come with the MPI installation that was used for configuring PETSc installation. If you have used the option --download-mpich=1 while installing PETSc, then the correct mpiexec to use is the one in ${PETSC_DIR}/${PETSC_ARCH}/bin. You can set this mpiexec to be your default by adding this line to your sh, bash, or zsh shell start-up file: ::

    $ export PATH="${PETSC_DIR}/${PETSC_ARCH}/bin:${PATH}"

or this line in case you are using csh or tcsh shells: ::

    $ setenv PATH "${PETSC_DIR}/${PETSC_ARCH}/bin:${PATH}"

You can test that you are using the right mpiexec by running a demonstration script that can be found in $PYCLAW/demo as follows: ::
    
    $ cd $PYCLAW
    $ mpiexec -n 4 python demo/petsc_hello_world.py 

and you should get an output that looks like follows: ::
    
    Hello World! From process 3 out of 4 process(es).
    Hello World! From process 1 out of 4 process(es).
    Hello World! From process 0 out of 4 process(es).
    Hello World! From process 2 out of 4 process(es).
    

**NOTE:** An alternative way to install petsc4py is simply using the python 
script setup.py inside petsc4py, i.e. ::
    
    $ cd petsc4py-dev
    $ python setup.py build 
    $ python setup.py install --user


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
platform on the `claw-users Google group <http://http://groups.google.com/group/claw-users>`_.


Running and plotting an example
================================
Next ::

    $ cd $PYCLAW/apps/advection/1d/constant
    $ python advection.py use_PETSc=True iplot=1

This will run the code and then place you in an interactive plotting shell.
To view the simulation output frames in sequence, simply press 'enter'
repeatedly.  To exit the shell, type 'q'.  For help, type '?' or see
this `Clawpack interactive python plotting help page <http://kingkong.amath.washington.edu/clawpack/users/plotting.html#interactive-plotting-with-iplotclaw>`_.


Tips for making your application run correctly in parallel
================================================================
Generally serial PyClaw code should "just work" in parallel, but if you are not
reasonably careful it is certainly possible to write serial code that will fail
in parallel.

Most importantly, use the appropriate grid attributes.  In serial, both `grid.n` and
`grid.ng` give you the dimensions of the grid (i.e., the number of cells in
each dimension).  In parallel, `grid.n` contains the size
of the whole grid, while `grid.ng` contains just the size of the part that a given
process deals with.  You should typically use only `grid.ng` (you can also use `q.shape[1:]`,
which is equal to `grid.ng`).

Similarly, `grid.lower` contains the lower bounds of the problem domain in the
computational coordinates, whereas `grid.lowerg` contains the lower bounds of the
part of the grid belonging to the current process.  Typically you should use
`grid.lowerg`.

Additionally, be aware that when a Grid object is instantiated in a parallel run,
it is not instantiated as a parallel object.  A typical code excerpt looks like

.. doctest::

    >>> import clawpack.petclaw as pyclaw # doctest: +SKIP
    >>> from clawpack import pyclaw
    >>> mx = 320; my = 80
    >>> x = pyclaw.Dimension('x',0.0,2.0,mx)
    >>> y = pyclaw.Dimension('y',0.0,0.5,my)
    >>> grid = pyclaw.Domain([x,y])

At this point, `grid.ng` is identically equal to `grid.n`, rather than containing
the size of the grid partition on the current process.  Before using it, you
should instantiate a State object

.. doctest::

    >>> num_eqn = 5
    >>> num_aux=1
    >>> state = pyclaw.State(grid,num_eqn,num_aux)

Now `state.grid.ng` contains appropriate information.

Passing options to PETSc
=========================
The built-in applications (see :ref:`apps`) are set up to automatically pass
command-line options starting with a dash ("-") to PETSc.
