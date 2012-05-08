.. _installation:

===============
Installation
===============
This page describes how to set install a serial version of PyClaw.  There are
some additional dependencies for installing a parallel-enabled code in
:ref:`parallel`.  If you encounter any difficulties in the installation
process, please send an email to claw-dev@googlegroups.com.

Dependencies: Python, gfortran, numpy, and matplotlib
--------------------------------------------------------

PyClaw depends on Python 2.5 or greater and a modern Fortran 95
compiler.  PyClaw is known to work with GNU gfortran 4.2 and higher and the IBM
XLF compiler.  In principle, PyClaw should work with other modern Fortran
compilers as well.   PyClaw heavily relies on numpy throughout the code.
matplotlib is used for visualization of results.  You will also need a git
client on your system to obtain PyClaw itself.

  * `Python <http://python.org>`_. PyClaw is written in Python!  The minimum
    recommended version is 2.5.

  * Fortran 95 compiler.  The computational kernels used by PyClaw are written
    in Fortran 95.  You will need to use a Fortran compiler compatible with
    your Python installation.   Most of the PyClaw developers use a gfortran
    from the 4.6.x series.

  * `numpy <http://numpy.scipy.org/>`_ is used both for handling
    arrays in Python and for interfacing between Python and Fortran
    (via f2py).  The minimum recommended version is 1.6.

  * `matplotlib <http://matplotlib.sourceforge.net/>`_ is
    used for plotting results.  The minimum recommended version is 1.0.1.

  * `git <http://git-scm.com/>`_ is the freely available distributed
    version control system used by the PyClaw developers to manage
    development.  The minimum recommended version is 1.7

  * `pip <http://www.pip-installer.org/en/latest/installing.html>`_ is a tool
    for installing and managing Python packages.  We strongly recommend using
    pip to install PyClaw and its dependencies.

Obtaining a freely available Fortran compiler
+++++++++++++++++++++++++++++++++++++++++++++++

Binary files for gfortran on OS X, Linux and Windows are available from
`GCC Wiki GFortranBinaries <http://gcc.gnu.org/wiki/GFortranBinaries>`_.  

Obtaining a git client
+++++++++++++++++++++++++++++++++++++++++++++++

git clients can be downloaded from the `git homepage <http://git-scm.com/download>`_.

Installing pip
+++++++++++++++++++++++++++++++++++++++++++++++

There are instructions for installing pip at `the pip homepage
<http://www.pip-installer.org/en/latest/installing.html>`_.  See `this question on StackOverflow about properly installing virtualenv and pip without administrative privileges <http://stackoverflow.com/questions/4324558/whats-the-proper-way-to-install-pip-virtualenv-and-distribute-for-python>`_ if you would liketo keep your software environment "sandboxed".

Obtaining Python and its dependencies
+++++++++++++++++++++++++++++++++++++++++++++++

The PyClaw developers recommend the `Enthought Python Distribution <http://enthought.com/products/epd.php>`_ to
obtain a modern Python with several important scientific computing libraries
installed "out-of-the-box".   The Enthought Python Distribution is freely
available for academic use.  PyClaw will also try to install the other
dependencies for you, but you will need to install numpy first with:

    pip install numpy

Similarly, you can install nose with:

   pip install nose

Installing PyClaw
-----------------------------------------------------------
PyClaw currently resides within the Clawpack meta-project, which contains:

*PyClaw*
    The actual package containing the PyClaw source code, tests, and examples
    
*Clawutil*
    A package containing important utilities for working with Clawpack projects
    
*Riemann*
    A package containing a collection of Riemann solvers for PyClaw and 
    Clawpack.
    
*VisClaw*
    A set of visualization tools built on top of Matplotlib    

You can install of the packages with:

    pip install clawpack

We recommend installing "in-place" with pip directly from the github repository:

    pip install -e git+git://github.com/clawpack/clawpack.git#egg=clawpack-dev

This will place the install either in a sub-directory named src or your
virtualenv/src directory, depending on if you are using a virtual environment.

You can also perform the install from an existing git clone:

    git clone git://github.com/clawpack/clawpack.git
    cd clawpack
    pip install -e .

Testing your installation with nose
-----------------------------------------------------------

If you have nose, you can test your installation with ::

    cd clawpack/pyclaw
    nosetests 

If everything is set up correctly, this will compile the Fortran source,
run several tests, and inform you that the tests passed.  Note that the
tests *must* be run from the main PyClaw directory.

You can also run the tests in parallel (if you have petsc4py properly installed)
::

    mpirun -n 4 nosetests

.. note::

    PyClaw automatically enables PETSc tests if it detects a proper petsc4py installation.

Running and plotting an example
-----------------------------------------------------------
Next ::

    cd clawpack/pyclaw/apps/advection_1d
    python advection.py iplot=1

This will run the code and then place you in an interactive plotting shell.
To view the simulation output frames in sequence, simply press 'enter'
repeatedly.  To exit the shell, type 'q'.  For help, type '?' or see
this `Clawpack interactive python plotting help page <http://depts.washington.edu/clawpack/users/plotting.html>`_.

Next steps
-----------------------------------------------------------
Now you're ready to set up your own PyClaw simulation.  Try the :ref:`pyclaw_tutorial`!
