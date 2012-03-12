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

Obtaining a freely available Fortran compiler
+++++++++++++++++++++++++++++++++++++++++++++++

Binary files for gfortran on OS X, Linux and Windows are available from
`GCC Wiki GFortranBinaries <http://gcc.gnu.org/wiki/GFortranBinaries>`_.  

Obtaining a git client
+++++++++++++++++++++++++++++++++++++++++++++++

git clients can be downloaded from the `git homepage <http://git-scm.com/download>`_.

Obtaining Python and its dependencies
+++++++++++++++++++++++++++++++++++++++++++++++

The PyClaw developers recommend the `Enthought Python Distribution <http://enthought.com/products/epd.php>`_ to
obtain a modern Python with several important scientific computing libraries
installed "out-of-the-box".   The Enthought Python Distribution is freely
available for academic use.

Installing PyClaw
-----------------------------------------------------------
PyClaw currently relies on the manual installation of four Clawpack projects:

*PyClaw*
    The actual package containing the PyClaw source code, tests, and examples
    
*Clawutil*
    A package containing important utilities for working with Clawpack projects
    
*Riemann*
    A package containing a collection of Riemann solvers for PyClaw and 
    Clawpack.
    
*VisClaw*
    A set of visualization tools built on top of Matplotlib
    
We recommend that you create a directory to contain all four 
packages.  If you wanted to call this directory ``clawpack``::

    $ mkdir clawpack
    $ cd clawpack

The current method for installing PyClaw, Clawutil, Riemann and VisClaw is to create
a local copy of their github-hosted repositories::

    $ git clone git@github.com:clawpack/clawutil.git
    $ git clone git@github.com:clawpack/pyclaw.git
    $ git clone git@github.com:clawpack/riemann.git
    $ git clone git@github.com:clawpack/visclaw.git

Environment variables
-----------------------------------------------------------
PyClaw currently installs 'in-place'.  That is, .  The following variables are used in PyClaw:

  * ``CLAW`` must point to the base directory you created and cloned the 
    repositories into, above we called this ``clawpack``.
  * ``CLAWUTIL`` must point to the path where you installed Clawutil (usually 
    ``$CLAW/clawutil``) 
  * ``PYCLAW`` must point to the path where you installed PyClaw (usually
    ``$CLAW/pyclaw``) 
  * ``RIEMANN`` must point to the path where you installed Riemann (usually 
    ``$CLAW/riemann``) 
  * Your ``PYTHONPATH`` must include PyClaw, Clawutil, Riemann, and VisClaw.

The easiest way to do this is to use a script provided in Clawutil that 
produces the appropriate environment variables for PyClaw (and the other
Clawpack projects).  To run the script, go into your base directory you 
created above and run ::

    $ python clawutil/src/python/clawutil/setenv.py
    
This script should produce two files that contain the shell script for setting
the above variables.  By default these files are called ``setenv.bash`` and 
``setenv.csh``.  These can be used by running ::

    $ source setenv.bash
    
or ::
    
    $ source setenv.csh
    
depending on your shell (this can be checked by typing ``printenv SHELL`` at
your command line).  The shell code in these files can be copied to your
.bashrc, .cshrc, or .profile file to be run automatically when you open a 
terminal.

Finally, compile the Fortran code for the solvers and Riemann solvers::

    $ cd $PYCLAW/src/pyclaw/clawpack
    $ make
    $ cd $PYCLAW/src/pyclaw/sharpclaw
    $ make
    $ cd $RIEMANN/src/python/riemann
    $ make


Testing your installation with nose
-----------------------------------------------------------

If you have nose, you can test a serial installation with ::

    $ cd $PYCLAW
    $ nosetests -a petsc=False

If everything is set up correctly, this will compile the Fortran source,
run several tests, and inform you that the tests passed.  Note that the
tests *must* be run from the main PyClaw directory.

.. note::

    The flag `-a petsc=False` tells nose not to run the tests that require PETSc.
    If you have installed PETSc and petsc4py, you can run all tests by omitting this
    flag.

Running and plotting an example
-----------------------------------------------------------
Next ::

    $ cd $PYCLAW/apps/advection_1d
    $ python advection.py iplot=1

This will run the code and then place you in an interactive plotting shell.
To view the simulation output frames in sequence, simply press 'enter'
repeatedly.  To exit the shell, type 'q'.  For help, type '?' or see
this `Clawpack interactive python plotting help page <http://depts.washington.edu/clawpack/users/plotting.html>`_.

Next steps
-----------------------------------------------------------
Now you're ready to set up your own PyClaw simulation.  Try the :ref:`pyclaw_tutorial`!
