.. _pyclaw:

*******
Welcome
*******

Pyclaw is a `Python <http://www.python.org>`_-based solver for hyperbolic PDEs that includes the algorithms
of `Clawpack <http://www.clawpack.org>`_ and 
`SharpClaw <http://numerics.kaust.edu.sa/sharpclaw/>`_.
It has been designed with easy extensibility, performance, and exploration in mind.
The high-level code is written in Python using numpy and based on the 
PyClaw package.
The low-level kernels used are the Clawpack and SharpClaw routines, which are written in Fortran.
The PyClaw package also includes PetClaw, a scalable parallel implementation of Clawpack and SharpClaw,
using PETSc.  If you are interested in installing and using PetClaw, we recommend that you
verify that you have a working PyClaw installation before following the specific
instructions in :ref:`petclaw_start`.

PyClaw features:

    * Solves general hyperbolic PDEs in 1D and 2D, including mapped grids and surfaces
    * Includes the full functionality of `Clawpack <http://www.clawpack.org>`_ and 
      `SharpClaw <http://numerics.kaust.edu.sa/sharpclaw/>`_
    * Has a simple and intuitive pythonic interface
    * Allows you to run your simulation on the world's biggest supercomputers with 
      the same simple script that runs it on your laptop
    * Makes it easy to access the powerful `PETSc <http://www.mcs.anl.gov/petsc/petsc-as/>`_ 
      library of solvers for handling stiff source terms or for implicit time stepping 
      (under construction)

PyClaw makes use of the additional Clawpack packages, 
`Riemann <http://github.com/clawpack/riemann>`_ and
`VisClaw <http://github.com/clawpack/visclaw>`_ for Riemann solvers and visualization, 
respectively.

Requirements:

    * `Python <http://www.python.org>`_ 2.6 or higher
    * `numpy <http://numpy.scipy.org/>`_ 1.5 or higher
    * Recommended: `matplotlib <http://matplotlib.sourceforge.net/>`_ version 1.0 or higher
      (for plotting solutions)

Additional requirements for PetClaw only: 
    * `PETSc <http://www.mcs.anl.gov/petsc/petsc-as/>`_  3.1 or higher
    * `petsc4py <http://code.google.com/p/petsc4py/>`_ 1.1.2

You can get the latest development version of PyClaw from http://github.com/clawpack/pyclaw.

If you have any issues or need help using PyClaw and PetClaw, `send e-mail <petclaw-dev@googlegroups.com>`_
to the `discussion group <http://groups.google.com/group/petclaw-dev/>`_.

******
PyClaw
******

.. toctree::
   :maxdepth: 1

   started
   tutorial
   plotting
   problem
   parallel
   output
   classes
   troubleshooting
   differences
   develop
   about
   future


.. _pyclaw_reference:

PyClaw Modules reference documentation
======================================
.. toctree::
   :maxdepth: 1
   
   controller
   evolve/solvers
   petclaw/solvers
   evolve/limiters
   io
   solution
   util

.. _riemann_reference:

Riemann Solvers reference documentation
========================================
The Riemann solvers now comprise a separate package.  For convenience,
documentation of the available pure python Riemann solvers is included
here.  Many other Fortran-based Riemann solvers are available.

.. toctree::
   :maxdepth: 3
   
   rp

.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

