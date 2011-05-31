.. _pyclaw:

******
PyClaw
******

Pyclaw is a Python-based solver for hyperbolic PDEs that includes the algorithms
of Clawpack and SharpClaw.  
It has been designed with easy extensibility, performance, and exploration in mind.
PyClaw is also the basis for PetClaw, a scalable parallel implementation of Clawpack
using PETSc.
At present, PyClaw also includes a python
toolkit containing various other routines to help users use Clawpack.  

You can get the latest development version of PyClaw from
http://github.com/clawpack/pyclaw.

.. toctree::
   :maxdepth: 2

   started
   tutorial
   differences
   develop
   about

*********
PetClaw
*********
PetClaw is a parallel version of Clawpack
The high-level code is written in Python using numpy and based on the 
PyClaw package.
The low-level kernels are written in Fortran and based on Clawpack.
The parallelism is implemented using PETSc and petsc4py.

PetClaw Contents:

.. toctree::
   :maxdepth: 3

   petclaw/started
   petclaw/solvers
   petclaw/plotting
   petclaw/tips
   petclaw/about
   petclaw/installationMac
#   petclaw/rulesProposal



PyClaw Modules reference documentation
======================================
In order to get the most out of Pyclaw, a brief primer into the structure of
the package is in order.  Pyclaw is broken into two main classes and a set of 
functions that operate with those classes.

.. toctree::
   :maxdepth: 3
   
   controller
   data
   evolve/solvers
   evolve/limiters
   io
   solution
   util

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

