.. _pyclaw:

******
Pyclaw
******

Pyclaw is a python toolkit containing the basic functionality of Clawpack 
along with various other routines to help users use clawpack.  It has been 
designed with easy extensibility, performance, and exploration in mind.

PyClaw is also the basis for PetClaw, a scalable parallel implementation of Clawpack
using PETSc.

You can get the latest development version of PyClaw from
http://github.com/clawpack/pyclaw.

In order to get the most out of Pyclaw, a brief primer into the structure of
the package is in order.  Pyclaw is broken into two main classes and a set of 
functions that operate with those classes.

.. toctree::
   :maxdepth: 3
   :glob:
   
   controller
   data
   evolve
   io
   plotting
   solution
   tutorial
   util
   petclaw

.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`

