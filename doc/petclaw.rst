.. PetClaw documentation master file, created by
   sphinx-quickstart on Mon Mar 21 10:43:50 2011.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PetClaw: A Parallel, Pythonic, PETSc-fied Clawpack
===================================================

PetClaw is a parallel version of Clawpack
The high-level code is written in Python using numpy and based on the 
PyClaw package.
The low-level kernels are written in Fortran and based on Clawpack.
The parallelism is implemented using PETSc and petsc4py.

The PetClaw user interface is based on the PyClaw package.  
For details, please read the 
`PyClaw tutorial <http://kingkong.amath.washington.edu/clawpack/users/pyclaw/tutorial.html>`_.
The `PyClaw reference documentation <http://kingkong.amath.washington.edu/clawpack/users/pyclaw/index.html>`_ may also be helpful.

Contents:

.. toctree::
   :maxdepth: 2

   petclaw/started
   plotting
   tips
   develop
   about

.. Indices and tables
.. ==================
.. 
.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`


