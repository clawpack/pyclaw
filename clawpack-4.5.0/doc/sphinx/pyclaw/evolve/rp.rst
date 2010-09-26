.. _pyclaw_rp:

======================
Riemann Solver Package
======================

This package contains all of the Riemann solvers provided with pyclaw.  Each
module solves a particular type of Riemann solver listed below with but have
common function signatures that look like::

    rp_<name>_<dim>d(q_l,q_r,aux_l,aux_r,aux_global)
    
with ``<name>`` replaced with the appropriate solver name and ``<dim>`` with
the appropriate dimension.

:Input:
 - *q_l* - (ndarray(...,meqn)) Contains the left states of the Riemann problem
 - *q_r* - (ndarray(...,meqn)) Contains the right states of the Riemann problem
 - *aux_l* - (ndarray(...,maux)) Contains the left values of the auxiliary array
 - *aux_r* - (ndarray(...,maux)) Contains the right values oft he auxiliary array
 - *aux_global* - (dict) Dictionary containing miscellaneous data which is 
    usually problem dependent.
 
:Output:
 - *wave* - (ndarray(...,meqn,mwaves)) Contains the resulting waves from the cell
    edge
 - *s* - (ndarray(...,mwaves)) Speeds of each wave
 - *amdq* - (ndarray(...,meqn)) Left going fluctuation
 - *apdq* - (ndarray(...,meqn)) Right going fluctuation

All of the input and output values are arrays except for *aux_global* which 
are located according to the following scheme
::

    Indexing works like this:  here mbc=2 as an example
     0     1     2     3     4     mx+mbc-2     mx+mbc      mx+mbc+2
                 |                        mx+mbc-1 |  mx+mbc+1
     |     |     |     |     |   ...   |     |     |     |     |
        0     1  |  2     3            mx+mbc-2    |mx+mbc       
                                              mx+mbc-1   mx+mbc+1

    The top indices represent the values that are located on the grid
    cell boundaries such as waves, s and other Riemann problem values, 
    the bottom for the cell centered values.  In particular the ith grid cell 
    boundary has the following related information:
                      i-1         i         i+1
                       |          |          |
                       |   i-1    |     i    |
                       |          |          |
    Again, grid cell boundary quantities are at the top, cell centered
    values are in the cell.

The arrays ``q_l[i]``, ``q_r`` are the left and right state respectively of 
the ``ith`` Riemann problem.  All of the return values are also indexed by
cell edge (Riemann problem being solved).

See [LeVeque_book_2002]_ for more details.

List of available Riemann solvers:

 * Acoustics_
 * Advection_
 * `Burgers Equation`_
 * `Euler Equations`_
 * `Shallow Water Equations`_

.. _Acoustics:

:mod:`Acoustics <pyclaw.evolve.rp.rp_acoustics>`
================================================

.. automodule:: pyclaw.evolve.rp.rp_acoustics
   :members:
   
.. _Advection:
   
:mod:`Advection <pyclaw.evolve.rp.rp_advection>`
================================================

.. automodule:: pyclaw.evolve.rp.rp_advection
   :members:
   
.. _`Burgers Equation`:

:mod:`Burgers Equation <pyclaw.evolve.rp.rp_burgers>`
=====================================================

.. automodule:: pyclaw.evolve.rp.rp_burgers
   :members:
   
.. _`Euler Equations`:

:mod:`Euler Equations <pyclaw.evolve.rp.rp_euler>`
==================================================

.. automodule:: pyclaw.evolve.rp.rp_euler
   :members:
   
.. _`Shallow Water Equations`:
   
:mod:`Shallow Water Equations <pyclaw.evolve.rp.rp_shallow>`
============================================================

.. automodule:: pyclaw.evolve.rp.rp_shallow
    :member-order: groupwise
    :members:
   