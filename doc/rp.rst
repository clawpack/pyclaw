.. _pyclaw_rp:

======================
Riemann Solver Package
======================

This package contains all of the Python-based Riemann solvers.  Each
module solves the Riemann solver for a particular system of hyperbolic
equations.  The solvers all have a common function signature::

    rp_<name>_<dim>d(q_l,q_r,aux_l,aux_r,aux_global)
    
with ``<name>`` replaced with the appropriate solver name and ``<dim>`` with
the appropriate dimension.

:Input:
 - *q_l* - (ndarray(...,num_eqn)) Contains the left states of the Riemann problem
 - *q_r* - (ndarray(...,num_eqn)) Contains the right states of the Riemann problem
 - *aux_l* - (ndarray(...,num_aux)) Contains the left values of the auxiliary array
 - *aux_r* - (ndarray(...,num_aux)) Contains the right values oft he auxiliary array
 - *aux_global* - (dict) Dictionary containing miscellaneous data which is 
    usually problem dependent.
 
:Output:
 - *wave* - (ndarray(...,num_eqn,num_waves)) Contains the resulting waves from the cell
    edge
 - *s* - (ndarray(...,num_waves)) Speeds of each wave
 - *amdq* - (ndarray(...,num_eqn)) Left going fluctuation
 - *apdq* - (ndarray(...,num_eqn)) Right going fluctuation

Except for *aux_global*, all of the input and output values are arrays whose
elements represent grid values with locations indicated by the following scheme
::

    Indexing works like this:  here num_ghost=2 as an example
     0     1     2     3     4     mx+num_ghost-2     mx+num_ghost      mx+num_ghost+2
                 |                        mx+num_ghost-1 |  mx+num_ghost+1
     |     |     |     |     |   ...   |     |     |     |     |
        0     1  |  2     3            mx+num_ghost-2    |mx+num_ghost       
                                              mx+num_ghost-1   mx+num_ghost+1

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

.. note::

    The values ``q_l[i]``, ``q_r[i]`` are the left and right states, respectively, of 
    the ``ith`` Riemann problem.  This convention is different than that used in
    the Fortran Riemann solvers, where ``q_l[i]``, ``q_r[i]`` are the values at the
    left and right edges of a cell.

All of the return values (waves, speeds, and fluctuations) are indexed by cell edge 
(Riemann problem being solved), with ``s[i]`` referring to the wave speed at interface
$i-1/2$.  This follows the same convention used in the Fortran solvers.

See [LeVeque_book_2002]_ for more details.

List of available Riemann solvers:

 * Acoustics_
 * Advection_
 * `Burgers Equation`_
 * `Euler Equations`_
 * `Shallow Water Equations`_

.. _Acoustics:

:mod:`Acoustics <riemann.rp_acoustics>`
================================================

.. automodule:: riemann.rp_acoustics
   :members:
   
.. _Advection:
   
:mod:`Advection <riemann.rp_advection>`
================================================

.. automodule:: riemann.rp_advection
   :members:
   
.. _`Burgers Equation`:

:mod:`Burgers Equation <riemann.rp_burgers>`
=====================================================

.. automodule:: riemann.rp_burgers
   :members:
   
.. _`Euler Equations`:

:mod:`Euler Equations <riemann.rp_euler>`
==================================================

.. automodule:: riemann.rp_euler
   :members:
   
.. _`Shallow Water Equations`:
   
:mod:`Shallow Water Equations <riemann.rp_shallow>`
============================================================

.. automodule:: riemann.rp_shallow
    :member-order: groupwise
    :members:
   
