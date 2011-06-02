.. _pyclaw_solution:

****************
Pyclaw Solutions
****************

All Pyclaw solutions work via the following hierarchy:

.. image:: images/pyclaw_solution_structure.*

Each solution contains a list of Grids which in turn contain a list of 
Dimensions, each containing higher level attributes.

List of classes:
 - Solution_
 - Grid_
 - Dimension_

.. _Solution:

:class:`pyclaw.solution.Solution`
=================================

.. autoclass:: pyclaw.solution.Solution
   :members:
   :member-order: groupwise

.. _Grid:

:class:`pyclaw.grid.Grid`
=============================
   
.. autoclass:: pyclaw.grid.Grid
   :members:
   :member-order: groupwise

.. _Dimension:

:class:`pyclaw.grid.Dimension`
==================================

.. autoclass:: pyclaw.grid.Dimension
   :members:
   :member-order: groupwise
