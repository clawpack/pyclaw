.. _pyclaw_solution:

****************
Pyclaw Solutions
****************

Pyclaw Solutions are containers for Grid and State objects:

.. image:: images/solution.pdf

Each solution contains a list of Grids which in turn contain a list of 
Dimensions, each containing higher level attributes.

List of classes:
 - Solution_
 - State_
 - Grid_
 - Dimension_

.. _Solution:

:class:`pyclaw.solution.Solution`
=================================

.. autoclass:: pyclaw.solution.Solution
   :members:
   :member-order: groupwise

.. _Grid:

:class:`pyclaw.state.State`
=============================
   
.. autoclass:: pyclaw.state.State
   :members:
   :member-order: groupwise

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
