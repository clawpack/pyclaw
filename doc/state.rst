.. _pyclaw_state:

************
PyClaw State
************

The :class:`~pyclaw.state.State` object records the fields that exist on a given 
:class:`~pyclaw.geometry.Patch`.  These fields include ``q`` and ``aux``.  The
:class:`~pyclaw.state.State` also includes references to the 
:class:`~pyclaw.geometry.Patch` that the state belongs to.

In parallel the :class:`~petclaw.state.State` 
object also handles some of the parallel communication required of the state on the 
given patch such that only the parts of the fields local to the process.  If you
are interested in the geometry of the local state you can find it through the 
:class:`~petclaw.geometry.Patch` object's reference to its own 
:class:`~petclaw.geometry.Grid`.

.. _State:

Serial :class:`pyclaw.state.State`
==================================
   
.. autoclass:: pyclaw.state.State
   :members:
   :member-order: groupwise

Parallel :class:`petclaw.state.State`
===========================================

.. autoclass:: petclaw.state.State
   :members:
   :member-order: groupwise
