.. _pyclaw_geometry:

***************
PyClaw Geometry
***************

The PyClaw geometry package contains the classes used to define the
geometry of a :class:`~pyclaw.solution.Solution` object.  The base container
for all other geometry is the :class:`~pyclaw.geometry.Domain` object.  It
contains a list of :class:`~pyclaw.geometry.Patch` objects that reside inside
of the :class:`~pyclaw.geometry.Domain`.  

.. image:: images/geometry/domain_structure_1.*

:class:`~pyclaw.geometry.Patch` 
represents a piece of the domain that could be a different resolution than
the others, have a different coordinate mapping, or be used to construct
complex domain shapes.  

.. image:: images/geometry/domain_structure_2.*

It contains :class:`~pyclaw.geometry.Dimension`
objects that define the extent of the :class:`~pyclaw.geometry.Patch` and the
number of grid cells in each dimension.  :class:`~pyclaw.geometry.Patch` also
contains a reference to a nearly identical :class:`~pyclaw.geometry.Grid`
object.  The :class:`~pyclaw.geometry.Grid` object also contains a set of
:class:`~pyclaw.geometry.Dimension` objects describing its extent and number
of grid cells.  The :class:`~pyclaw.geometry.Grid` is meant to represent the 
geometry of the data local to the process in the case of a parallel run.  In
a serial simulation the :class:`~pyclaw.geometry.Patch` and
:class:`~pyclaw.geometry.Grid` share the same dimensions.

In the case where only one :class:`~pyclaw.geometry.Patch` object exists in
a :class:`~pyclaw.geometry.Domain` but it is run with four processes in
parallel, the :class:`~pyclaw.geometry.Domain` hierarchy could look like:

.. image:: images/geometry/domain_structure_3.*

In the most complex case with multiple patches and a parallel run we may
have the following:

.. image:: images/geometry/domain_structure_5.*

.. _pyclaw_serial_geometry:

Serial Geometry Objects
#######################


:class:`pyclaw.geometry.Domain`
=================================

.. autoclass:: pyclaw.geometry.Domain
   :members:
   :member-order: groupwise


:class:`pyclaw.geometry.Patch`
=================================

.. autoclass:: pyclaw.geometry.Patch
   :members:
   :member-order: groupwise
   :show-inheritance:


:class:`pyclaw.geometry.Grid`
=================================

.. autoclass:: pyclaw.geometry.Grid
   :members:
   :member-order: groupwise
   

:class:`pyclaw.geometry.Dimension`
==================================

.. autoclass:: pyclaw.geometry.Dimension
   :members:
   :member-order: groupwise

.. _pyclaw_parallel_geometry:

Parallel Geometry Objects
#########################

:class:`petclaw.geometry.Domain`
===================================

.. autoclass:: petclaw.geometry.Domain
   :members:
   :member-order: groupwise
   :show-inheritance:
   
:class:`petclaw.geometry.Patch`
===============================

.. autoclass:: petclaw.geometry.Patch
  :members:
  :member-order: groupwise
  :show-inheritance:
   
