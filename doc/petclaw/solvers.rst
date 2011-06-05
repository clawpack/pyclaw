.. _petclaw_solvers:

==========================
PetClaw Solvers
==========================
The PetClaw solvers function identically to their Pyclaw counterparts (see :ref:`pyclaw_solvers`).
The only difference is that they additionally inherit from the PetSolver 
:class:`petclaw.solver.PetSolver` class, which adds parallel capability.

.. inheritance-diagram:: petclaw.ClawSolver1D petclaw.ClawSolver2D petclaw.SharpClawSolver1D petclaw.SharpClawSolver2D
   :parts: 4

.. toctree::
   :maxdepth: 3
   :glob:
   
   solver
   clawpack
   sharpclaw

