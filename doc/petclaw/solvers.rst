.. _petclaw_solvers:

==========================
PetClaw Solvers
==========================
The PetClaw solvers function identically to their Pyclaw counterparts (see :ref:`pyclaw_solvers`).
The only difference is that they additionally inherit from the PetSolver 
:class:`petclaw.evolve.solver.PetSolver` class, which adds parallel capability.

.. inheritance-diagram:: petclaw.evolve.clawpack.PetClawSolver1D petclaw.evolve.clawpack.PetClawSolver2D petclaw.evolve.sharpclaw.PetSharpClawSolver1D
   :parts: 2

.. toctree::
   :maxdepth: 3
   :glob:
   
   solver
   clawpack
   sharpclaw

