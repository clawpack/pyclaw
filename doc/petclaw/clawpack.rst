.. _petclaw_clawpack_solvers:

================================
PetClaw Classic Clawpack Solvers
================================

The PetClaw classic clawpack solvers are a collection of solvers that represent
the functionality of the older versions of clawpack.  It comes in two forms, a
pure Python version and a Python wrapping the fortran libraries.  All of the
solvers available provide the same basic interface and provide the same 
options as the old versions of clawpack.  The superclass solvers are not meant
to be used separately but there to provide common routines for all the
Clawpack solvers.  Please refer to each of the inherited classes for more info
about the methods and attributes they provide each class.  The inheritance
structure is:

.. inheritance-diagram:: petclaw.evolve.clawpack.PetClawSolver1D petclaw.evolve.clawpack.PetClawSolver2D

:Example:
    This is a simple example of how to instantiate and evolve a solution to a
    later time ``t_end`` using the linearized 1d acoustics Riemann solver::
    
        >>> from petclaw.evolve.clawpack import PetClawSolver1D
        >>> solver = PetClawSolver1D()    # Instantiate a default, 1d solver
        >>> solver.mthlim = [3,3]      # Use the van Leer limiter
        >>> solver.dt = 0.0001         # Set the initial time step
        >>> solver.max_steps = 500     # Set the maximum number of time steps
        
        >>> solver.evolve_to_time(solution,t_end)  # Evolve the solution to t_end


:mod:`petclaw.evolve.clawpack`
==============================

.. automodule:: petclaw.evolve.clawpack
   :members:
   :show-inheritance:
