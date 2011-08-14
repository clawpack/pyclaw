r"""
Module containg the PetClaw solvers

This module contains the pure and wrapped PetClaw solvers.  All 
PetClaw solvers inherit from the :class:`PetClawSolver` superclass which in turn 
inherits from the :class:`~petclaw.solver.PetSolver` superclass.  As such, 
the only solver classes that should be directly used should be the 
dimensionally dependent ones such as :class:`petclaw.clawpack.ClawSolver1D`.

:Authors:
    Amal Alghamdi
    David Ketcheson
"""
from petclaw.solver import PetSolver
import pyclaw.clawpack

# ============================================================================
#  PetClaw 1d Solver Class
# ============================================================================
class ClawSolver1D(PetSolver,pyclaw.clawpack.ClawSolver1D):
    r"""
    PetClaw solver for 1D problems using classic Clawpack algorithms.

    This class implements nothing; it just inherits from PetClawSolver and
    ClawSolver1D.
    """
    def setup(self,solution):
        r"""
        Set Fortran data structures (for Clawpack) and set up a DA with
        the appropriate stencil width.
        """
        # This is a hack to deal with the fact that petsc4py
        # doesn't allow us to change the stencil_width (mbc)
        state = solution.state
        state.set_stencil_width(self.mbc)
        # End hack

        self.set_mthlim()
        
        if(self.kernel_language == 'Fortran'):
            self.set_fortran_parameters(solution)

# ============================================================================
#  PetClaw 2d Solver Class
# ============================================================================
class ClawSolver2D(PetSolver,pyclaw.clawpack.ClawSolver2D):
    r"""
    PetClaw solver for 2D problems using classic Clawpack algorithms.

    This class implements nothing; it just inherits from PetClawSolver and
    ClawSolver2D.
    
    Note that only the fortran routines are supported for now in 2D.
    """
    def setup(self,solution):
        r"""
        See setup doc string in the super class.
        We are initializing (allocating) the working arrays needed by fortran kernels 
        in this routine. These arrays are passed in each call to the fortran kernel dimsp2.
        """
        # This is a hack to deal with the fact that petsc4py
        # doesn't allow us to change the stencil_width (mbc)
        state = solution.state
        state.set_stencil_width(self.mbc)
        # End hack

        self.set_mthlim()

        if(self.kernel_language == 'Fortran'):
            self.set_fortran_parameters(solution)
        else: raise Exception('Only Fortran kernels are supported in 2D.')


