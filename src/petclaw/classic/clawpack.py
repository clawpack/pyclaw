r"""
Module containing the PetClaw solvers

This file currently only exists so that these solvers have a different
__module__ property, used by pyclaw.solver.Solver.__init__ to
determine the containing claw_package to use.
"""

from __future__ import absolute_import
from clawpack.pyclaw.classic import clawpack

class ClawSolver1D(clawpack.ClawSolver1D):
    r"""
    PetClaw solver for 1D problems using classic Clawpack algorithms.

    This class implements nothing; it just inherits from ClawSolver1D.
    """

class ClawSolver2D(clawpack.ClawSolver2D):
    r"""
    PetClaw solver for 2D problems using classic Clawpack algorithms.

    This class implements nothing; it just inherits from ClawSolver2D.
    
    Note that only the fortran routines are supported for now in 2D.
    """
    
class ClawSolver3D(clawpack.ClawSolver3D):
    r"""
    PetClaw solver for 3D problems using classic Clawpack algorithms.

    This class implements nothing; it just inherits from ClawSolver3D.
    
    Note that only fortran routines are supported in 3D.
    """
