r"""
Module containing the PetClaw solvers

This file currently only exists so that these solvers have a different
__module__ property, used by pyclaw.solver.Solver.__init__ to
determine the containing claw_package to use.
"""

from __future__ import absolute_import
from clawpack import pyclaw

class ClawSolver1D(pyclaw.ClawSolver1D):
    r"""
    Parallel solver for 1D problems using classic Clawpack algorithms.
    """

    __doc__ += pyclaw.util.add_parent_doc(pyclaw.ClawSolver1D)

class ClawSolver2D(pyclaw.ClawSolver2D):
    r"""
    Parallel solver for 2D problems using classic Clawpack algorithms.
    """

    __doc__ += pyclaw.util.add_parent_doc(pyclaw.ClawSolver2D)
    
class ClawSolver3D(pyclaw.ClawSolver3D):
    r"""
    Parallel solver for 3D problems using classic Clawpack algorithms.
    """

    __doc__ += pyclaw.util.add_parent_doc(pyclaw.ClawSolver3D)
