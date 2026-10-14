#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing SharpClaw solvers for PetClaw
"""

from clawpack import pyclaw

class SharpClawSolver1D(pyclaw.SharpClawSolver1D):
    """1D parallel SharpClaw solver.
    """

    __doc__ += pyclaw.util.add_parent_doc(pyclaw.SharpClawSolver2D)

    
class SharpClawSolver2D(pyclaw.SharpClawSolver2D):
    """2D parallel SharpClaw solver. 
    """

    __doc__ += pyclaw.util.add_parent_doc(pyclaw.SharpClawSolver2D)


class SharpClawSolver3D(pyclaw.SharpClawSolver3D):
    """3D parallel SharpClaw solver. 
    """

    __doc__ += pyclaw.util.add_parent_doc(pyclaw.SharpClawSolver3D)
