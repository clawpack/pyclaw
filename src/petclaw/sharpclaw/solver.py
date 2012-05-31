#!/usr/bin/env python
# encoding: utf-8
r"""
Module containing SharpClaw solvers for PetClaw
"""

from __future__ import absolute_import
from clawpack.pyclaw.sharpclaw import solver

class SharpClawSolver1D(solver.SharpClawSolver1D):
    """
    1D parallel SharpClaw solver.  See the corresponding PyClaw class for documentation.
    """

class SharpClawSolver2D(solver.SharpClawSolver2D):
    """
    2D parallel SharpClaw solver.  See the corresponding PyClaw class for documentation.
    """
