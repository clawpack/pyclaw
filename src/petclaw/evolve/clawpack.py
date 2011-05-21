#!/usr/bin/env python
# encoding: utf-8
r"""
Module containg the PetClaw solvers

This module contains the pure and wrapped PetClaw solvers.  All 
PetClaw solvers inherit from the :class:`PetClawSolver` superclass which in turn 
inherits from the :class:`~petclaw.evolve.solver.PetSolver` superclass.  As such, 
the only solver classes that should be directly used should be the 
dimensionally dependent ones such as :class:`PetClawSolver1D`.

:Authors:
    Amal Alghamdi
    David Ketcheson
"""
# ============================================================================
#      Copyright (C) 2010 David I. Ketcheson <david.ketcheson@kaust.edu.sa>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

from petclaw.evolve.solver import PetSolver
from pyclaw.evolve.clawpack import ClawSolver, ClawSolver1D, ClawSolver2D

# ============================================================================
#  PetClaw 1d Solver Class
# ============================================================================
class PetClawSolver1D(PetSolver,ClawSolver1D):
    r"""
    PetClaw solver for 1D problems using classic Clawpack algorithms.

    This class implements nothing; it just inherits from PetClawSolver and
    ClawSolver1D.
    """

# ============================================================================
#  PetClaw 2d Solver Class
# ============================================================================
class PetClawSolver2D(PetSolver,ClawSolver2D):
    r"""
    PetClaw solver for 2D problems using classic Clawpack algorithms.

    This class implements nothing; it just inherits from PetClawSolver and
    ClawSolver2D.
    
    Note that only the fortran routines are supported for now in 2D.
    """
