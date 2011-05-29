#!/usr/bin/env python
# encoding: utf-8
"""
__init__.py

Created by Kyle Mandli on 2008-08-21.
Copyright (c) 2008 University of Washington. All rights reserved.
"""

# This __init__ script only imports common utilities, most of the import 
# should be done depending on the solver needed

__all__ = ['ClawSolver1D','ClawSolver2D','SharpClawSolver1D','SharpClawSolver2D']
from clawpack import ClawSolver1D, ClawSolver2D
from sharpclaw import SharpClawSolver1D, SharpClawSolver2D

__all__.append('limiters')
import limiters

__all__.append('BC')
from solver import BC
