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
from clawpack import PetClawSolver1D as ClawSolver1D
from clawpack import PetClawSolver2D as ClawSolver2D
from sharpclaw import PetSharpClawSolver1D as SharpClawSolver1D
from sharpclaw import PetSharpClawSolver2D as SharpClawSolver2D
from riemann import *

__all__.append('limiters')
from pyclaw.evolve import limiters
