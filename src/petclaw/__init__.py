"""Main petclaw package"""

import os
import logging, logging.config

# Default logging configuration file
_DEFAULT_LOG_CONFIG_PATH = os.path.join(os.path.dirname(__file__),'log.config')
del os

# Setup loggers
logging.config.fileConfig(_DEFAULT_LOG_CONFIG_PATH)

__all__ = []

# Module imports
__all__.extend(['Controller','Dimension','Patch','Domain','Solution','State','CFL','riemann'])
from petclaw.controller import Controller
from petclaw.geometry import Patch, Domain 
from pyclaw.geometry import Dimension
from petclaw.solution import Solution
from petclaw.state import State
from petclaw.cfl import CFL

__all__.extend(['ClawSolver1D','ClawSolver2D','ClawSolver3D','SharpClawSolver1D','SharpClawSolver2D'])
from petclaw.clawpack.clawpack import ClawSolver1D,ClawSolver2D,ClawSolver3D
from petclaw.sharpclaw.sharpclaw import SharpClawSolver1D,SharpClawSolver2D

__all__.append('BC')
from pyclaw.solver import BC

# Sub-packages
import limiters
from limiters import *
__all__.extend(limiters.__all__)

import plot
__all__.append('plot')
