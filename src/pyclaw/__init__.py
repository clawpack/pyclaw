"""Main pyclaw package"""

import os
import logging, logging.config

# Default logging configuration file
_DEFAULT_LOG_CONFIG_PATH = os.path.join(os.path.dirname(__file__),'log.config')
del os

# Setup loggers
logging.config.fileConfig(_DEFAULT_LOG_CONFIG_PATH)

__all__ = []

# Module imports
__all__.extend(['Controller','Dimension','Patch','Domain','Solution','State','CFL','riemann','plot','classic1','classic2','classic3','sharpclaw1','sharpclaw2'])
from pyclaw.controller import Controller
from pyclaw.solution import Solution
from pyclaw.geometry import Dimension, Patch, Domain
from pyclaw.state import State
from pyclaw.cfl import CFL
from pyclaw.clawpack import classic1, classic2,classic3
from pyclaw.sharpclaw import sharpclaw1, sharpclaw2

import riemann
import plot

__all__.extend(['ClawSolver1D','ClawSolver2D','ClawSolver3D','SharpClawSolver1D','SharpClawSolver2D'])
from pyclaw.clawpack.clawpack import ClawSolver1D, ClawSolver2D, ClawSolver3D
from pyclaw.sharpclaw.sharpclaw import SharpClawSolver1D, SharpClawSolver2D


# Sub-packages
import limiters
from pyclaw.limiters import *
__all__.extend(limiters.__all__)

__all__.append('BC')
from pyclaw.solver import BC
