#  Authors:      Kyle Mandli
#  ======================================================================
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
__all__.extend(['Controller','Data','Dimension','Grid','Solution','State','riemann','plot'])
from controller import Controller
from data import Data
from solution import Solution
from grid import Dimension, Grid
from state import State
import riemann
import plot

__all__.extend(['ClawSolver1D','ClawSolver2D','SharpClawSolver1D','SharpClawSolver2D'])
from clawpack import ClawSolver1D, ClawSolver2D
from sharpclaw import SharpClawSolver1D, SharpClawSolver2D


# Sub-packages
import limiters
from limiters import *
__all__.extend(limiters.__all__)

__all__.append('BC')
from solver import BC
