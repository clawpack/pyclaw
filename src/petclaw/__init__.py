#  =====================================================================
#  Package:     petclaw
#  File:        __init__.py
#  Authors:     Amal Alghamdi
#               David Ketcheson
#               Aron Ahmadia
#  ======================================================================
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
__all__.extend(['Controller','Data','Dimension','Grid','Solution','State','CFL','riemann'])
from controller import Controller
from grid import Dimension
from pyclaw.grid import Grid 
from pyclaw.data import Data
from pyclaw.solution import Solution
from state import State
from cfl import CFL

__all__.extend(['ClawSolver1D','ClawSolver2D','ClawSolver3D','SharpClawSolver1D','SharpClawSolver2D'])
from clawpack import ClawSolver1D,ClawSolver2D,ClawSolver3D
from sharpclaw import SharpClawSolver1D,SharpClawSolver2D
from implicitclawpack import ImplicitClawSolver1D

__all__.append('BC')
from pyclaw.solver import BC

# Sub-packages
import limiters
from limiters import *
__all__.extend(limiters.__all__)

import plot
__all__.append('plot')
