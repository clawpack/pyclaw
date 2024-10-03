"""Main petclaw package"""
import os
import logging, logging.config

# Default logging configuration file
_DEFAULT_LOG_CONFIG_PATH = os.path.join(os.path.dirname(__file__),'../pyclaw/log.config')
del os

# Setup loggers
logging.config.fileConfig(_DEFAULT_LOG_CONFIG_PATH)

__all__ = []

# Module imports
__all__.extend(['Controller','Dimension','Patch','Domain','Solution','State','CFL','riemann'])
from .controller import Controller
from clawpack.petclaw.geometry import Patch, Domain 
from clawpack.pyclaw.geometry import Dimension
from .solution import Solution
from .state import State
from .cfl import CFL

__all__.extend(['ClawSolver1D','ClawSolver2D','ClawSolver3D','SharpClawSolver1D','SharpClawSolver2D','SharpClawSolver3D'])
from .classic.solver import ClawSolver1D,ClawSolver2D,ClawSolver3D
from .sharpclaw.solver import SharpClawSolver1D,SharpClawSolver2D,SharpClawSolver3D

__all__.append('BC')
from clawpack.pyclaw.solver import BC

# Sub-packages
from . import limiters
from .limiters import *
__all__.extend(limiters.__all__)

from . import plot
__all__.append('plot')
