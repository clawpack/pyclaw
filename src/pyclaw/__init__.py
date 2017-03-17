"""Main pyclaw package"""

from __future__ import absolute_import
import os
import logging, logging.config

# Default logging configuration file
_DEFAULT_LOG_CONFIG_PATH = os.path.join(os.path.dirname(__file__),'log.config')
del os

# Setup loggers
logging.config.fileConfig(_DEFAULT_LOG_CONFIG_PATH)

__all__ = []

# Module imports
__all__.extend(['Controller','Dimension','Patch','Domain','Solution','State','CFL','plot'])
from .controller import Controller
from .solution import Solution
from .geometry import Dimension, Patch, Domain
from .state import State
from .cfl import CFL

from clawpack.pyclaw import plot

__all__.extend(['ClawSolver1D','ClawSolver2D','ClawSolver3D','SharpClawSolver1D','SharpClawSolver2D','SharpClawSolver3D'])
from .classic.solver import ClawSolver1D, ClawSolver2D, ClawSolver3D
from .sharpclaw.solver import SharpClawSolver1D, SharpClawSolver2D,SharpClawSolver2D,SharpClawSolver3D


# Sub-packages
from . import limiters
from .limiters import *
__all__.extend(limiters.__all__)

__all__.append('BC')
from .solver import BC
__all__.extend('IOTest')
from .tests.test_io import IOTest
