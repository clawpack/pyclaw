"""Additions for ForestClaw support"""

from __future__ import absolute_import
import os
import logging
import logging.config

# Default logging configuration file
_DEFAULT_LOG_CONFIG_PATH = os.path.join(os.path.dirname(__file__),
                                        '../pyclaw/log.config')
del os

# Setup loggers
logging.config.fileConfig(_DEFAULT_LOG_CONFIG_PATH)

__all__ = []

# Module imports - Note the only difference here is the geometry module
__all__.extend(['Controller', 'Dimension', 'Patch', 'Domain', 'Solution',
                'State', 'CFL', 'plot'])
from clawpack.forestclaw.controller import Controller
from clawpack.forestclaw.solution import Solution
from clawpack.forestclaw.geometry import Dimension, Patch, Domain
from clawpack.forestclaw.state import State
from clawpack.forestclaw.cfl import CFL

from clawpack.pyclaw import plot

__all__.extend(['ClawSolver1D', 'ClawSolver2D', 'ClawSolver3D',
                'SharpClawSolver1D', 'SharpClawSolver2D', 'SharpClawSolver3D'])
from .classic.solver import ClawSolver1D, ClawSolver2D, ClawSolver3D
from .sharpclaw.solver import SharpClawSolver1D, SharpClawSolver2D
from .sharpclaw.solver import SharpClawSolver2D, SharpClawSolver3D


# Sub-packages
# Not added yet
# from . import limiters
# from .limiters import *
# __all__.extend(limiters.__all__)

# __all__.append('BC')
# from .solver import BC