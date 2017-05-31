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
                'State', 'CFL'])
from clawpack.pyclaw.controller import Controller
from clawpack.pyclaw.solution import Solution
from clawpack.pyclaw.state import State
from clawpack.pyclaw.cfl import CFL
from clawpack.pyclaw.geometry import Dimension
from clawpack.pyclaw.geometry import Domain
from .geometry import Patch

# The below are not yet needed
# __all__.extend(['ClawSolver1D', 'ClawSolver2D', 'ClawSolver3D',
#                 'SharpClawSolver1D', 'SharpClawSolver2D', 'SharpClawSolver3D'])
# from clawpack.pyclaw.classic.solver import ClawSolver1D, ClawSolver2D, ClawSolver3D
# from clawpack.pyclaw.sharpclaw.solver import SharpClawSolver1D, SharpClawSolver2D
# from clawpack.pyclaw.sharpclaw.solver import SharpClawSolver2D, SharpClawSolver3D


# Sub-packages
# from clawpack.pyclaw import limiters
# from clawpack.pyclaw.limiters import *
# __all__.extend(limiters.__all__)

# __all__.append('BC')
# from clawpack.pyclaw.solver import BC
