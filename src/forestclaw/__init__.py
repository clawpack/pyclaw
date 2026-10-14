"""Additions for ForestClaw support"""
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
