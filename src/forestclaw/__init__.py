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

# Module imports
__all__.extend(['Patch'])
from clawpack.forestclaw.geometry import Patch
