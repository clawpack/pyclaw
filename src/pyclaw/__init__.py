#  =====================================================================
#  Package:     pyclaw
#  File:        __init__.py
#  Created:     Feb 19, 2008
#  Author:      Kyle Mandli
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
__all__.extend(['Controller','Data','Dimension','Grid','Solution','riemann'])
from controller import Controller
from data import Data
from solution import Solution
from grid import Dimension, Grid
import riemann

# Sub-packages
import evolve
from evolve import *
__all__.extend(evolve.__all__)
