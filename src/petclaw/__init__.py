#!/usr/bin/env python
# encoding: utf-8
#  =====================================================================
#  Package:     petclaw
#  File:        __init__.py
#  Created:     Feb 19, 2008
#  Author:      Kyle Mandli
#  ======================================================================
"""Main petclaw package"""

import os
import logging, logging.config

# Default logging configuration file
_DEFAULT_LOG_CONFIG_PATH = os.path.join(os.path.dirname(__file__),'log.config')

# Setup loggers
logging.config.fileConfig(_DEFAULT_LOG_CONFIG_PATH)

__all__ = []

# Module imports
__all__.extend(['Controller','Data','Dimension','Grid','Solution'])
from petclaw.controller import Controller
from petclaw.data import Data
from petclaw.solution import Dimension, Grid, Solution

# Sub-packages
import evolve
from petclaw.evolve import *
__all__.extend(evolve.__all__)
