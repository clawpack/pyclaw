#!/usr/bin/env python
# encoding: utf-8

"""I/O support for ForestClaw"""

from __future__ import absolute_import
import logging
logger = logging.getLogger('pyclaw.fileio')

from . import ascii
__all__ = ["ascii.read", "ascii.write"]
