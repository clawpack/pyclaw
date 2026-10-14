#!/usr/bin/env python
# encoding: utf-8

"""I/O support for ForestClaw"""
import logging
logger = logging.getLogger('pyclaw.fileio')

from . import ascii as forestascii
from clawpack.pyclaw.fileio import ascii

ascii.read_patch_header = forestascii.read_patch_header
ascii.write_patch_header = forestascii.write_patch_header

__all__ = ["ascii.read", "ascii.write"]
