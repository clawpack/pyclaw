#!/usr/bin/env python
# encoding: utf-8
#  ======================================================================
#  Package:     pyclaw.io
#  File:        __init__.py
#  Created:     Feb 10, 2008
#  Author:      Kyle Mandli
#  ======================================================================
"""Output package for Pyclaw"""

from __future__ import absolute_import
import logging
logger = logging.getLogger('pyclaw.io')

import h5py

from clawpack.pyclaw.io import ascii 
__all__ = ['ascii.read','ascii.write']

from clawpack.pyclaw.io import binary
__all__ += ['binary.read']

from clawpack.pyclaw.io import hdf5
__all__ += ['hdf5.read','hdf5.write']
    
# Check for netcdf support
try:
    import netCDF4
    from clawpack.pyclaw.io import netcdf
    __all__ += ['netcdf.read','netcdf.write']
except(ImportError):
    logger.debug("No netcdf4 support found.")
