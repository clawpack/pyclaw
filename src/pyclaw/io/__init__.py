#!/usr/bin/env python
# encoding: utf-8
#  ======================================================================
#  Package:     pyclaw.io
#  File:        __init__.py
#  Created:     Feb 10, 2008
#  Author:      Kyle Mandli
#  ======================================================================
"""Output package for Pyclaw"""

import logging

import ascii 
__all__ = ['ascii.read','ascii.write']

from .binary import read_binary
__all__ += ['read_binary']

# Check for HDF 5 support
try:
    import h5py
    import hdf5.read,hdf5.write
    __all__ += ['hdf5.read','hdf5.write']
except ImportError:
    logging.debug("No hdf5 support found.")
    
# Check for netcdf support
try:
    import netCDF4
    import netcdf.read, netcdf.write
    __all__ += ['netcdf.read','netcdf.write']
except(ImportError):
    logging.debug("No netcdf4 support found.")
