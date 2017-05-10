#!/usr/bin/env python
# encoding: utf-8
"""Output package for PetClaw"""

from __future__ import absolute_import
import logging
from clawpack.pyclaw.fileio import ascii
__all__ = ['ascii.read','ascii.write']

# Check for HDF 5 support
try:
    import h5py
    from . import hdf5
    __all__ += ['hdf5.read','hdf5.write']
except:
    logging.debug("No hdf5 support found.")
    
# Check for netcdf support
try:
    import netCDF4
    from clawpack.pyclaw.fileio import netcdf 
    __all__ += ['netcdf.read','netcdf.write']
except(ImportError):
    logging.debug("No netcdf4 support found.")

# Check for petsc4py support
try:
    from . import petsc
    __all__ += ['petsc.read','petsc.write']
except(ImportError):
    logging.debug("No petsc support found.")
 
