"""Output package for Pyclaw"""

import logging
from clawpack.pyclaw.io import ascii

__all__ = ['ascii.read','ascii.write']

try:
    import h5py
    from clawpack.pyclaw.io import hdf5
    __all__ += ['hdf5.read','hdf5.write']
except:
    logging.debug("No hdf5 support found.")

try:
    import netCDF4
    from clawpack.pyclaw.io import netcdf
    __all__ += ['netcdf.read','netcdf.write']
except(ImportError):
    logging.debug("No netcdf4 support found.")

try:
    import multifab
    __all__ += ['multifab.write']
except(ImportError):
    logging.debug("No multifab io support found.")

