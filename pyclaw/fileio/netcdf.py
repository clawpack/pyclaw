#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing a NetCDF output file

Routines for reading and writing a NetCDF output file via either
    - netcdf4-python - http://code.google.com/p/netcdf4-python/
    - pupynere - http://pypi.python.org/pypi/pupynere/
    
These interfaces are very similar so if a different module needs to be used,
it can more than likely be inserted with a minimal of effort.

This module will first try to import the netcdf4-python module which is based
on the compiled libraries and failing that will attempt to import the pure
python interface pupynere which requires no libraries.

To install the netCDF 4 library, please see:
    http://www.unidata.ucar.edu/software/netcdf/
    
:Authors:
    Kyle T. Mandli (2009-02-17) Initial version
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

from __future__ import absolute_import
from __future__ import print_function
import os,sys
import logging

import clawpack.pyclaw.solution
import six
from six.moves import range

logger = logging.getLogger('pyclaw.fileio')

# Import appropriate netcdf package
use_netcdf4 = False
use_pupynere = False
try:
    import netCDF4
    use_netcdf4 = True
except:
    pass
if not use_netcdf4:
    try:
        import pupynere
        use_pupynere = True
    except:
        error_msg = ("Could not import netCDF4 or Pupynere, please install " +
            "one of the available modules for netcdf files.  Refer to this " +
            "modules doc_string for more information.")
        #raise Exception(error_msg)
        print(error_msg)

def write(solution,frame,path,file_prefix='claw',write_aux=False,
                    options={},write_p=False):
    r"""
    Write out a NetCDF data file representation of solution
    
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Pyclaw object to be 
       output
     - *frame* - (int) Frame number
     - *path* - (string) Root path
     - *file_prefix* - (string) Prefix for the file name. ``default = 'claw'``
     - *write_aux* - (bool) Boolean controlling whether the associated 
       auxiliary array should be written out. ``default = False``     
     - *options* - (dict) Optional argument dictionary, see 
       `NetCDF Option Table`_
    
    .. _`NetCDF Option Table`:
    
    +-------------------------+----------------------------------------------+
    | Key                     | Value                                        |
    +=========================+==============================================+
    | description             | Dictionary of key/value pairs that will be   |
    |                         | attached to the root group as attributes,    |
    |                         | i.e. {'time':3}                              |
    +-------------------------+----------------------------------------------+
    | format                  | Can be one of the following netCDF flavors:  |
    |                         | NETCDF3_CLASSIC, NETCDF3_64BIT,              |
    |                         | NETCDF4_CLASSIC, and NETCDF4                 |
    |                         | ``default = NETCDF4``                        |
    +-------------------------+----------------------------------------------+
    | clobber                 | if True (Default), file will be overwritten, |
    |                         | if False an exception will be raised         |
    +-------------------------+----------------------------------------------+
    | zlib                    | if True, data assigned to the Variable       |
    |                         | instance is compressed on disk.              |
    |                         | ``default = False``                          |
    +-------------------------+----------------------------------------------+
    | complevel               | the level of zlib compression to use (1 is   |
    |                         | the fastest, but poorest compression, 9 is   |
    |                         | the slowest but best compression).  Ignored  |
    |                         | if zlib=False.  ``default = 6``              |
    +-------------------------+----------------------------------------------+
    | shuffle                 | if True, the HDF5 shuffle filter is applied  |
    |                         | to improve compression. Ignored if           |
    |                         | zlib=False. ``default = True``               |
    +-------------------------+----------------------------------------------+
    | fletcher32              | if True (default False), the Fletcher32      |
    |                         | checksum algorithm is used for error         |
    |                         | detection.                                   |
    +-------------------------+----------------------------------------------+
    | contiguous              | if True (default False), the variable data   |
    |                         | is stored contiguously on disk.  Setting to  |
    |                         | True for a variable with an unlimited        |
    |                         | dimension will trigger an error.             |
    |                         | ``default = False``                          |
    +-------------------------+----------------------------------------------+
    | chunksizes              | Can be used to specify the HDF5 chunksizes   |
    |                         | for each dimension of the variable. A        |
    |                         | detailed discussion of HDF chunking and I/O  |
    |                         | performance is available here. Basically,    |
    |                         | you want the chunk size for each dimension   |
    |                         | to match as closely as possible the size of  |
    |                         | the data block that users will read from the |
    |                         | file. chunksizes cannot be set if            |
    |                         | contiguous=True.                             |
    +-------------------------+----------------------------------------------+
    | least_significant_digit | If specified, variable data will be          |
    |                         | truncated (quantized). In conjunction with   |
    |                         | zlib=True this produces 'lossy', but         |
    |                         | significantly more efficient compression.    |
    |                         | For example, if least_significant_digit=1,   |
    |                         | data will be quantized using around          |
    |                         | (scale*data)/scale, where scale = 2**bits,   |
    |                         | and bits is determined so that a precision   |
    |                         | of 0.1 is retained (in this case bits=4).    |
    |                         | ``default = None``, or no quantization.      |
    +-------------------------+----------------------------------------------+
    | endian                  | Can be used to control whether the data is   |
    |                         | stored in little or big endian format on     | 
    |                         | disk. Possible values are little, big or     |
    |                         | native (default). The library will           |
    |                         | automatically handle endian conversions when |
    |                         | the data is read, but if the data is always  |
    |                         | going to be read on a computer with the      |
    |                         | opposite format as the one used to create    |
    |                         | the file, there may be some performance      |
    |                         | advantage to be gained by setting the        |
    |                         | endian-ness.                                 |
    +-------------------------+----------------------------------------------+
    | fill_value              | If specified, the default netCDF _FillValue  |
    |                         | (the value that the variable gets filled     |
    |                         | with before any data is written to it) is    |
    |                         | replaced with this value. If fill_value is   |
    |                         | set to False, then the variable is not       |
    |                         | pre-filled.                                  |
    +-------------------------+----------------------------------------------+
    
    .. note:: 
        The zlib, complevel, shuffle, fletcher32, contiguous, chunksizes and
        endian keywords are silently ignored for netCDF 3 files that do not 
        use HDF5.
        
    """
    
    # Option parsing
    option_defaults = {'format':'NETCDF4','zlib':False,'complevel':6,
                       'shuffle':True,'fletcher32':False,'contiguous':False,
                       'chunksizes':None,'endian':'native',
                       'least_significant_digit':None,'fill_value':None,
                       'clobber':True,'description':{}}
    for (k,v) in six.iteritems(option_defaults):
        if k in options:
            exec("%s = options['%s']" % (k,k))
        else:
            exec('%s = v' % k)
            
    # Filename
    filename = os.path.join(path,"%s%s.nc" % (file_prefix,str(frame).zfill(4)))
        
    if use_netcdf4:
        # Open new file
        f = netCDF4.Dataset(filename,'w',clobber=clobber,format=format)
        
        # Loop through description dictionary and add the attributes to the
        # root group
        for (k,v) in six.iteritems(description):
            exec('f.%s = %s' % (k,v))
        
        # For each patch, write out attributes
        for state in solution.states:
            patch = solution.patch
            # Create group for this patch
            subgroup = f.createGroup('patch%s' % patch.patch_index)
        
            # General patch properties
            for attr in ['t','num_eqn']:
                setattr(subgroup,attr,getattr(state,attr))
            for attr in ['patch_index','level']:
                setattr(subgroup,attr,getattr(patch,attr))
            
            # Write out dimension names
            setattr(subgroup,'dim_names',patch.name)
            
            # Create dimensions for q (and aux)
            for dim in patch.dimensions:
                subgroup.createDimension(dim.name,dim.num_cells)
                # Write other dimension attributes
                for attr in ['num_cells','lower','delta','upper','bc_lower',
                             'bc_upper','units']:
                    if hasattr(dim,attr):
                        if getattr(dim,attr) is not None:
                            attr_name = '%s.%s' % (dim.name,attr)
                            setattr(subgroup,attr_name,getattr(dim,attr))
            subgroup.createDimension('num_eqn',state.num_eqn)
            
            # Write q array
            from copy import copy
            dim_names = copy(patch.name)
            dim_names.append('num_eqn')
            index_str = ','.join( [':' for name in dim_names] )
            q = subgroup.createVariable('q','f8',dim_names,zlib,
                                            complevel,shuffle,fletcher32,
                                            contiguous,chunksizes,endian,
                                            least_significant_digit,fill_value)
            exec("q[%s] = state.q" % index_str)
            
            # Write out aux
            if state.num_aux > 0 and write_aux:
                dim_names[-1] = 'num_aux'
                subgroup.createDimension('num_aux',state.num_aux)
                aux = subgroup.createVariable('aux','f8',dim_names,
                                            zlib,complevel,shuffle,fletcher32,
                                            contiguous,chunksizes,endian,
                                            least_significant_digit,fill_value)
                exec("aux[%s] = state.aux" % index_str)
        
        f.close()
    elif use_pupynere:
        logging.critical("Pupynere support has not been implemented yet.")
        raise IOError("Pupynere support has not been implemented yet.")
    else:
        err_msg = "No netcdf python modules available."
        logging.critical(err_msg)
        raise Exception(err_msg)

    
def read(solution,frame,path='./',file_prefix='claw',read_aux=True,
                options={}):
    r"""
    Read in a NetCDF data files into solution
    
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Pyclaw object to be 
       output
     - *frame* - (int) Frame number
     - *path* - (string) Root path
     - *file_prefix* - (string) Prefix for the file name.  ``default = 'claw'``
     - *write_aux* - (bool) Boolean controlling whether the associated 
       auxiliary array should be written out.  ``default = False``     
     - *options* - (dict) Optional argument dictionary, unused for reading.
    """
    
    # Option parsing
    option_defaults = {}
    for (k,v) in six.iteritems(option_defaults):
        if k in options:
            exec("%s = options['%s']" % (k,k))
        else:
            exec('%s = v' % k)
            
    # Filename
    filename = os.path.join(path,"%s%s.nc" % (file_prefix,str(frame).zfill(4)))
        
    if use_netcdf4:
        # Open file
        f = netCDF4.Dataset(filename,'r')
        
        # We only expect subgroups of patches, otherwise we need to put some
        # sort of conditional here
        for subgroup in six.itervalues(f.groups):
            # Construct each dimension
            dimensions = []
            
            # Read in dimension attribute to keep dimension order
            dim_names = getattr(subgroup,'dim_names')
            for dim_name in dim_names:
                dim = pyclaw.solution.Dimension(
                                      getattr(subgroup,'%s.lower' % dim_name),
                                      getattr(subgroup,'%s.upper' % dim_name),
                                      getattr(subgroup,'%s.n' % dim_name),
                                      name = dim_name)
                 # Optional attributes
                for attr in ['bc_lower','bc_upper','units']:
                    attr_name = "%s.%s" % (dim_name,attr)
                    if hasattr(subgroup,attr_name):
                        setattr(dim,attr,getattr(subgroup, "%s.%s" % (dim_name,attr)))
                dimensions.append(dim)
            
            # Create patch
            patch = pyclaw.solution.Patch(dimensions)
            
            # General patch properties
            for attr in ['t','num_eqn','patch_index','level']:
                setattr(patch,attr,getattr(subgroup,attr))
                
            # Read in q
            index_str = ','.join( [':' for i in range(patch.num_dim+1)] )
            exec("patch.q = subgroup.variables['q'][%s]" % index_str)
            
            # Read in aux if applicable
            if read_aux and 'num_aux' in subgroup.dimensions:
                exec("patch.aux = subgroup.variables['aux'][%s]" % index_str)
        
            solution.patches.append(patch)
            
        f.close()
    elif use_pupynere:
        logging.critical("Pupynere support has not been implemented yet.")
        raise IOError("Pupynere support has not been implemented yet.")
    else:
        err_msg = "No netcdf python modules available."
        logging.critical(err_msg)
        raise Exception(err_msg)
