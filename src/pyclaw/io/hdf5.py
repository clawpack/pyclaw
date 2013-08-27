#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing a HDF5 output file

This module reads and writes hdf5 files via either of the following modules:
    h5py - http://code.google.com/p/h5py/
    PyTables - http://www.pytables.org/moin

It will first try h5py and then PyTables and use the correct calls
according to whichever is present on the system.  We recommend that you use
h5py as it is a minimal wrapper to the HDF5 library and will create 

To install either, you must also install the hdf5 library from the website:
    http://www.hdfgroup.org/HDF5/release/obtain5.html
    
:Authors:
    Kyle T. Mandli (2009-02-13) Initial version
"""
# ============================================================================
#      Copyright (C) 2009 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import os
import logging

import clawpack.pyclaw.solution

logger = logging.getLogger('io')

# Import appropriate hdf5 package
use_h5py = False
use_PyTables = False
try:
    import h5py
    use_h5py = True
except:
    pass
if not use_h5py:
    try:
        import tables
        use_PyTables = True
    except:
        error_msg = ("Could not import h5py or PyTables, please install " +
            "either h5py or PyTables.  See the doc_string for more " +
            "information.")
        raise Exception(error_msg)

if not use_h5py and not use_PyTables:
    logging.critical("Could not import h5py or PyTables!")

def write(solution,frame,path,file_prefix='claw',write_aux=False,
                options={},write_p=False):
    r"""
    Write out a Solution to a HDF5 file.
    
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Pyclaw solution 
       object to input into
     - *frame* - (int) Frame number
     - *path* - (string) Root path
     - *file_prefix* - (string) Prefix for the file name.  ``default = 'claw'``
     - *write_aux* - (bool) Boolean controlling whether the associated 
       auxiliary array should be written out.  ``default = False``     
     - *options* - (dict) Optional argument dictionary, see 
       `HDF5 Option Table`_
    
    .. _`HDF5 Option Table`:
    
    +-----------------+------------------------------------------------------+
    | Key             | Value                                                |
    +=================+======================================================+
    | compression     | (None, string ["gzip" | "lzf" | "szip"] or int 0-9)  |
    |                 | Enable dataset compression. DEFLATE, LZF and (where  |
    |                 | available) SZIP are supported. An integer is         |
    |                 | interpreted as a GZIP level for backwards            |
    |                 | compatibility.                                       |
    +-----------------+------------------------------------------------------+
    |compression_opts | (None, or special value) Setting for compression     |
    |                 | filter; legal values for each filter type are:       |
    |                 |                                                      |
    |                 | - *gzip* - (int) 0-9                                 |
    |                 | - *lzf* - None allowed                               |
    |                 | - *szip* - (tuple) 2-tuple ('ec'|'nn', even integer  |
    |                 |     0-32)                                            |
    |                 |                                                      |
    |                 | See the filters module for a detailed description of |
    |                 | each of these filters.                               |
    +-----------------+------------------------------------------------------+
    | chunks          | (None, True or shape tuple) Store the dataset in     |
    |                 | chunked format. Automatically selected if any of the |
    |                 | other keyword options are given. If you don't provide|
    |                 | a shape tuple, the library will guess one for you.   |
    +-----------------+------------------------------------------------------+
    | shuffle         | (True/False) Enable/disable data shuffling, which can|
    |                 | improve compression performance. Automatically       |
    |                 | enabled when compression is used.                    |
    +-----------------+------------------------------------------------------+
    | fletcher32      | (True/False) Enable Fletcher32 error detection; may  |
    |                 | be used with or without compression.                 |
    +-----------------+------------------------------------------------------+
    """
    
    # Option parsing
    option_defaults = {'compression':None,'compression_opts':None,
                       'chunks':None,'shuffle':False,'fletcher32':False}
    for (k,v) in option_defaults.iteritems():
        if options.has_key(k):
            exec("%s = options['%s']" % (k,k))
        else:
            exec('%s = v' % k)
    
    # File name
    filename = os.path.join(path,'%s%s.hdf' % 
                                (file_prefix,str(frame).zfill(4)))
    
    # Write out using h5py
    if use_h5py:
        f = h5py.File(filename,'w')
        
        # For each patch, write out attributes
        for state in solution.states:
            patch = state.patch
            # Create group for this patch
            subgroup = f.create_group('patch%s' % patch.patch_index)
            
            # General patch properties
            for attr in ['t','num_eqn','num_ghost','patch_index','level']:
                if hasattr(patch,attr):
                    if getattr(patch,attr) is not None:
                        subgroup.attrs[attr] = getattr(patch,attr)
                    
            # Add the dimension names as a attribute
            subgroup.attrs['dimensions'] = patch.get_dim_attribute('name')
            # Dimension properties
            for dim in patch.dimensions:
                for attr in ['n','lower','d','upper','bc_lower',
                             'bc_upper','units','num_cells']:
                    if hasattr(dim,attr):
                        if getattr(dim,attr) is not None:
                            attr_name = '%s.%s' % (dim.name,attr)
                            subgroup.attrs[attr_name] = getattr(dim,attr)
            
            # Write out q
            if write_p:
                q = state.p
            else:
                q = state.q
            subgroup.create_dataset('q',data=q,
                                        compression=compression,
                                        compression_opts=compression_opts,
                                        chunks=chunks,shuffle=shuffle,
                                        fletcher32=fletcher32)
            if write_aux and patch.num_aux > 0:
                subgroup.create_dataset('aux',data=patch.aux,
                                        compression=compression,
                                        compression_opts=compression_opts,
                                        chunks=chunks,shuffle=shuffle,
                                        fletcher32=fletcher32)
    
        # Flush and close the file
        f.close()
        
    # Write out using PyTables
    elif use_PyTables:
        # f = tables.openFile(filename, mode = "w", title = options['title'])
        logging.critical("PyTables has not been implemented yet.")
        raise IOError("PyTables has not been implemented yet.")
    else:
        err_msg = "No hdf5 python modules available."
        logging.critical(err_msg)
        raise Exception(err_msg)

def read(solution,frame,path='./',file_prefix='claw',read_aux=True,
                options={}):
    r"""
    Read in a HDF5 file into a Solution
    
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
    for (k,v) in option_defaults.iteritems():
        if options.has_key(k):
            exec("%s = options['%s']" % (k,k))
        else:
            exec('%s = v' % k)
    
    # File name
    filename = os.path.join(path,'%s%s.hdf' % 
                                (file_prefix,str(frame).zfill(4)))

    if use_h5py:
        f = h5py.File(filename,'r')
        
        for subgroup in f.iterobjects():

            # Construct each dimension
            dimensions = []
            dim_names = subgroup.attrs['dimensions']
            for dim_name in dim_names:
                # Create dimension
                dim = pyclaw.solution.Dimension(dim_name,
                                    subgroup.attrs["%s.lower" % dim_name],
                                    subgroup.attrs["%s.upper" % dim_name],
                                    subgroup.attrs["%s.n" % dim_name])                    
                # Optional attributes
                for attr in ['bc_lower','bc_upper','units']:
                    attr_name = "%s.%s" % (dim_name,attr)
                    if subgroup.attrs.get(attr_name, None):
                        setattr(dim,attr,subgroup.attrs["%s.%s" % (dim_name,attr)])
                dimensions.append(dim)
            
            # Create patch
            patch = pyclaw.solution.Patch(dimensions)
                
            # Fetch general patch properties
            for attr in ['t','num_eqn','patch_index','level']:
                setattr(patch,attr,subgroup.attrs[attr])
            
            # Read in q
            index_str = ','.join( [':' for i in xrange(len(subgroup['q'].shape))] )
            exec("patch.q = subgroup['q'][%s]" % index_str)
            
            # Read in aux if applicable
            if read_aux and subgroup.get('aux',None) is not None:
                index_str = ','.join( [':' for i in xrange(len(subgroup['aux'].shape))] )
                exec("patch.aux = subgroup['aux'][%s]" % index_str)
                
            solution.patchs.append(patch)
            
        # Flush and close the file
        f.close()
            
    elif use_PyTables:
        # f = tables.openFile(filename, mode = "r", title = options['title'])
        logging.critical("PyTables has not been implemented yet.")
        raise IOError("PyTables has not been implemented yet.")
    else:
        err_msg = "No hdf5 python modules available."
        logging.critical(err_msg)
        raise Exception(err_msg)
        
