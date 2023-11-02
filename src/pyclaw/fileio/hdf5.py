#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing a HDF5 output file

This module reads and writes hdf5 files via either of the following modules:
    h5py - http://code.google.com/p/h5py/
    PyTables - http://www.pytables.org/moin

It will first try h5py and then PyTables and use the correct calls
according to whichever is present on the system.  We recommend that you use
h5py as it is a minimal wrapper to the HDF5 library.

To install either, you must also install the hdf5 library from the website:
    http://www.hdfgroup.org/HDF5/release/obtain5.html
"""

import os
import logging

from clawpack import pyclaw
import numpy as np

logger = logging.getLogger('pyclaw.fileio')

# Import appropriate hdf5 package
use_h5py = False
use_PyTables = False
try:
    import h5py
    use_h5py = True
except:
    try:
        import tables
        use_PyTables = True
    except:
        logging.critical("Could not import h5py or PyTables!")
        error_msg = ("Could not import h5py or PyTables, please install " +
            "either h5py or PyTables.  See the doc_string for more " +
            "information.")
        raise Exception(error_msg)


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
    option_defaults = {'compression':None,'compression_opts':None,
                       'chunks':None,'shuffle':False,'fletcher32':False}
    for (k,v) in option_defaults.items():
        options[k] = options.get(k,v)
    
    filename = os.path.join(path,'%s%s.hdf' % 
                                (file_prefix,str(frame).zfill(4)))
    
    if use_h5py:
        with h5py.File(filename,'w') as f:
        
            # For each patch, write out attributes
            for state in solution.states:
                patch = state.patch
                # Create group for this patch
                subgroup = f.create_group('patch%s' % patch.patch_index)

                # General patch properties
                subgroup.attrs['t'] = state.t
                subgroup.attrs['num_eqn'] = state.num_eqn
                subgroup.attrs['num_aux'] = state.num_aux
                for attr in ['num_ghost','patch_index','level']:
                    if hasattr(patch,attr):
                        if getattr(patch,attr) is not None:
                            subgroup.attrs[attr] = getattr(patch,attr)

                # Add the dimension names as a attribute
                subgroup.attrs['dimensions'] = [name.encode('utf-8') 
                            for name in patch.get_dim_attribute('name')]
                # Dimension properties
                for dim in patch.dimensions:
                    for attr in ['num_cells','lower','delta','upper',
                                 'units']:
                        if hasattr(dim,attr):
                            if getattr(dim,attr) is not None:
                                attr_name = '%s.%s' % (dim.name,attr)
                                subgroup.attrs[attr_name] = getattr(dim,attr)

                if write_p:
                    q = state.p
                else:
                    q = state.q
                subgroup.create_dataset('q',data=q,**options)
                if write_aux and state.num_aux > 0:
                    subgroup.create_dataset('aux',data=state.aux,**options)
        
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
     - *options* - (dict) Optional argument dictionary, not used for reading.
    """
    filename = os.path.join(path,'%s%s.hdf' % 
                                (file_prefix,str(frame).zfill(4)))
    patches = []

    if use_h5py:
        with h5py.File(filename,'r') as f:
        
            for patch in f.values():
                # Construct each dimension
                dimensions = []
                dim_names = np.array(patch.attrs["dimensions"]).astype(str)
                for dim_name in dim_names:
                    dim = pyclaw.solution.Dimension(
                                        patch.attrs["%s.lower" % dim_name],
                                        patch.attrs["%s.upper" % dim_name],
                                        patch.attrs["%s.num_cells" % dim_name],
                                        name = dim_name)
                    # Optional attributes
                    for attr in ['units']:
                        attr_name = "%s.%s" % (dim_name,attr)
                        if patch.attrs.get(attr_name, None):
                            setattr(dim,attr,patch.attrs["%s.%s" % (dim_name,attr)])
                    dimensions.append(dim)

                pyclaw_patch = pyclaw.solution.Patch(dimensions)

                # Fetch general patch properties
                for attr in ['t','num_eqn','patch_index','level']:
                    setattr(pyclaw_patch,attr,patch.attrs[attr])

                state = pyclaw.state.State(pyclaw_patch, \
                         patch.attrs['num_eqn'],patch.attrs['num_aux'])
                state.t = patch.attrs['t']
                state.q = patch['q'][:].ravel(order='F').reshape(state.q.shape,order='F')

                # Read in aux if applicable
                if read_aux and patch.get('aux',None) is not None:
                    state.aux = patch['aux'][:].ravel(order='F').reshape(state.aux.shape,order='F')

                solution.states.append(state)
                patches.append(pyclaw_patch)
                
            solution.domain = pyclaw.geometry.Domain(patches)
                
    elif use_PyTables:
        # f = tables.openFile(filename, mode = "r", title = options['title'])
        logging.critical("PyTables has not been implemented yet.")
        raise IOError("PyTables has not been implemented yet.")
    else:
        err_msg = "No hdf5 python modules available."
        logging.critical(err_msg)
        raise Exception(err_msg)
