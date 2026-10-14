#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing a HDF5 output file

This module reads and writes hdf5 files via the following module:
    h5py - http://code.google.com/p/h5py/

To install h5py, you must also install the hdf5 library from the website:
    http://www.hdfgroup.org/HDF5/release/obtain5.html
"""
from mpi4py import MPI
import os
import logging

from clawpack.petclaw import geometry
from clawpack import petclaw
import six

logger = logging.getLogger('pyclaw.fileio')

try:
    import h5py
except:
    logging.critical("Could not import h5py!")
    error_msg = ("Could not import h5py, please install " +
                 "either h5py.  See the doc_string for more " +
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
    for (k,v) in six.iteritems(option_defaults):
        options[k] = options.get(k,v)
    
    filename = os.path.join(path,'%s%s.hdf' % 
                                (file_prefix,str(frame).zfill(4)))

    if options['compression'] is not None:
        err_msg = "Compression (filters) are not available for parallel h5py yet."
        logging.critical(err_msg)
        raise Exception(err_msg)

    with h5py.File(filename,'w',driver='mpio',comm=MPI.COMM_WORLD) as f:
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
            subgroup.attrs['dimensions'] = [ name.encode('utf-8')
                        for name in patch.get_dim_attribute('name') ]
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
            r = patch._da.getRanges()
            globalSize = []
            globalSize.append(q.shape[0])
            globalSize.extend(patch.num_cells_global)
            dset = subgroup.create_dataset('q',globalSize,dtype='float',**options)
            to_hdf5_dataset(q, dset, len(patch.dimensions), r)

            if write_aux and state.num_aux > 0:
                r = patch._da.getRanges()
                globalSize = []
                globalSize.append(state.num_aux)
                globalSize.extend(patch.num_cells_global)
                dset = subgroup.create_dataset('aux',globalSize,dtype='float',**options)
                to_hdf5_dataset(state.aux, dset, len(patch.dimensions), r)


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

    with h5py.File(filename,'r',driver='mpio',comm=MPI.COMM_WORLD) as f:
        for patch in six.itervalues(f):
            # Construct each dimension
            dimensions = []
            dim_names = [ name.decode('ascii')
                          for name in patch.attrs['dimensions'] ]
            for dim_name in dim_names:
                dim = geometry.Dimension(
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

            pyclaw_patch = petclaw.Patch(dimensions)

            # Fetch general patch properties
            for attr in ['t','num_eqn','patch_index','level']:
                setattr(pyclaw_patch,attr,patch.attrs[attr])

            state = petclaw.state.State(pyclaw_patch, \
                                        patch.attrs['num_eqn'],patch.attrs['num_aux'])
            state.t = patch.attrs['t']

            globalSize = []
            globalSize.append(state.q.shape[0])
            globalSize.extend(pyclaw_patch.num_cells_global)
            r = pyclaw_patch._da.getRanges()

            dset = patch['q'][:].reshape(globalSize)
            state.q = from_hdf5_dataset(dset, len(pyclaw_patch.dimensions), r, state.q.shape)

            # Read in aux if applicable
            if read_aux and patch.get('aux',None) is not None:
                dset = patch['aux'][:]
                state.aux = from_hdf5_dataset(dset, len(pyclaw_patch.dimensions), r, state.aux.shape)

            solution.states.append(state)
            patches.append(pyclaw_patch)

        solution.domain = geometry.Domain(patches)

def to_hdf5_dataset(arr, dset, ndim, ranges):
    if ndim == 1:
        dset[:,ranges[0][0]:ranges[0][1]] = arr
    elif ndim == 2:
        dset[:,ranges[0][0]:ranges[0][1],ranges[1][0]:ranges[1][1]] = arr
    elif ndim == 3:
        dset[:,ranges[0][0]:ranges[0][1],ranges[1][0]:ranges[1][1],ranges[2][0]:ranges[2][1]] = arr

def from_hdf5_dataset(dset, ndim, ranges, shape):
    if ndim == 1:
        return dset[:,ranges[0][0]:ranges[0][1]].reshape(shape,order='F')
    elif ndim == 2:
        return dset[:,ranges[0][0]:ranges[0][1],ranges[1][0]:ranges[1][1]].reshape(shape,order='F')
    elif ndim == 3:
        return dset[:,ranges[0][0]:ranges[0][1],ranges[1][0]:ranges[1][1],ranges[2][0]:ranges[2][1]].reshape(shape,order='F')
