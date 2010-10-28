#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing a petsc-style output file.

These routines preserve petclaw/pyclaw syntax for i/o while taking advantage of PETSc's parallel i/o capabilities to allow for parallel reads and writes of frame data.
    
:Authors:
    Aron J. Ahmadia (2010-10-26) Initial version
"""

# ============================================================================
#      Copyright (C) 2010 Aron J Ahmadia <aron@ahmadia.net>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import os,sys
import logging

import petclaw.solution

logger = logging.getLogger('io')

from petsc4py import PETSc
from mpi4py import MPI

try:
    import cPickle as pickle
except:
    import pickle
    

def write_petscio(solution,frame,path='./',file_prefix='claw',write_aux=False,options={}):
    r"""
        Write out pickle and PETSc data files representing the
        solution.  Common data is written from process 0 in pickle
        files.  Shared data is written from all processes into PETSc
        data files.
    
    :Input:
     - *solution* - (:class:`~petclaw.solution.Solution`) petclaw
       object to be output
     - *frame* - (int) Frame number
     - *path* - (string) Root path
     - *file_prefix* - (string) Prefix for the file name. ``default =
        'claw'``
     - *write_aux* - (bool) Boolean controlling whether the associated 
       auxiliary array should be written out. ``default = False``     
     - *options* - (dict) Optional argument dictionary, see 
       `PETScIO Option Table`_
    
    .. _`PETScIO Option Table`:
    
    format   : one of 'ascii' or 'binary'
    clobber  : if True (Default), files will be overwritten
    """

    # Option parsing
    option_defaults = {'format':'ascii','clobber':True}

    for (k,v) in option_defaults.iteritems():
        if options.has_key(k):
            pass
        else:
            options[k] = option_defaults[k]

    clobber = options['clobber']

    if frame < 0:
        # Don't construct file names with negative frameno values.
        raise IOError("Frame " + str(frame) + " does not exist ***")

    pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    viewer_filename = os.path.join(path, '%s.ptc' % file_prefix) + str(frame).zfill(4)

    if solution.maux > 0 and write_aux:
        write_aux = True
        aux_filename = os.path.join(path, '%s_aux.ptc' % file_prefix) + str(frame).zfill(4)
    else:
        write_aux = False
        
    if not clobber:
        if os.path.exists(pickle_filename):
            raise IOError('Cowardly refusing to clobber %s!' % pickle_filename)
        if os.path.exists(viewer_filename):
            raise IOError('Cowardly refusing to clobber %s!' % viewer_filename)
        if write_aux and os.path.exists(aux_filename):
            raise IOError('Cowardly refusing to clobber %s!' % aux_filename)

    rank = MPI.COMM_WORLD.rank
    if rank==0:
        pickle_file = open(pickle_filename,'w')
        # explicitly dumping a dictionary here to help out anybody trying to read the pickle file
        pickle.dump({'t':solution.t,'meqn':solution.meqn,'ngrids':len(solution.grids),
                     'maux':solution.maux,'ndim':solution.ndim, 'write_aux':write_aux}, pickle_file)
        
    # now set up the PETSc viewers
    if options['format'] == 'ascii':
        viewer = PETSc.Viewer().createASCII(viewer_filename, PETSc.Viewer.Mode.WRITE)
        if write_aux:
            aux_viewer = PETSc.Viewer().createASCII(aux_filename, PETSc.Viewer.Mode.WRITE) 
    elif options['format'] == 'binary':
        viewer = PETSc.Viewer().createBinary(viewer_filename, PETSc.Viewer.Mode.WRITE)
        if write_aux:
            aux_viewer = PETSc.Viewer().createBinary(aux_filename, PETSc.Viewer.Mode.WRITE)
    else:
        raise IOError('format type %s not supported' % options['format'])
    
    for grid in solution.grids:
        if rank==0:
            pickle.dump(grid, pickle_file)

        grid.gqVec.view(viewer)
        
        if write_aux:
            grid.gauxVec.view(aux_viewer)
    
    viewer.flush()
    viewer.destroy()
    if rank==0:
        pickle_file.close()

def read_petscio(solution,frame,path='./',file_prefix='claw',read_aux=False,options={}):
    r"""
    Read in pickles and PETSc data files representing the solution
    
    :Input:
     - *solution* - (:class:`~petclaw.solution.Solution`) Solution object to 
       read the data into.
     - *frame* - (int) Frame number to be read in
     - *path* - (string) Path to the current directory of the file
     - *file_prefix* - (string) Prefix of the files to be read in.  
       ``default = 'fort'``
     - *read_aux* (bool) Whether or not an auxillary file will try to be read 
       in.  ``default = False``
     - *options* - (dict) Optional argument dictionary, see 
       `PETScIO Option Table`_
    
    .. _`PETScIO Option Table`:
    
    format   : one of 'ascii' or 'binary'
     
    """

    # Option parsing
    option_defaults = {'format':'ascii'}

    for (k,v) in option_defaults.iteritems():
        if options.has_key(k):
            pass
        else:
            options[k] = option_defaults[k]
    
    pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    viewer_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    aux_filename = os.path.join(path, '%s_aux.ptc' % file_prefix) + str(frame).zfill(4)

    in_file = open(filename,'r')
    in_dict = pickle.load(in_file)

    if frame < 0:
        # Don't construct file names with negative frameno values.
        raise IOError("Frame " + str(frame) + " does not exist ***")

    pickle_file = open(pickle_filename,'r')

    # this dictionary is mostly holding debugging information, only ngrids is needed
    # most of this information is explicitly saved in the individual grids
    value_dict = pickle.load(pickle_file)
    ngrids = value_dict['ngrids']                    
    read_aux = value_dict['write_aux']

    # now set up the PETSc viewer
    if options['format'] == 'ascii':
        viewer = PETSc.Viewer().createASCII(viewer_filename, PETSc.Viewer.Mode.READ)
        if read_aux:
            aux_viewer = PETSc.Viewer().createASCII(aux_viewer_filename, PETSc.Viewer.Mode.READ)
    elif options['format'] == 'binary':
        viewer = PETSc.Viewer().createBinary(viewer_filename, PETSc.Viewer.Mode.READ)
        if read_aux:
            aux_viewer = PETSc.Viewer().createBinary(aux_viewer_filename, PETSc.Viewer.Mode.READ)
    else:
        raise IOError('format type %s not supported' % options['format'])

    for m in xrange(ngrids):
        grid = pickle.load(pickle_file)
        grid.gqVec = PETSc.Vec.load(viewer)
        if read_aux:
            grid.gauxVec = PETSc.Vec.load(aux_viewer)

        # still need to rebuild the DAs on the fly
        solution.grids.append(grid)
