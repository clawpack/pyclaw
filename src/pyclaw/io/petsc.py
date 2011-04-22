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

import pyclaw.solution

logger = logging.getLogger('io')

from petsc4py import PETSc
#from mpi4py import MPI
import pickle
import numpy as np
    

def write_petsc(solution,frame,path='./',file_prefix='claw',write_aux=False,options={}):
    r"""
        Write out pickle and PETSc data files representing the
        solution.  Common data is written from process 0 in pickle
        files.  Shared data is written from all processes into PETSc
        data files.
    
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) pyclaw
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
    option_defaults = {'format':'binary','clobber':True}

    for (k,v) in option_defaults.iteritems():
        if options.has_key(k):
            pass
        else:
            options[k] = option_defaults[k]

    clobber = options['clobber']
    
    pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    if options['format']=='vtk':
        viewer_filename = os.path.join(path, file_prefix+str(frame).zfill(4)+'.vtk')
    else:
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
    rank =  PETSc.Comm.getRank(PETSc.COMM_WORLD)
    if rank==0:
        pickle_file = open(pickle_filename,'wb')
        # explicitly dumping a dictionary here to help out anybody trying to read the pickle file
        pickle.dump({'t':solution.t,'meqn':solution.meqn,'ngrids':len(solution.grids),
                     'maux':solution.maux,'ndim':solution.ndim,'write_aux':write_aux}, pickle_file)

    # now set up the PETSc viewers
    if options['format'] == 'ascii':
        viewer = PETSc.Viewer().createASCII(viewer_filename, PETSc.Viewer.Mode.WRITE)
        if write_aux:
            aux_viewer = PETSc.Viewer().createASCII(aux_filename, PETSc.Viewer.Mode.WRITE) 
    elif options['format'] == 'binary':
        viewer = PETSc.Viewer().createBinary(viewer_filename, PETSc.Viewer.Mode.WRITE)
        if write_aux:
            aux_viewer = PETSc.Viewer().createBinary(aux_filename, PETSc.Viewer.Mode.WRITE)
    elif options['format'] == 'vtk':
        viewer = PETSc.Viewer().createASCII(viewer_filename, PETSc.Viewer.Mode.WRITE, format=PETSc.Viewer.Format.ASCII_VTK)
        if write_aux:
            aux_viewer = PETSc.Viewer().createASCII(aux_filename, PETSc.Viewer.Mode.WRITE) 
    else:
        raise IOError('format type %s not supported' % options['format'])
    
    for grid in solution.grids:
        if rank==0:
            pickle.dump({'gridno':grid.gridno,'level':grid.level,
                         'names':grid.name,'lower':grid.lower,
                         'n':grid.n,'d':grid.d}, pickle_file)

        grid.gqVec.view(viewer)
        
        if write_aux:
            grid.gauxVec.view(aux_viewer)
    
    viewer.flush()
    viewer.destroy()
    if rank==0:
        pickle_file.close()

def read_petsc(solution,frame,path='./',file_prefix='claw',read_aux=False,options={}):
    r"""
    Read in pickles and PETSc data files representing the solution
    
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Solution object to 
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
    option_defaults = {'format':'binary'}

    for (k,v) in option_defaults.iteritems():
        if options.has_key(k):
            pass
        else:
            options[k] = option_defaults[k]
    
    pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    viewer_filename = os.path.join(path, '%s.ptc' % file_prefix) + str(frame).zfill(4)
    aux_filename = os.path.join(path, '%s_aux.ptc' % file_prefix) + str(frame).zfill(4)

    if frame < 0:
        # Don't construct file names with negative frameno values.
        raise IOError("Frame " + str(frame) + " does not exist ***")

    pickle_file = open(pickle_filename,'rb')

    # this dictionary is mostly holding debugging information, only ngrids is needed
    # most of this information is explicitly saved in the individual grids
    value_dict = pickle.load(pickle_file)
    ngrids   = value_dict['ngrids']                    
    read_aux = value_dict['write_aux']
    ndim     = value_dict['ndim']

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
        grid_dict = pickle.load(pickle_file)

        gridno  = grid_dict['gridno']
        level   = grid_dict['level']
        names   = grid_dict['names']
        lower   = grid_dict['lower']
        n       = grid_dict['n']
        d       = grid_dict['d']
                
        dimensions = []
        for i in xrange(ndim):
            dimensions.append(
                pyclaw.solution.Dimension(names[i],lower[i],lower[i] + n[i]*d[i],n[i]))
        grid = pyclaw.solution.Grid(dimensions)

        grid.t = value_dict['t']
        grid.meqn = value_dict['meqn']

        grid.q_da = PETSc.DA().create(
            dim=grid.ndim,
            dof=grid.meqn, # should be modified to reflect the update
            sizes=grid.n, 
            #periodic_type = PETSc.DA.PeriodicType.X,
            #periodic_type=grid.PERIODIC,
            #stencil_type=grid.STENCIL,
            stencil_width=grid.mbc,
            comm=PETSc.COMM_WORLD)

        grid.gqVec = PETSc.Vec().load(viewer)
        q = grid.gqVec.getArray().copy()
        q_dim=[grid.meqn]
        if len(grid.n)==2: 
            q_dim.append(grid.n[1])
            q_dim.append(grid.n[0])
            q = q.reshape(q_dim,order='F')
            #q = np.rollaxis(q,2,1)
        else: 
            q_dim.append(grid.n[0])
            q = q.reshape(q_dim,order='F')
        grid.q = q
        
        if read_aux:
            nbc = [x+(2*grid.mbc) for x in grid.n]
            grid.aux_da = PETSc.DA().create(dim=grid.ndim,
                dof=maux, # should be modified to reflect the update
                sizes=nbc,  #Amal: what about for 2D, 3D
                #periodic_type = PETSc.DA.PeriodicType.X,
                #periodic_type=grid.PERIODIC,
                #stencil_type=grid.STENCIL,
                stencil_width=grid.mbc,
                comm=PETSc.COMM_WORLD)
            grid.gauxVec = PETSc.Vec().load(aux_viewer)
            grid.aux = grid.gauxVec.getArray().copy()
            grid.aux.shape = (grid.aux.size/grid.meqn,grid.meqn)
            
        # Add AMR attributes:
        grid.gridno = gridno
        grid.level = level 
        solution.grids.append(grid)

    pickle_file.close()
    viewer.destroy()
    if read_aux:
        aux_viewer.destroy()

def read_petsc_t(frame,path='./',file_prefix='claw'):
    r"""Read only the petsc.pkl file and return the data
    
    :Input:
     - *frame* - (int) Frame number to be read in
     - *path* - (string) Path to the current directory of the file
     - *file_prefix* - (string) Prefix of the files to be read in.  
       ``default = 'claw'``
     
    :Output:
     - (list) List of output variables
      - *t* - (int) Time of frame
      - *meqn* - (int) Number of equations in the frame
      - *ngrids* - (int) Number of grids
      - *maux* - (int) Auxillary value in the frame
      - *ndim* - (int) Number of dimensions in q and aux
    
    """

    base_path = os.path.join(path,)
    path = os.path.join(base_path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    try:
        f = open(path,'rb')
        logger.debug("Opening %s file." % path)
        grid_dict = pickle.load(f)

        t = grid_dict['t']
        meqn = grid_dict['meqn']
        ngrids   = grid_dict['ngrids']                    
        maux   = grid_dict['maux']                    
        ndim     = grid_dict['ndim']

        f.close()
    except(IOError):
        raise
    except:
        logger.error("File " + t_fname + " should contain t, meqn, ngrids, maux, ndim")
        print "File " + t_fname + " should contain t, meqn, ngrids, maux, ndim"
        raise
        
    return t,meqn,ngrids,maux,ndim
