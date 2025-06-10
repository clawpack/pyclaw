#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading a raw binary output file from AMRClaw.
Note that there is no corresponding output option in PyClaw,
which is why there is no "write" function here (the code that
writes these files is in AMRClaw, in Fortran).
"""

import os
import logging

from ..util import read_data_line
import numpy as np
import clawpack.pyclaw as pyclaw
from clawpack.pyclaw.fileio.ascii import read_t

logger = logging.getLogger('pyclaw.fileio')


def read(solution,frame,path='./',file_prefix='fort',read_aux=False,
                options={}):
    r"""
    Read in a set of raw binary files
    
    This routine reads the binary formatted files 
    fort.txxxx contains info about frame
    fort.qxxxx still contains headers for each grid patch
    fort.bxxxx is binary dump of data from all patches.
    fort.axxxx is binary dump of aux arrays from all patches.

    Note that the fort prefix can be changed.
    
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Solution object to 
       read the data into.
     - *frame* - (int) Frame number to be read in
     - *path* - (string) Path to the current directory of the file
     - *file_prefix* - (string) Prefix of the files to be read in.  
       ``default = 'fort'``
     - *read_aux* (bool) Whether or not to read auxiliary data from fort.axxxx.
       ``default = False``
     - *options* - (dict) Dictionary of optional arguments dependent on 
       the format being read in.  ``default = {}``
    """
    
    
    # Construct path names
    base_path = os.path.join(path,)
    q_fname = os.path.join(base_path, '%s.q' % file_prefix) + str(frame).zfill(4)
    b_fname = os.path.join(base_path, '%s.b' % file_prefix) + str(frame).zfill(4)

    # Read in values from fort.t file:
    [t,num_eqn,nstates,num_aux,num_dim,num_ghost,file_format] = \
        read_t(frame,path,file_prefix)

    patches = []
    
    # Read in values from fort.b file:
    file_format = options.get('format','binary64')
    
    with open(b_fname,'rb') as b_file:
        if file_format in ['binary', 'binary64']:
            qdata = np.fromfile(file=b_file, dtype=np.float64)
        elif file_format == 'binary32':
            qdata = np.fromfile(file=b_file, dtype=np.float32)
        else:
            msg = "Unrecognized file_format: %s" % file_format
            logger.critical(msg)
            raise Exception(msg)

    i_start_patch = 0  # index into qdata for start of next patch
    n     = np.zeros((num_dim),dtype=int)
    lower = np.zeros((num_dim))
    d     = np.zeros((num_dim))
    # Since we do not have dimension names here, we will construct
    # patches with dimensions named x,y,z
    names = ['x','y','z']

    with open(q_fname,'r') as f:
        # Loop through patches, setting the appropriate information
        for m in range(nstates):
        
            # Read in header for this patch
            patch_index = read_data_line(f,data_type=int)
            level       = read_data_line(f,data_type=int)
            for i in range(num_dim):
                n[i] = read_data_line(f,data_type=int)
            for i in range(num_dim):
                lower[i] = read_data_line(f)
            for i in range(num_dim):
                d[i] = read_data_line(f)

            blank = f.readline()

            # Construct the patch
            dimensions = []
            for i in range(num_dim):
                dimensions.append(
                    pyclaw.geometry.Dimension(lower[i],lower[i] + n[i]*d[i],n[i],name=names[i]))
            patch = pyclaw.geometry.Patch(dimensions)
            state = pyclaw.state.State(patch,num_eqn,num_aux)
            state.t = t

            if num_aux > 0:
                # Write NaNs for now to indicate this is uninitialized
                state.aux[:] = np.nan

            # Fill in q values
            meqn = state.num_eqn
            mbc = num_ghost
            # For some reason, this file format includes the ghost cell values!
            q_shape = [m+2*mbc for m in patch.grid.num_cells]
            q_shape.insert(0,meqn)
            q_size = np.prod(q_shape)
            i_end_patch = i_start_patch + q_size
            qpatch = qdata[i_start_patch:i_end_patch]
            qpatch = np.reshape(qpatch, q_shape, order='F')

            if patch.num_dim == 1:
                ##  NOT YET TESTED ##
                if mbc==0:
                    state.q = qpatch
                else:
                    state.q = qpatch[:,mbc:-mbc]
            elif patch.num_dim == 2:
                ## FIXED FOR BINARY ##
                if mbc==0:
                    state.q = qpatch
                else:
                    state.q = qpatch[:,mbc:-mbc,mbc:-mbc]
            elif patch.num_dim == 3:
                ##  NOT YET TESTED ##
                if mbc==0:
                    state.q = qpatch
                else:
                    state.q = qpatch[:,mbc:-mbc,mbc:-mbc,mbc:-mbc]
            else:
                msg = "Read only supported up to 3d."
                logger.critical(msg)
                raise Exception(msg)

            i_start_patch = i_end_patch  # prepare for next patch

            # Add AMR attributes:
            patch.patch_index = patch_index
            patch.level = level

            # Add new patch to solution
            solution.states.append(state)
            patches.append(state.patch)
        solution.domain = pyclaw.geometry.Domain(patches)


    # Read auxiliary file if available and requested
    if solution.states[0].num_aux > 0 and read_aux:
        # Check for aux file
        fname1 = os.path.join(base_path,'%s.a' % file_prefix)+str(frame).zfill(4)
        fname2 = os.path.join(base_path,'%s.a' % file_prefix)+str(0).zfill(4)
        if os.path.exists(fname1):
            fname = fname1
        elif os.path.exists(fname2):
            fname = fname2
        else:
            logger.debug("Unable to open auxiliary file %s or %s" % (fname1,fname2))
            return
            
        # Found a valid path, try to open and read it
        with open(fname,'rb') as b_file:
            if file_format in ['binary', 'binary64']:
                auxdata = np.fromfile(file=b_file, dtype=np.float64)
            elif file_format == 'binary32':
                auxdata = np.fromfile(file=b_file, dtype=np.float32)
            else:
                msg = "Unrecognized file_format: %s" % file_format
                logger.critical(msg)
                raise Exception(msg)

        i_start_patch = 0  # index into auxdata for start of next patch
        for state in solution.states:
            patch = state.patch
            maux = state.num_aux
            mbc = num_ghost

            aux_shape = [m+2*mbc for m in patch.grid.num_cells]
            aux_shape.insert(0,maux)
            aux_size = np.prod(aux_shape)
            i_end_patch = i_start_patch + aux_size
            auxpatch = auxdata[i_start_patch:i_end_patch]
            auxpatch = np.reshape(auxpatch, aux_shape, order='F')

            # Fill in aux values
            if patch.num_dim == 1:
                ##  NOT YET TESTED ##
                state.aux = auxpatch[:,mbc:-mbc]
            elif patch.num_dim == 2:
                ## FIXED FOR BINARY ##
                state.aux = auxpatch[:,mbc:-mbc,mbc:-mbc]
            elif patch.num_dim == 3:
                ##  NOT YET TESTED ##
                state.aux = auxpatch[:,mbc:-mbc,mbc:-mbc,mbc:-mbc]
            else:
                msg = "Read only supported up to 3d."
                logger.critical(msg)
                raise Exception(msg)

            i_start_patch = i_end_patch  # prepare for next patch

