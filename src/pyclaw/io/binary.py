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

logger = logging.getLogger('pyclaw.io')


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
     - *read_aux* (bool) Whether or not an auxillary file will try to be read 
       in.  ``default = False``
     - *options* - (dict) Dictionary of optional arguments dependent on 
       the format being read in.  ``default = {}``
    """
    
    # Construct path names
    base_path = os.path.join(path,)
    q_fname = os.path.join(base_path, '%s.q' % file_prefix) + str(frame).zfill(4)
    b_fname = os.path.join(base_path, '%s.b' % file_prefix) + str(frame).zfill(4)

    # Read in values from fort.t file:
    [t,num_eqn,nstates,num_aux,num_dim,num_ghost] = read_t(frame,path,file_prefix)

    patches = []
    
    # Read in values from fort.q file:
    try:
        b_file = open(b_fname,'rb')
    except IOError:
        print "Error: file " + b_fname + " does not exist or is unreadable."
        raise IOError("Could not read binary file %s" % b_fname)

    qdata = np.fromfile(file=b_file, dtype=np.float64)

    i_start_patch = 0  # index into qdata for start of next patch

    try:
        f = open(q_fname,'r')
    except(IOError):
        print "Error: file " + q_fname + " does not exist or is unreadable."
        raise

   
    # Loop through every patch setting the appropriate information
    # for ng in range(len(solution.patchs)):
    for m in xrange(nstates):
    
        # Read in base header for this patch
        patch_index = read_data_line(f,data_type=int)
        level = read_data_line(f,data_type=int)
        n = np.zeros((num_dim))
        lower = np.zeros((num_dim))
        d = np.zeros((num_dim))
        for i in xrange(num_dim):
            n[i] = read_data_line(f,data_type=int)
        for i in xrange(num_dim):
            lower[i] = read_data_line(f)
        for i in xrange(num_dim):
            d[i] = read_data_line(f)
    
        blank = f.readline()
    
        # Construct the patch
        # Since we do not have names here, we will construct the patch with
        # the assumed dimensions x,y,z
        names = ['x','y','z']
        dimensions = []
        for i in xrange(num_dim):
            dimensions.append(
                pyclaw.geometry.Dimension(names[i],lower[i],lower[i] + n[i]*d[i],n[i]))
        patch = pyclaw.geometry.Patch(dimensions)
        state= pyclaw.state.State(patch,num_eqn,num_aux)
        state.t = t

        if num_aux > 0:   
            state.aux[:]=0.
        
        # Fill in q values
        if patch.num_dim == 1:
            ##  NOT YET TESTED ##
            mx = patch.dimensions[0].num_cells
            meqn = state.num_eqn
            mbc = num_ghost
            i_end_patch = i_start_patch + meqn*(mx+2*mbc)
            qpatch = qdata[i_start_patch:i_end_patch]
            qpatch = np.reshape(qpatch, (meqn,mx+2*mbc), \
                        order='F')
            state.q = qpatch[:,mbc:-mbc]
            i_start_patch = i_end_patch  # prepare for next patch

        elif patch.num_dim == 2:
            ## FIXED FOR BINARY ##
            mx = patch.dimensions[0].num_cells
            my = patch.dimensions[1].num_cells
            meqn = state.num_eqn
            mbc = num_ghost
            i_end_patch = i_start_patch + meqn*(mx+2*mbc)*(my+2*mbc)
            qpatch = qdata[i_start_patch:i_end_patch]
            qpatch = np.reshape(qpatch, (meqn,mx+2*mbc,my+2*mbc), \
                        order='F')
            state.q = qpatch[:,mbc:-mbc,mbc:-mbc]
            i_start_patch = i_end_patch  # prepare for next patch

        elif patch.num_dim == 3:
            ##  NOT YET TESTED ##
            mx = patch.dimensions[0].num_cells
            my = patch.dimensions[1].num_cells
            mz = patch.dimensions[2].num_cells
            meqn = state.num_eqn
            mbc = num_ghost
            i_end_patch = i_start_patch + \
                        meqn*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)
            qpatch = qdata[i_start_patch:i_end_patch]
            qpatch = np.reshape(qpatch, \
                        (meqn,mx+2*mbc,my+2*mbc,mz+2*mbc), \
                        order='F')
            state.q = qpatch[:,mbc:-mbc,mbc:-mbc,mbc:-mbc]
            i_start_patch = i_end_patch  # prepare for next patch

        else:
            msg = "Read only supported up to 3d."
            logger.critical(msg)
            raise Exception(msg)
    
        # Add AMR attributes:
        patch.patch_index = patch_index
        patch.level = level

        # Add new patch to solution
        solution.states.append(state)
        patches.append(state.patch)
    solution.domain = pyclaw.geometry.Domain(patches)
        
    #-------------
    # aux file:
    #-------------

    # Read auxillary file if available and requested
    if solution.states[0].num_aux > 0 and read_aux:
        # Check for aux file
        fname1 = os.path.join(base_path,'%s.a' % file_prefix)+str(frame).zfill(4)
        fname2 = os.path.join(base_path,'%s.a' % file_prefix)+str(0).zfill(4)
        if os.path.exists(fname1):
            fname = fname1
        elif os.path.exists(fname2):
            fname = fname2
        else:
            logger.debug("Unable to open auxillary file %s or %s" % (fname1,fname2))
            return
            
        # Found a valid path, try to open and read it
        try:
            b_file = open(fname,'rb')
            auxdata = np.fromfile(file=b_file, dtype=np.float64)
        except IOError:
            print "Error: file " + fname + " does not exist or is unreadable."
            raise IOError("Could not read binary file %s" % fname)

        i_start_patch = 0  # index into auxdata for start of next patch
        for state in solution.states:
            patch = state.patch

            # Fill in aux values
            if patch.num_dim == 1:
                ##  NOT YET TESTED ##
                mx = patch.dimensions[0].num_cells
                maux = state.num_aux
                mbc = num_ghost
                i_end_patch = i_start_patch + maux*(mx+2*mbc)
                auxpatch = auxdata[i_start_patch:i_end_patch]
                auxpatch = np.reshape(auxpatch, (maux,mx+2*mbc), \
                            order='F')
                state.aux = auxpatch[:,mbc:-mbc]
                i_start_patch = i_end_patch  # prepare for next patch

            elif patch.num_dim == 2:
                ## FIXED FOR BINARY ##
                mx = patch.dimensions[0].num_cells
                my = patch.dimensions[1].num_cells
                maux = state.num_aux
                mbc = num_ghost
                i_end_patch = i_start_patch + maux*(mx+2*mbc)*(my+2*mbc)
                auxpatch = auxdata[i_start_patch:i_end_patch]
                auxpatch = np.reshape(auxpatch, (maux,mx+2*mbc,my+2*mbc), \
                            order='F')
                state.aux = auxpatch[:,mbc:-mbc,mbc:-mbc]
                i_start_patch = i_end_patch  # prepare for next patch

            elif patch.num_dim == 3:
                ##  NOT YET TESTED ##
                mx = patch.dimensions[0].num_cells
                my = patch.dimensions[1].num_cells
                mz = patch.dimensions[2].num_cells
                maux = state.num_aux
                mbc = num_ghost
                i_end_patch = i_start_patch + \
                            maux*(mx+2*mbc)*(my+2*mbc)*(mz+2*mbc)
                auxpatch = auxdata[i_start_patch:i_end_patch]
                auxpatch = np.reshape(auxpatch, \
                            (maux,mx+2*mbc,my+2*mbc,mz+2*mbc), \
                            order='F')
                state.aux = auxpatch[:,mbc:-mbc,mbc:-mbc,mbc:-mbc]
                i_start_patch = i_end_patch  # prepare for next patch

            else:
                logger.critical("Read aux only up to 3d is supported.")
                raise Exception("Read aux only up to 3d is supported.")

            
def read_t(frame,path='./',file_prefix='fort'):
    r"""Read only the fort.t file and return the data


    Note that this version reads in the extra value for num_ghost so that we
    can extract only the data that's relevant.
    
    :Input:
     - *frame* - (int) Frame number to be read in
     - *path* - (string) Path to the current directory of the file
     - *file_prefix* - (string) Prefix of the files to be read in.  
       ``default = 'fort'``
     
    :Output:
     - (list) List of output variables
     - *t* - (int) Time of frame
     - *num_eqn* - (int) Number of equations in the frame
     - *nstates* - (int) Number of states
     - *num_aux* - (int) Auxillary value in the frame
     - *num_dim* - (int) Number of dimensions in q and aux
     - *num_ghost* - (int) Number of ghost cells on each side
    
    """

    base_path = os.path.join(path,)
    path = os.path.join(base_path, '%s.t' % file_prefix) + str(frame).zfill(4)
    try:
        logger.debug("Opening %s file." % path)
        f = open(path,'r')
        
        t = read_data_line(f)
        num_eqn = read_data_line(f, data_type=int)
        nstates = read_data_line(f, data_type=int)
        num_aux = read_data_line(f, data_type=int)
        num_dim = read_data_line(f, data_type=int)
        num_ghost = read_data_line(f, data_type=int)
        
        f.close()
    except(IOError):
        raise
    except:
        logger.error("File " + path + " should contain t, num_eqn, nstates, num_aux, num_dim")
        print "File " + path + " should contain t, num_eqn, nstates, num_aux, num_dim"
        raise
        
    return t,num_eqn,nstates,num_aux,num_dim,num_ghost

