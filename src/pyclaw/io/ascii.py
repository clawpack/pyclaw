#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing an ascii output file

:Authors:
    Kyle T. Mandli (2008-02-07) Initial version
"""
# ============================================================================
#      Copyright (C) 2008 Kyle T. Mandli <mandli@amath.washington.edu>
#
#  Distributed under the terms of the Berkeley Software Distribution (BSD) 
#  license
#                     http://www.opensource.org/licenses/
# ============================================================================

import os,sys
import logging

from pyclaw.util import read_data_line
import pyclaw.solution

logger = logging.getLogger('io')

def write_ascii(solution,frame,path,file_prefix='fort',write_aux=False,
                    options={},write_p=False):
    r"""
    Write out ascii data file
    
    Write out an ascii file formatted identical to the fortran clawpack files
    including writing out fort.t, fort.q, and fort.aux if necessary.  Note
    that there are some parameters that assumed to be the same for every grid
    in this format which is not necessarily true for the actual data objects.
    Make sure that if you use this output format that all of you grids share
    the appropriate values of ndim, meqn, maux, and t.  Only supports up to
    3 dimensions.
    
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Pyclaw object to be 
       output.
     - *frame* - (int) Frame number
     - *path* - (string) Root path
     - *file_prefix* - (string) Prefix for the file name.  ``default = 'fort'``
     - *write_aux* - (bool) Boolean controlling whether the associated 
       auxiliary array should be written out.  ``default = False``
     - *options* - (dict) Dictionary of optional arguments dependent on 
       the format being written.  ``default = {}``
    """
    try:
        # Create file name
        file_name = '%s.t%s' % (file_prefix,str(frame).zfill(4))
        f = open(os.path.join(path,file_name),'w')
        
        # Header for fort.txxxx file
        f.write("%18.8e     time\n" % solution.t)
        f.write("%5i                  meqn\n" % solution.meqn)
        f.write("%5i                  nstates\n" % len(solution.states))
        f.write("%5i                  maux\n" % solution.maux)
        f.write("%5i                  ndim\n" % solution.ndim)
        f.close()
        
        # Open fort.qxxxx for writing
        file_name = 'fort.q%s' % str(frame).zfill(4)
        q_file = open(os.path.join(path,file_name),'w')
        
        # If maux != 0 then we open up a file to write it out as well
        if solution.maux > 0 and write_aux:
            file_name = 'fort.a%s' % str(frame).zfill(4)
            aux_file = open(os.path.join(path,file_name),'w')
        
        # for i in range(0,len(solution.grids)):
        for state in solution.states:
            grid = state.grid
            # Header for fort.qxxxx file
            q_file.write("%5i                  grid_number\n" % grid.gridno)
            q_file.write("%5i                  AMR_level\n" % grid.level)
            for dim in grid.dimensions:
                q_file.write("%5i                  m%s\n" % (dim.n,dim.name))
            for dim in grid.dimensions:
                q_file.write("%18.8e     %slow\n" % (dim.lower,dim.name))
            for dim in grid.dimensions:
                q_file.write("%18.8e     d%s\n" % (dim.d,dim.name))
            
            q_file.write("\n")
            
            # Write data from q
            if write_p:
                q = state.p
            else:
                q = state.q

            dims = grid.dimensions
            if grid.ndim == 1:
                for k in xrange(dims[0].n):
                    for m in xrange(solution.meqn):
                        q_file.write("%18.8e" % q[m,k])
                    q_file.write('\n')
            elif grid.ndim == 2:
                for j in xrange(dims[1].n):
                    for k in xrange(dims[0].n):
                        for m in xrange(solution.meqn):
                            q_file.write("%18.8e" % q[m,k,j])
                        q_file.write('\n')
                    q_file.write('\n')
            elif grid.ndim == 3:
                for l in xrange(dims[2].n):
                    for j in xrange(dims[1].n):
                        for k in xrange(dims[0].n):
                            for m in range(solution.meqn):
                                q_file.write("%18.8e" % q[m,k,j,l])
                            q_file.write('\n')
                    q_file.write('\n')
                q_file.write('\n')
            else:
                raise Exception("Dimension Exception in writing fort file.")
            
            if state.maux > 0 and write_aux:
                aux = state.aux
                
                aux_file.write("%5i                  grid_number\n" % grid.gridno)
                aux_file.write("%5i                  AMR_level\n" % grid.level)
                
                for dim in grid.dimensions:
                    aux_file.write("%5i                  m%s\n" % (dim.n,dim.name))
                for dim in grid.dimensions:
                    aux_file.write("%18.8e     %slow\n" % (dim.lower,dim.name))
                for dim in grid.dimensions:
                    aux_file.write("%18.8e     d%s\n" % (dim.d,dim.name))

                aux_file.write("\n")
                if grid.ndim == 1:
                    for k in xrange(grid.dimensions[0]):
                        for m in xrange(state.maux):
                            aux_file.write("%18.8e" % aux[m,k])
                        aux_file.write('\n')
                elif grid.ndim == 2:
                    for j in xrange(grid.dimensions[1].n):
                        for k in xrange(grid.dimension[0].n):
                            for m in xrange(state.maux):
                                aux_file.write("%18.8e" % aux[m,k,j])
                            aux_file.write('\n')    
                        aux_file.write('\n')
                elif grid.ndim == 3:
                    for l in xrange(grid.dimensions[2].n):
                        for j in xrange(grid.dimensions[1].n):
                            for k in xrange(grid.dimensions[0].n):
                                for m in xrange(grid.maux):
                                    aux_file.write("%18.8e" % aux[m,k,j,l])
                                aux_file.write('\n')
                            aux_file.write('\n')    
                        aux_file.write('\n')
    
        q_file.close()
        if state.maux > 0 and write_aux:
            aux_file.close()
    except IOError, (errno, strerror):
        logger.error("Error writing file: %s" % os.path.join(path,file_name))
        logger.error("I/O error(%s): %s" % (errno, strerror))
        raise 
    except:
        logger.error("Unexpected error:", sys.exc_info()[0])
        raise

def read_ascii(solution,frame,path='./',file_prefix='fort',read_aux=False,
                options={}):
    r"""
    Read in a set of ascii formatted files
    
    This routine reads the ascii formatted files corresponding to the classic
    clawpack format 'fort.txxxx', 'fort.qxxxx', and 'fort.axxxx' or 'fort.aux'
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
    
    import numpy as np

    if frame < 0:
        # Don't construct file names with negative frameno values.
        raise IOError("Frame " + str(frame) + " does not exist ***")

    # Construct path names
    base_path = os.path.join(path,)
    # t_fname = os.path.join(base_path, '%s.t' % file_prefix) + str(frame).zfill(4)
    q_fname = os.path.join(base_path, '%s.q' % file_prefix) + str(frame).zfill(4)

    # Read in values from fort.t file:
    [t,meqn,nstates,maux,ndim] = read_ascii_t(frame,path,file_prefix)
    
    # Read in values from fort.q file:
    try:
        f = open(q_fname,'r')
    
        # Loop through every grid setting the appropriate information
        # for ng in range(len(solution.grids)):
        for m in xrange(nstates):
        
            # Read in base header for this grid
            gridno = read_data_line(f,type='int')
            level = read_data_line(f,type='int')
            n = np.zeros((ndim))
            lower = np.zeros((ndim))
            d = np.zeros((ndim))
            for i in xrange(ndim):
                n[i] = read_data_line(f,type='int')
            for i in xrange(ndim):
                lower[i] = read_data_line(f)
            for i in xrange(ndim):
                d[i] = read_data_line(f)
        
            blank = f.readline()
        
            # Construct the grid
            # Since we do not have names here, we will construct the grid with
            # the assumed dimensions x,y,z
            names = ['x','y','z']
            dimensions = []
            for i in xrange(ndim):
                dimensions.append(
                    pyclaw.grid.Dimension(names[i],lower[i],lower[i] + n[i]*d[i],n[i]))
            grid = pyclaw.grid.Grid(dimensions)
            state= pyclaw.state.State(grid,meqn,maux)
            state.t = t


            # RJL 1/8/10:  Changed empty_aux to zeros_aux below so aux won't
            # be filled with random values if aux arrays not read in.  Would
            # like to delete this and initialize grid.aux only if it will be
            # read in below, but for some reason that doesn't work.

            if maux > 0:   
                state.aux[:]=0.
            
            # Fill in q values
            if grid.ndim == 1:
                for i in xrange(grid.dimensions[0].n):
                    l = []
                    while len(l)<state.meqn:
                        line = f.readline()
                        l = l + line.split()
                    for m in xrange(state.meqn):
                        state.q[m,i] = float(l[m])
            elif grid.ndim == 2:
                for j in xrange(grid.dimensions[1].n):
                    for i in xrange(grid.dimensions[0].n):
                        l = []
                        while len(l)<state.meqn:
                            line = f.readline()
                            l = l + line.split()
                        for m in xrange(state.meqn):
                            state.q[m,i,j] = float(l[m])
                    blank = f.readline()
            elif grid.ndim == 3:
                for k in xrange(grid.dimensions[2].n):
                    for j in xrange(grid.dimensions[1].n):
                        for i in xrange(grid.dimensions[0].n):
                            l=[]
                            while len(l) < state.meqn:
                                line = f.readline()
                                l = l + line.split()
                            for m in xrange(state.meqn):
                                state.q[m,i,j,k] = float(l[m])
                        blank = f.readline()
                    blank = f.readline()
            else:
                msg = "Read only supported up to 3d."
                logger.critical(msg)
                raise Exception(msg)
        
            # Add AMR attributes:
            grid.gridno = gridno
            grid.level = level

            # Add new grid to solution
            solution.states.append(state)
            
    except(IOError):
        raise
    except:
        logger.error("File %s was not able to be read." % q_fname)
        raise

    # Read auxillary file if available and requested
    if solution.states[0].maux > 0 and read_aux:
        # Check for aux file
        fname1 = os.path.join(base_path,'%s.a' % file_prefix)+str(frame).zfill(4)
        fname2 = os.path.join(base_path,'%s.a' % file_prefix)
        if os.path.exists(fname1):
            fname = fname1
        elif os.path.exists(fname2):
            fname = fname2
        else:
            logger.info("Unable to open auxillary file %s or %s" % (fname1,fname2))
            return
            
        # Found a valid path, try to open and read it
        try:
            f = open(fname,'r')
            
            # Read in aux file
            for n in xrange(len(solution.states)):
                # Fetch correct grid
                gridno = read_data_line(f,type='int')
                grid = solution.states[gridno-1].grid
        
                # These should match this grid already, raise exception otherwise
                if not (grid.level == read_data_line(f,type='int')):
                    raise IOError("Grid level in aux file header did not match grid no %s." % grid.gridno)
                for dim in grid.dimensions:
                    if not (dim.n == read_data_line(f,type='int')):
                        raise IOError("Dimension %s's n in aux file header did not match grid no %s." % (dim.name,grid.gridno))
                for dim in grid.dimensions:
                    if not (dim.lower == read_data_line(f,type='float')):
                        raise IOError("Dimension %s's lower in aux file header did not match grid no %s." % (dim.name,grid.gridno))
                for dim in grid.dimensions:
                    if not (dim.d == read_data_line(f,type='float')):
                        raise IOError("Dimension %s's d in aux file header did not match grid no %s." % (dim.name,grid.gridno))

                f.readline()
        
                # Read in auxillary array
                if grid.ndim == 1:
                    for i in xrange(grid.dimensions[0].n):
                        l = []
                        while len(l)<state.maux:
                            line = f.readline()
                            l = l + line.split()
                        for m in xrange(state.maux):
                            state.aux[i,m] = float(l[m])
                elif grid.ndim == 2:
                    for j in xrange(grid.dimensions[1].n):
                        for i in xrange(grid.dimensions[0].n):
                            l = []
                            while len(l)<state.maux:
                                line = f.readline()
                                l = l + line.split()
                            for m in xrange(state.maux):
                                state.aux[i,j,m] = float(l[m])
                        blank = f.readline()
                elif grid.ndim == 3:
                    for k in xrange(grid.dimensions[2].n):
                        for j in xrange(grid.dimensions[1].n):
                            for i in xrange(grid.dimensions[0].n):
                                l = []
                                while len(l)<state.maux:
                                    line = f.readline()
                                    l = l + line.split()
                                for m in xrange(state.maux):
                                    state.aux[i,j,k,m] = float(l[m])
                            blank = f.readline()
                        blank = f.readline()
                else:
                    logger.critical("Read aux only up to 3d is supported.")
                    raise Exception("Read aux only up to 3d is supported.")
        except(IOError):
            raise
        except:
            logger.error("File %s was not able to be read." % q_fname)
            raise
            
            
def read_ascii_t(frame,path='./',file_prefix='fort'):
    r"""Read only the fort.t file and return the data
    
    :Input:
     - *frame* - (int) Frame number to be read in
     - *path* - (string) Path to the current directory of the file
     - *file_prefix* - (string) Prefix of the files to be read in.  
       ``default = 'fort'``
     
    :Output:
     - (list) List of output variables
     - *t* - (int) Time of frame
     - *meqn* - (int) Number of equations in the frame
     - *nstates* - (int) Number of states
     - *maux* - (int) Auxillary value in the frame
     - *ndim* - (int) Number of dimensions in q and aux
    
    """

    base_path = os.path.join(path,)
    path = os.path.join(base_path, '%s.t' % file_prefix) + str(frame).zfill(4)
    try:
        logger.debug("Opening %s file." % path)
        f = open(path,'r')
        
        t = read_data_line(f)
        meqn = read_data_line(f,type='int')
        nstates = read_data_line(f,type='int')
        maux = read_data_line(f,type='int')
        ndim = read_data_line(f,type='int')
        
        f.close()
    except(IOError):
        raise
    except:
        logger.error("File " + path + " should contain t, meqn, nstates, maux, ndim")
        print "File " + path + " should contain t, meqn, nstates, maux, ndim"
        raise
        
    return t,meqn,nstates,maux,ndim
