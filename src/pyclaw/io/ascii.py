#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing an ascii output file
"""

import os,sys
import logging

from ..util import read_data_line

logger = logging.getLogger('io')

def write_ascii(solution,frame,path,file_prefix='fort',write_aux=False,
                    options={},write_p=False):
    r"""
    Write out ascii data file
    
    Write out an ascii file formatted identical to the fortran clawpack files
    including writing out fort.t, fort.q, and fort.aux if necessary.  Note
    that there are some parameters that assumed to be the same for every patch
    in this format which is not necessarily true for the actual data objects.
    Make sure that if you use this output format that all of you patchs share
    the appropriate values of num_dim, num_eqn, num_aux, and t.  Only supports up to
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
        f.write("%5i                  num_eqn\n" % solution.num_eqn)
        f.write("%5i                  nstates\n" % len(solution.states))
        f.write("%5i                  num_aux\n" % solution.num_aux)
        f.write("%5i                  num_dim\n" % solution.domain.num_dim)
        f.close()
        
        # Open fort.qxxxx for writing
        file_name = 'fort.q%s' % str(frame).zfill(4)
        q_file = open(os.path.join(path,file_name),'w')
        
        # If num_aux != 0 then we open up a file to write it out as well
        if solution.num_aux > 0 and write_aux:
            file_name = 'fort.a%s' % str(frame).zfill(4)
            aux_file = open(os.path.join(path,file_name),'w')
        
        # for i in range(0,len(solution.patchs)):
        for state in solution.states:
            patch = state.patch
            # Header for fort.qxxxx file
            q_file.write("%5i                  patch_number\n" % patch.patch_index)
            q_file.write("%5i                  AMR_level\n" % patch.level)
            for dim in patch.dimensions:
                q_file.write("%5i                  m%s\n" % (dim.num_cells,dim.name))
            for dim in patch.dimensions:
                q_file.write("%18.8e     %slow\n" % (dim.lower,dim.name))
            for dim in patch.dimensions:
                q_file.write("%18.8e     d%s\n" % (dim.delta,dim.name))
            
            q_file.write("\n")
            
            # Write data from q
            if write_p:
                q = state.p
            else:
                q = state.q

            dims = patch.dimensions
            if patch.num_dim == 1:
                for k in xrange(dims[0].num_cells):
                    for m in xrange(state.num_eqn):
                        q_file.write("%18.8e" % q[m,k])
                    q_file.write('\n')
            elif patch.num_dim == 2:
                for j in xrange(dims[1].num_cells):
                    for k in xrange(dims[0].num_cells):
                        for m in xrange(state.num_eqn):
                            q_file.write("%18.8e" % q[m,k,j])
                        q_file.write('\n')
                    q_file.write('\n')
            elif patch.num_dim == 3:
                for l in xrange(dims[2].num_cells):
                    for j in xrange(dims[1].num_cells):
                        for k in xrange(dims[0].num_cells):
                            for m in range(state.num_eqn):
                                q_file.write("%18.8e" % q[m,k,j,l])
                            q_file.write('\n')
                    q_file.write('\n')
                q_file.write('\n')
            else:
                raise Exception("Dimension Exception in writing fort file.")
            
            if state.num_aux > 0 and write_aux:
                aux = state.aux
                
                aux_file.write("%5i                  patch_number\n" % patch.patch_index)
                aux_file.write("%5i                  AMR_level\n" % patch.level)
                
                for dim in patch.dimensions:
                    aux_file.write("%5i                  m%s\n" % (dim.num_cells,dim.name))
                for dim in patch.dimensions:
                    aux_file.write("%18.8e     %slow\n" % (dim.lower,dim.name))
                for dim in patch.dimensions:
                    aux_file.write("%18.8e     d%s\n" % (dim.delta,dim.name))

                aux_file.write("\n")
                dims = patch.dimensions
                if patch.num_dim == 1:
                    for k in xrange(dims[0].num_cells):
                        for m in xrange(state.num_aux):
                            aux_file.write("%18.8e" % aux[m,k])
                        aux_file.write('\n')
                elif patch.num_dim == 2:
                    for j in xrange(dims[1].num_cells):
                        for k in xrange(dims[0].num_cells):
                            for m in xrange(state.num_aux):
                                aux_file.write("%18.8e" % aux[m,k,j])
                            aux_file.write('\n')    
                        aux_file.write('\n')
                elif patch.num_dim == 3:
                    for l in xrange(dims[2].num_cells):
                        for j in xrange(dims[1].num_cells):
                            for k in xrange(dims[0].num_cells):
                                for m in xrange(state.num_aux):
                                    aux_file.write("%18.8e" % aux[m,k,j,l])
                                aux_file.write('\n')
                            aux_file.write('\n')    
                        aux_file.write('\n')
    
        q_file.close()
        if state.num_aux > 0 and write_aux:
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
    [t,num_eqn,nstates,num_aux,num_dim] = read_ascii_t(frame,path,file_prefix)

    patches = []
    
    # Read in values from fort.q file:
    try:
        f = open(q_fname,'r')
    
        # Loop through every patch setting the appropriate information
        # for ng in range(len(solution.patchs)):
        for m in xrange(nstates):
        
            # Read in base header for this patch
            patch_index = read_data_line(f,type='int')
            level = read_data_line(f,type='int')
            n = np.zeros((num_dim))
            lower = np.zeros((num_dim))
            d = np.zeros((num_dim))
            for i in xrange(num_dim):
                n[i] = read_data_line(f,type='int')
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
            import clawpack.pyclaw as pyclaw
            for i in xrange(num_dim):
                dimensions.append(
                    pyclaw.geometry.Dimension(names[i],lower[i],lower[i] + n[i]*d[i],n[i]))
            patch = pyclaw.geometry.Patch(dimensions)
            state= pyclaw.state.State(patch,num_eqn,num_aux)
            state.t = t


            # RJL 1/8/10:  Changed empty_aux to zeros_aux below so aux won't
            # be filled with random values if aux arrays not read in.  Would
            # like to delete this and initialize patch.aux only if it will be
            # read in below, but for some reason that doesn't work.

            if num_aux > 0:   
                state.aux[:]=0.
            
            # Fill in q values
            if patch.num_dim == 1:
                for i in xrange(patch.dimensions[0].num_cells):
                    l = []
                    while len(l)<state.num_eqn:
                        line = f.readline()
                        l = l + line.split()
                    for m in xrange(state.num_eqn):
                        state.q[m,i] = float(l[m])
            elif patch.num_dim == 2:
                for j in xrange(patch.dimensions[1].num_cells):
                    for i in xrange(patch.dimensions[0].num_cells):
                        l = []
                        while len(l)<state.num_eqn:
                            line = f.readline()
                            l = l + line.split()
                        for m in xrange(state.num_eqn):
                            state.q[m,i,j] = float(l[m])
                    blank = f.readline()
            elif patch.num_dim == 3:
                for k in xrange(patch.dimensions[2].num_cells):
                    for j in xrange(patch.dimensions[1].num_cells):
                        for i in xrange(patch.dimensions[0].num_cells):
                            l=[]
                            while len(l) < state.num_eqn:
                                line = f.readline()
                                l = l + line.split()
                            for m in xrange(state.num_eqn):
                                state.q[m,i,j,k] = float(l[m])
                        blank = f.readline()
                    blank = f.readline()
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
            
    except(IOError):
        raise
    except:
        logger.error("File %s was not able to be read." % q_fname)
        raise

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
            logger.info("Unable to open auxillary file %s or %s" % (fname1,fname2))
            return
            
        # Found a valid path, try to open and read it
        try:
            f = open(fname,'r')
            
            # Read in aux file
            for n in xrange(len(solution.states)):
                # Fetch correct patch
                patch_index = read_data_line(f,type='int')
                patch = solution.states[patch_index-1].patch
        
                # These should match this patch already, raise exception otherwise
                if not (patch.level == read_data_line(f,type='int')):
                    raise IOError("Patch level in aux file header did not match patch no %s." % patch.patch_index)
                for dim in patch.dimensions:
                    if not (dim.num_cells == read_data_line(f,type='int')):
                        raise IOError("Dimension %s's n in aux file header did not match patch no %s." % (dim.name,patch.patch_index))
                for dim in patch.dimensions:
                    if not (dim.lower == read_data_line(f,type='float')):
                        raise IOError("Dimension %s's lower in aux file header did not match patch no %s." % (dim.name,patch.patch_index))
                for dim in patch.dimensions:
                    if not (dim.delta == read_data_line(f,type='float')):
                        raise IOError("Dimension %s's d in aux file header did not match patch no %s." % (dim.name,patch.patch_index))

                f.readline()
        
                # Read in auxillary array
                if patch.num_dim == 1:
                    for i in xrange(patch.dimensions[0].num_cells):
                        l = []
                        while len(l)<state.num_aux:
                            line = f.readline()
                            l = l + line.split()
                        for m in xrange(state.num_aux):
                            state.aux[m,i] = float(l[m])
                elif patch.num_dim == 2:
                    for j in xrange(patch.dimensions[1].num_cells):
                        for i in xrange(patch.dimensions[0].num_cells):
                            l = []
                            while len(l)<state.num_aux:
                                line = f.readline()
                                l = l + line.split()
                            for m in xrange(state.num_aux):
                                state.aux[m,i,j] = float(l[m])
                        blank = f.readline()
                elif patch.num_dim == 3:
                    for k in xrange(patch.dimensions[2].num_cells):
                        for j in xrange(patch.dimensions[1].num_cells):
                            for i in xrange(patch.dimensions[0].num_cells):
                                l = []
                                while len(l)<state.num_aux:
                                    line = f.readline()
                                    l = l + line.split()
                                for m in xrange(state.num_aux):
                                    state.aux[m,i,j,k] = float(l[m])
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
     - *num_eqn* - (int) Number of equations in the frame
     - *nstates* - (int) Number of states
     - *num_aux* - (int) Auxillary value in the frame
     - *num_dim* - (int) Number of dimensions in q and aux
    
    """

    base_path = os.path.join(path,)
    path = os.path.join(base_path, '%s.t' % file_prefix) + str(frame).zfill(4)
    try:
        logger.debug("Opening %s file." % path)
        f = open(path,'r')
        
        t = read_data_line(f)
        num_eqn = read_data_line(f,type='int')
        nstates = read_data_line(f,type='int')
        num_aux = read_data_line(f,type='int')
        num_dim = read_data_line(f,type='int')
        
        f.close()
    except(IOError):
        raise
    except:
        logger.error("File " + path + " should contain t, num_eqn, nstates, num_aux, num_dim")
        print "File " + path + " should contain t, num_eqn, nstates, num_aux, num_dim"
        raise
        
    return t,num_eqn,nstates,num_aux,num_dim
