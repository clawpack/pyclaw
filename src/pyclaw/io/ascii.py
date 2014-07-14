#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing an ascii output file
"""

import os,sys
import logging
import numpy as np
import pickle

from ..util import read_data_line

logger = logging.getLogger('pyclaw.io')

def write(solution,frame,path,file_prefix='fort',write_aux=False,
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
        
        for state in solution.states:
            patch = state.patch

            write_patch_header(q_file,patch)
            
            if write_p:
                q = state.p
            else:
                q = state.q

            write_array(q_file, patch, q)

            if state.num_aux > 0 and write_aux:
                write_patch_header(aux_file,state.patch)
                write_array(aux_file,patch,state.aux)

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

    pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    pickle_file = open(pickle_filename,'wb')
    sol_dict = {'t':solution.t,'num_eqn':solution.num_eqn,'nstates':len(solution.states),
                     'num_aux':solution.num_aux,'num_dim':solution.domain.num_dim,
                     'write_aux':write_aux,
                     'problem_data' : solution.problem_data,
                     'mapc2p': solution.state.grid.mapc2p}
    if write_p:
        sol_dict['num_eqn'] = solution.mp

    pickle.dump(sol_dict, pickle_file)
    pickle_file.close()


def write_patch_header(f,patch):
    f.write("%5i                  patch_number\n" % patch.patch_index)
    f.write("%5i                  AMR_level\n" % patch.level)
    for dim in patch.dimensions:
        f.write("%5i                  m%s\n" % (dim.num_cells,dim.name))
    for dim in patch.dimensions:
        f.write("%18.8e     %slow\n" % (dim.lower,dim.name))
    for dim in patch.dimensions:
        f.write("%18.8e     d%s\n" % (dim.delta,dim.name))
    
    f.write("\n")


def write_array(f,patch,q):
    """
    Write a single array to output file f as ASCII text.

    The variable q here may in fact refer to q or to aux.
    """

    dims = patch.dimensions
    if patch.num_dim == 1:
        for k in xrange(dims[0].num_cells):
            for m in xrange(q.shape[0]):
                f.write("%18.8e" % q[m,k])
            f.write('\n')
    elif patch.num_dim == 2:
        for j in xrange(dims[1].num_cells):
            for k in xrange(dims[0].num_cells):
                for m in xrange(q.shape[0]):
                    f.write("%18.8e" % q[m,k,j])
                f.write('\n')    
            f.write('\n')
    elif patch.num_dim == 3:
        for l in xrange(dims[2].num_cells):
            for j in xrange(dims[1].num_cells):
                for k in xrange(dims[0].num_cells):
                    for m in xrange(q.shape[0]):
                        f.write("%18.8e" % q[m,k,j,l])
                    f.write('\n')
                f.write('\n')    
            f.write('\n')
    else:
        raise Exception("Dimension Exception in writing fort file.")


def read(solution,frame,path='./',file_prefix='fort',read_aux=False,
                options={}):
    r"""
    Read in a frame of ascii formatted files, and enter the data into the
    solution object.
    
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

    pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    problem_data = None
    mapc2p = None
    try:
        if os.path.exists(pickle_filename):
            pickle_file = open(pickle_filename,'rb')
            value_dict = pickle.load(pickle_file)
            problem_data = value_dict.get('problem_data',None)
            mapc2p       = value_dict.get('mapc2p',None)
    except IOError:
        logger.info("Unable to open pickle file %s" % (pickle_filename))


    # Construct path names
    base_path = os.path.join(path,)
    q_fname = os.path.join(base_path, '%s.q' % file_prefix) + str(frame).zfill(4)

    # Read in values from fort.t file:
    [t,num_eqn,nstates,num_aux,num_dim] = read_t(frame,path,file_prefix)

    patches = []
    
    # Read in values from fort.q file:
    try:
        f = open(q_fname,'r')
    except IOError:
        print "Error: file " + q_fname + " does not exist or is unreadable."
        raise
    
    n = np.zeros((num_dim))
    d = np.zeros((num_dim))
    lower = np.zeros((num_dim))

    # Loop through every patch setting the appropriate information
    for m in xrange(nstates):
    
        # Read in base header for this patch
        patch_index = read_data_line(f,data_type=int)
        level       = read_data_line(f,data_type=int)
        for i in xrange(num_dim):
            n[i] = read_data_line(f,data_type=int)
        for i in xrange(num_dim):
            lower[i] = read_data_line(f)
        for i in xrange(num_dim):
            d[i] = read_data_line(f)
    
        blank = f.readline()
    
        # Construct the patch
        # Since we do not have names here, we will construct the patch with
        # dimension names x,y,z
        names = ['x','y','z']
        import clawpack.pyclaw as pyclaw
        Dim = pyclaw.Dimension
        dimensions = [Dim(names[i],lower[i],lower[i] + n[i]*d[i],n[i]) for i in xrange(num_dim)]
        patch = pyclaw.geometry.Patch(dimensions)
        state= pyclaw.state.State(patch,num_eqn,num_aux)
        state.t = t
        state.problem_data = problem_data
        if mapc2p is not None:
            # If no mapc2p the default in the identity map in grid will be used
            state.grid.mapc2p = mapc2p

        if num_aux > 0:   
            state.aux[:]=0.
        
        # Fill in q values
        state.q = read_array(f, state, num_eqn)

        # Add AMR attributes:
        patch.patch_index = patch_index
        patch.level = level

        # Add new patch to solution
        solution.states.append(state)
        patches.append(state.patch)
    solution.domain = pyclaw.geometry.Domain(patches)
    f.close()

    # Read auxillary file if available and requested
    # Matching dimension parameter tolerances
    ABS_TOL = 1e-8
    REL_TOL = 1e-15
    if solution.states[0].num_aux > 0 and read_aux:
        # Check for aux file
        fname1 = os.path.join(base_path,'%s.a' % file_prefix)+str(frame).zfill(4)
        fname2 = os.path.join(base_path,'%s.a' % file_prefix)+str(0).zfill(4)
        if os.path.exists(fname1):
            # aux file specific to this frame:
            fname = fname1
        elif os.path.exists(fname2):
            # Assume that aux data from initial time is valid for all frames:
            fname = fname2
            # Note that this is generally not true when AMR is used.
            # Should give a better warning message in line below where
            # IOError exception is raised.
        else:
            logger.debug("Unable to open auxillary file %s or %s" % (fname1,fname2))
            return
            
        # Found a valid path, try to open and read it
        try:
            f = open(fname,'r')
        except IOError:
            logger.error("File %s was not able to be read." % fname)
            raise
            
        # Read in aux file
        for state in solution.states:
            patch = state.patch
            # Fetch correct patch
            patch_index = read_data_line(f,data_type=int)
    
            # These should match this patch already, raise exception otherwise
            if not (patch.level == read_data_line(f,data_type=int)):
                raise IOError("Patch level in aux file header did not match patch no %s." % patch.patch_index)
            for dim in patch.dimensions:
                num_cells = read_data_line(f,data_type=int)
                if not dim.num_cells == num_cells:
                    raise Exception("Dimension %s's num_cells in aux file header did not match patch no %s." % (dim.name,patch.patch_index))
            for dim in patch.dimensions:
                lower = read_data_line(f,data_type=float)
                if np.abs(lower - dim.lower) > ABS_TOL + REL_TOL * np.abs(dim.lower):
                    raise Exception('Value of lower in aux file does not match.')
            for dim in patch.dimensions:
                delta = read_data_line(f,data_type=float)
                if np.abs(delta - dim.delta) > ABS_TOL + REL_TOL * np.abs(dim.delta):
                    raise Exception('Value of delta in aux file does not match.')

            blank = f.readline()
    
            state.aux = read_array(f, state, num_aux)

        f.close()
        
            
def read_t(frame,path='./',file_prefix='fort'):
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
    logger.debug("Opening %s file." % path)
    try:
        f = open(path,'r')
    except(IOError):
        print "Error: file " + path + " does not exist or is unreadable."
        raise
        
    t = read_data_line(f)
    num_eqn = read_data_line(f,data_type=int)
    nstates = read_data_line(f,data_type=int)
    num_aux = read_data_line(f,data_type=int)
    num_dim = read_data_line(f,data_type=int)
    
    f.close()
        
    return t,num_eqn,nstates,num_aux,num_dim


def read_array(f, state, num_var):
    """
    Read in an array from an ASCII output file f.

    The variable q here may in fact refer to q or to aux.

    This routine supports the possibility that the values
    q[:,i,j,k] (for a fixed i,j,k) have been split over multiple lines, because
    some packages write just 4 values per line.
    For Clawpack 6.0, we plan to make all packages write
    q[:,i,j,k] on a single line.  This routine can then be simplified.
    """
    patch = state.patch
    q_shape = [num_var] + patch.num_cells_global
    q = np.zeros(q_shape)


    if patch.num_dim == 1:
        for i in xrange(patch.dimensions[0].num_cells):
            l = []
            while len(l)<num_var:
                line = f.readline()
                l = l + line.split()
            for m in xrange(num_var):
                q[m,i] = float(l[m])
    elif patch.num_dim == 2:
        for j in xrange(patch.dimensions[1].num_cells):
            for i in xrange(patch.dimensions[0].num_cells):
                l = []
                while len(l)<num_var:
                    line = f.readline()
                    l = l + line.split()
                for m in xrange(num_var):
                    q[m,i,j] = float(l[m])
            blank = f.readline()
    elif patch.num_dim == 3:
        for k in xrange(patch.dimensions[2].num_cells):
            for j in xrange(patch.dimensions[1].num_cells):
                for i in xrange(patch.dimensions[0].num_cells):
                    l=[]
                    while len(l) < num_var:
                        line = f.readline()
                        l = l + line.split()
                    for m in xrange(num_var):
                        q[m,i,j,k] = float(l[m])
                blank = f.readline()
            blank = f.readline()
    else:
        msg = "Read only supported up to 3d."
        logger.critical(msg)
        raise Exception(msg)

    return q

