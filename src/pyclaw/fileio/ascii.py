#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing an ascii output file
"""

import os
import logging
import numpy as np
import pickle
import clawpack.pyclaw as pyclaw

from clawpack.pyclaw.util import read_data_line

logger = logging.getLogger('pyclaw.fileio')

def write(solution, frame, path, file_prefix='fort', write_aux=False,
                    options={}, write_p=False):
    r"""
    Write out ascii data file
    
    Write out an ascii file formatted identical to the fortran clawpack files
    including writing out fort.t, fort.q, and fort.aux if necessary.  Note
    that there are some parameters that assumed to be the same for every patch
    in this format which is not necessarily true for the actual data objects.
    Make sure that if you use this output format that all of your patches share
    the appropriate values of num_dim, num_eqn, num_aux, and t.  Only supports
    up to 3 dimensions.
    
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
    # Write fort.txxxx file
    file_name = '%s.t%s' % (file_prefix,str(frame).zfill(4))
    with open(os.path.join(path,file_name),'w') as f:
        f.write("%18.8e     time\n" % solution.t)
        f.write("%5i                  num_eqn\n" % solution.num_eqn)
        f.write("%5i                  nstates\n" % len(solution.states))
        f.write("%5i                  num_aux\n" % solution.num_aux)
        f.write("%5i                  num_dim\n" % solution.domain.num_dim)
        f.write("%5i                  num_ghost\n" % 0)
        f.write("%s                  file_format\n" % "ascii")

    # Write fort.qxxxx file
    file_name = '%s.q%s' % (file_prefix,str(frame).zfill(4))
    with open(os.path.join(path,file_name),'w') as q_file:
        for state in solution.states:
            write_patch_header(q_file,state.patch)
            if write_p:
                q = state.p
            else:
                q = state.q
            write_array(q_file, state.patch, q)

    # Write fort.auxxxxx file if required
    if solution.num_aux > 0 and write_aux:
        file_name = '%s.a%s' % (file_prefix,str(frame).zfill(4))
        with open(os.path.join(path,file_name),'w') as aux_file:
            for state in solution.states:
                write_patch_header(aux_file,state.patch)
                write_array(aux_file,state.patch,state.aux)

    # Write fort.pklxxxx file
    sol_dict = {'t':solution.t,'num_eqn':solution.num_eqn,'nstates':len(solution.states),
                     'num_aux':solution.num_aux,'num_dim':solution.domain.num_dim,
                     'write_aux':write_aux,
                     'problem_data' : solution.problem_data,
                     'mapc2p': solution.state.grid.mapc2p}
    if write_p:
        sol_dict['num_eqn'] = solution.mp
    pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    with open(pickle_filename,'wb') as pickle_file:
        pickle.dump(sol_dict, pickle_file)


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
        for k in range(dims[0].num_cells):
            for m in range(q.shape[0]):
                f.write("%18.8e" % q[m,k])
            f.write('\n')
    elif patch.num_dim == 2:
        for j in range(dims[1].num_cells):
            for k in range(dims[0].num_cells):
                for m in range(q.shape[0]):
                    f.write("%18.8e" % q[m,k,j])
                f.write('\n')    
            f.write('\n')
    elif patch.num_dim == 3:
        for l in range(dims[2].num_cells):
            for j in range(dims[1].num_cells):
                for k in range(dims[0].num_cells):
                    for m in range(q.shape[0]):
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
     - *read_aux* (bool) Whether or not an auxiliary file will try to be read 
       in.  ``default = False``
     - *options* - (dict) Dictionary of optional arguments dependent on 
       the format being read in.  ``default = {}``
    """

    pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    problem_data = None
    mapc2p = None
    try:
        if os.path.exists(pickle_filename):
            with open(pickle_filename,'rb') as pickle_file:
                value_dict = pickle.load(pickle_file)
            problem_data = value_dict.get('problem_data',None)
            mapc2p       = value_dict.get('mapc2p',None)
    except IOError:
        logger.info("Unable to open pickle file %s" % (pickle_filename))


    # Construct path names
    base_path = os.path.join(path,)
    q_fname = os.path.join(base_path, '%s.q' % file_prefix) + str(frame).zfill(4)

    # Read in values from fort.t file:
    [t,num_eqn,nstates,num_aux,num_dim,num_ghost,file_format] = \
         read_t(frame,path,file_prefix)

    patches = []

    # Read in values from fort.q file:
    with open(q_fname,'r') as f:
        # Loop through every patch setting the appropriate information
        for m in range(nstates):
            # Read header for this patch
            patch = read_patch_header(f, num_dim)

            # Construct state
            state= pyclaw.state.State(patch,num_eqn,num_aux)
            state.t = t
            state.problem_data = problem_data
            if mapc2p is not None:
                # If no mapc2p the default identity map in grid will be used
                state.grid.mapc2p = mapc2p

            if num_aux > 0:
                # Write NaNs for now to indicate this is uninitialized
                state.aux[:] = np.nan

            # Fill in q values
            state.q = read_array(f, state, num_eqn)

            # Add new patch to solution
            solution.states.append(state)
            patches.append(state.patch)

    solution.domain = pyclaw.geometry.Domain(patches)

    # Read auxiliary file if available and requested
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
        else:
            logger.debug("Unable to open auxiliary file %s or %s" % (fname1,fname2))
            return
            
        # Read in fort.auxxxxx file
        with open(fname,'r') as f:
            for state in solution.states:
                patch = state.patch
                aux_patch = read_patch_header(f, num_dim)

                # Read patch header and check that it matches that from fort.qxxxx
                assert patch.level == aux_patch.level, \
                        "Patch level in aux file header did not match patch no %s." % patch.patch_index
                for i in range(len(patch.dimensions)):
                    assert patch.dimensions[i].num_cells == aux_patch.dimensions[i].num_cells, \
                        "Dimension %s's num_cells in aux file header did not match patch no %s." % (patch.dimensions[i].name, patch.patch_index)
                    assert np.abs(patch.dimensions[i].lower - aux_patch.dimensions[i].lower) <= ABS_TOL + REL_TOL * np.abs(patch.dimensions[i].lower), \
                            'Value of lower in aux file does not match.'
                    assert np.abs(patch.dimensions[i].delta - aux_patch.dimensions[i].delta) <= ABS_TOL + REL_TOL * np.abs(patch.dimensions[i].delta), \
                            'Value of delta in aux file does not match.'

                state.aux = read_array(f, state, num_aux)

        
def read_t(frame,path='./',file_prefix='fort'):
    r"""Read only the fort.t file and return the data.

    Note this file is always ascii and now contains a line that tells
    the file_format, so we can read this file before importing the 
    appropriate read function for the solution data.

    For backward compatibility, if file_format line is missing then
    return None and handle this where it is called.

    This version also reads in num_ghost so that if the data is binary,
    we can extract only the data that's relevant (since ghost cells are
    included).
    
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
     - *num_aux* - (int) Auxiliary value in the frame
     - *num_dim* - (int) Number of dimensions in q and aux
     - *num_ghost* - (int) Number of ghost cells on each side
     - *file_format* - (str) 'ascii', 'binary32', 'binary64'
    
    """

    from clawpack.pyclaw.util import read_data_line
    import logging
    logger = logging.getLogger('pyclaw.fileio')

    base_path = os.path.join(path,)
    path = os.path.join(base_path, '%s.t' % file_prefix) + str(frame).zfill(4)
    logger.debug("Opening %s file." % path)
    with open(path,'r') as f:
        t = read_data_line(f)
        num_eqn = read_data_line(f, data_type=int)
        nstates = read_data_line(f, data_type=int)
        num_aux = read_data_line(f, data_type=int)
        num_dim = read_data_line(f, data_type=int)
        try:
            num_ghost = read_data_line(f, data_type=int)
        except:
            num_ghost = 0
        try:
            file_format = read_data_line(f, data_type=str)
        except:
            file_format = None
        
    return t,num_eqn,nstates,num_aux,num_dim,num_ghost,file_format



def read_patch_header(f, num_dim):
    r"""Read header describing the next patch
    
    :Input:
     - *f* - (file) Handle to open file
     - *num_dim* - (int) Number of dimensions
     
    :Output:
     - *patch* - (clawpack.pyclaw.geometry.Patch) Initialized patch represented
       by the header data.
    
    """

    n = np.zeros((num_dim), dtype=int)
    d = np.zeros((num_dim))
    lower = np.zeros((num_dim))
    patch_index = read_data_line(f, data_type=int)
    level       = read_data_line(f, data_type=int)
    for i in range(num_dim):
        n[i] = read_data_line(f, data_type=int)
    for i in range(num_dim):
        lower[i] = read_data_line(f)
    for i in range(num_dim):
        d[i] = read_data_line(f)

    blank = f.readline()

    # Construct the patch
    # Since we do not have names here, we will construct each patch with
    # dimension names x,y,z
    names = ['x', 'y', 'z']
    dimensions = [pyclaw.Dimension(lower[i], lower[i] + n[i] * d[i],
                                  n[i], name=names[i]) for i in range(num_dim)]
    patch = pyclaw.geometry.Patch(dimensions)

    # Add AMR attributes:
    patch.patch_index = patch_index
    patch.level = level

    return patch


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
    q = np.zeros(q_shape, order='F')


    try:
        if patch.num_dim == 1:
            for i in range(patch.dimensions[0].num_cells):
                l = []
                while len(l)<num_var:
                    line = f.readline()
                    if line == '':
                        raise IOError('Unexpected EOF in %s' % f.name)
                    l = l + line.split()
                for m in range(num_var):
                    q[m,i] = float(l[m])
        elif patch.num_dim == 2:
            for j in range(patch.dimensions[1].num_cells):
                for i in range(patch.dimensions[0].num_cells):
                    l = []
                    while len(l) < num_var:
                        line = f.readline()
                        if line == '':
                            raise IOError('Unexpected EOF in %s' % f.name)
                        l = l + line.split()
                    for m in range(num_var):
                        q[m,i,j] = float(l[m])
                blank = f.readline()
        elif patch.num_dim == 3:
            for k in range(patch.dimensions[2].num_cells):
                for j in range(patch.dimensions[1].num_cells):
                    for i in range(patch.dimensions[0].num_cells):
                        l=[]
                        while len(l) < num_var:
                            line = f.readline()
                            if line == '':
                                raise IOError('Unexpected EOF in %s' % f.name)
                            l = l + line.split()
                        for m in range(num_var):
                            q[m,i,j,k] = float(l[m])
                    blank = f.readline()
                blank = f.readline()
        else:
            msg = "Read only supported up to 3d."
            logger.critical(msg)
            raise Exception(msg)
    except:
        msg = '*** Problem reading patch data'
        msg = msg + '\n*** Format might be binary, is plotdata.format set properly?'
        logger.critical(msg)
        raise IOError(msg)

    return q
