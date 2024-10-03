#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing a petsc-style output file.

These routines preserve petclaw/pyclaw syntax for i/o while taking advantage of
PETSc's parallel i/o capabilities to allow for parallel reads and writes of
frame data.
"""
from petsc4py import PETSc
import pickle
import os
from six.moves import range
    

def write(solution,frame,path='./',file_prefix='claw',write_aux=False,
          options={},write_p=False):
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
    # Use defaults for options not passed in 
    clobber     = options.get('clobber', True)
    file_format = options.get('format', 'binary')

    if solution.num_aux == 0:
        write_aux = False

    filenames = set_filenames(frame,path,file_format,file_prefix,write_aux)
       
    if not clobber:
        for f in filenames.values():
            if os.path.exists(f):
                raise IOError('Cowardly refusing to clobber %s!' % f)

    rank =  PETSc.Comm.getRank(PETSc.COMM_WORLD)
    if rank==0:
        # Write metadata
        with open(filenames['metadata'],'wb') as metadata_file:
            # explicitly dumping a dictionary here to help out anybody trying to read the pickle file
            sol_dict = {'t':solution.t,'num_eqn':solution.num_eqn,'nstates':len(solution.states),
                             'num_aux':solution.num_aux,'num_dim':solution.domain.num_dim,
                             'write_aux':write_aux,
                             'problem_data' : solution.problem_data,
                             'mapc2p': solution.state.grid.mapc2p,
                             'file_format':file_format}
            if write_p:
                sol_dict['num_eqn'] = solution.mp

            pickle.dump(sol_dict, metadata_file)

            for state in solution.states:
                patch = state.patch
                pickle.dump({'level':patch.level,'names':patch.name,
                             'lower':patch.lower_global,
                             'num_cells':patch.num_cells_global,
                             'delta':patch.delta}, metadata_file)


    # Now write data
    q_viewer = set_up_viewers(filenames['q'],file_format.lower(),PETSc.Viewer.Mode.WRITE)
    if write_aux:
        aux_viewer = set_up_viewers(filenames['aux'],file_format.lower(),PETSc.Viewer.Mode.WRITE)
    
    for state in solution.states:
        patch = state.patch
        # we will reenable this bad boy when we switch over to petsc-dev:
        # state.q_da.view(q_viewer)
        if write_p:
            state.gpVec.view(q_viewer)
        else:
            state.gqVec.view(q_viewer)
        
        if write_aux:
            state.gauxVec.view(aux_viewer)
    
    q_viewer.flush()
    if write_aux:
        aux_viewer.flush()
    q_viewer.destroy() # Destroys aux_viewer also


def read(solution,frame,path='./',file_prefix='claw',read_aux=False,options={}):
    r"""
    Read in pickles and PETSc data files representing the solution
    
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Solution object to 
       read the data into.
     - *frame* - (int) Frame number to be read in
     - *path* - (string) Path to the current directory of the file
     - *file_prefix* - (string) Prefix of the files to be read in.  
       ``default = 'fort'``
     - *read_aux* (bool) Whether or not an auxiliary file will try to be read 
       in.  ``default = False``
     - *options* - (dict) Optional argument dictionary, see 
       `PETScIO Option Table`_
    
    .. _`PETScIO Option Table`:
    
    format   : one of 'ascii' or 'binary'
     
    """
    # Default format is binary
    file_format = options.get('format', 'binary')

    filenames = set_filenames(frame,path,file_format,file_prefix,read_aux)

    if read_aux:
        if not os.path.exists(filenames['aux']):
            # If no aux file for this frame, assume it is time-independent
            filenames['aux'] = os.path.join(path, '%s_aux.ptc' % file_prefix) + str(0).zfill(4)

    try:
        metadata_file = open(filenames['metadata'],'rb')
    except IOError:
        print("Error: file " + filenames['metadata'] + " does not exist or is unreadable.")
        raise

    # this dictionary is mostly holding debugging information, only nstates is needed
    # most of this information is explicitly saved in the individual patches
    value_dict = pickle.load(metadata_file)
    nstates    = value_dict['nstates']                    
    num_dim       = value_dict['num_dim']
    num_aux       = value_dict['num_aux']
    num_eqn       = value_dict['num_eqn']

    if read_aux and not os.path.exists(filenames['aux']):
        # Don't abort if aux file is missing
        from warnings import warn
        aux_file_path = os.path.join(path,filenames['aux'])
        warn('read_aux=True but aux file %s does not exist' % aux_file_path)
        read_aux = False

    q_viewer = set_up_viewers(filenames['q'],file_format.lower(),PETSc.Viewer.Mode.READ)
    if read_aux:
        aux_viewer = set_up_viewers(filenames['aux'],file_format.lower(),PETSc.Viewer.Mode.READ)

    patches = []
    for m in range(nstates):
        patch_dict = pickle.load(metadata_file)

        level   = patch_dict['level']
        names   = patch_dict['names']
        lower   = patch_dict['lower']
        n       = patch_dict['num_cells']
        d       = patch_dict['delta']

        from clawpack import petclaw
        dimensions = []
        for i in range(num_dim):
            dimensions.append(
                petclaw.Dimension(lower[i],lower[i] + n[i]*d[i],n[i],name=names[i]))
        patch = petclaw.Patch(dimensions)
        patch.level = level 
        state = petclaw.State(patch,num_eqn,num_aux)
        state.t = value_dict['t']
        state.problem_data = value_dict.get('problem_data',{})
        if 'mapc2p' in value_dict:
            # If no mapc2p is provided, leave the default identity map in grid
            state.grid.mapc2p = value_dict['mapc2p']

#       DA View/Load is broken in Petsc-3.1.8, we can load/view the DA if needed in petsc-3.2
#       state.q_da.load(q_viewer)
        state.gqVec.load(q_viewer)
        
        if read_aux:
            state.gauxVec.load(aux_viewer)
        
        solution.states.append(state)
        patches.append(state.patch)
    solution.domain = petclaw.geometry.Domain(patches)

    metadata_file.close()
    q_viewer.destroy() # Destroys aux_viewer also


def read_t(frame,path='./',file_prefix='claw'):
    r"""Read only the petsc.pkl file and return the data
    
    :Input:
     - *frame* - (int) Frame number to be read in
     - *path* - (string) Path to the current directory of the file
     - *file_prefix* - (string) Prefix of the files to be read in.  
       ``default = 'claw'``
     
    :Output:
     - (list) List of output variables
      - *t* - (int) Time of frame
      - *num_eqn* - (int) Number of equations in the frame
      - *npatches* - (int) Number of patches
      - *num_aux* - (int) Auxiliary value in the frame
      - *num_dim* - (int) Number of dimensions in q and aux
    
    """
    import logging
    logger = logging.getLogger('io')

    base_path = os.path.join(path,)
    path = os.path.join(base_path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    try:
        f = open(path,'rb')
    except IOError:
        print("Error: file " + path + " does not exist or is unreadable.")
        raise
    logger.debug("Opening %s file." % path)
    patch_dict = pickle.load(f)

    t      = patch_dict['t']
    num_eqn   = patch_dict['num_eqn']
    nstates = patch_dict['nstates']                    
    num_aux   = patch_dict['num_aux']                    
    num_dim   = patch_dict['num_dim']

    f.close()
        
    return t,num_eqn,nstates,num_aux,num_dim


def set_up_viewers(filename,file_format,mode):
    v = PETSc.Viewer()
    opts = {}
    if file_format == 'ascii':
        create_viewer = v.createASCII
    elif file_format == 'vtk':
        create_viewer = v.createASCII
        opts['format'] = PETSc.Viewer.Format.ASCII_VTK
    elif file_format == 'hdf5':
        create_viewer = v.createHDF5
    elif file_format == 'netcdf':
        create_viewer = v.createNetCDF
    elif file_format == 'binary':
        if hasattr(PETSc.Viewer,'createMPIIO'):
            create_viewer = v.createMPIIO
        else:
            create_viewer = v.createBinary
    else:
        raise IOError('PETSc has no viewer for the output format %s ' % file_format)

    viewer = create_viewer(filename, mode, **opts)
    return viewer


def set_filenames(frame,path,file_format,file_prefix,do_aux):
    filenames = {}
    filenames['metadata'] = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    if file_format == 'vtk':
        filenames['q'] = os.path.join(path, file_prefix+str(frame).zfill(4)+'.vtk')
    else:
        filenames['q'] = os.path.join(path, '%s.ptc' % file_prefix) + str(frame).zfill(4)

    if do_aux:
        filenames['aux'] = os.path.join(path, '%s_aux.ptc' % file_prefix) + str(frame).zfill(4)

    return filenames
