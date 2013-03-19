#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing a petsc-style output file.

These routines preserve petclaw/pyclaw syntax for i/o while taking advantage of PETSc's parallel i/o capabilities to allow for parallel reads and writes of frame data.
"""

import os
import logging
logger = logging.getLogger('io')

from petsc4py import PETSc
import pickle
    

def write_petsc(solution,frame,path='./',file_prefix='claw',write_aux=False,options={},write_p=False):
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

    for k in option_defaults.iterkeys():
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

    if solution.num_aux > 0 and write_aux:
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
        if write_p:
            pickle.dump({'t':solution.t,'num_eqn':solution.mp,'nstates':len(solution.states),
                         'num_aux':solution.num_aux,'num_dim':solution.domain.num_dim,'write_aux':write_aux,
                         'problem_data' : solution.problem_data}, pickle_file)
        else:
            pickle.dump({'t':solution.t,'num_eqn':solution.num_eqn,'nstates':len(solution.states),
                         'num_aux':solution.num_aux,'num_dim':solution.domain.num_dim,'write_aux':write_aux,
                         'problem_data' : solution.problem_data}, pickle_file)

    # now set up the PETSc viewers
    if options['format'] == 'ascii':
        viewer = PETSc.Viewer().createASCII(viewer_filename, PETSc.Viewer.Mode.WRITE)
        if write_aux:
            aux_viewer = PETSc.Viewer().createASCII(aux_filename, PETSc.Viewer.Mode.WRITE) 
    elif options['format'] == 'binary':
        if hasattr(PETSc.Viewer,'createMPIIO'):
            viewer = PETSc.Viewer().createMPIIO(viewer_filename, PETSc.Viewer.Mode.WRITE)
        else:
            viewer = PETSc.Viewer().createBinary(viewer_filename, PETSc.Viewer.Mode.WRITE)
        if write_aux:
            if hasattr(PETSc.Viewer,'createMPIIO'):
                aux_viewer = PETSc.Viewer().createMPIIO(aux_filename, PETSc.Viewer.Mode.WRITE)
            else:
                aux_viewer = PETSc.Viewer().createBinary(aux_filename, PETSc.Viewer.Mode.WRITE)
    elif options['format'] == 'vtk':
        viewer = PETSc.Viewer().createASCII(viewer_filename, PETSc.Viewer.Mode.WRITE, format=PETSc.Viewer.Format.ASCII_VTK)
        if write_aux:
            aux_viewer = PETSc.Viewer().createASCII(aux_filename, PETSc.Viewer.Mode.WRITE) 
    else:
        raise IOError('format type %s not supported' % options['format'])
    
    for state in solution.states:
        patch = state.patch
        if rank==0:
            pickle.dump({'level':patch.level,
                         'names':patch.name,'lower':patch.lower_global,
                         'num_cells':patch.num_cells_global,'delta':patch.delta}, pickle_file)
#       we will reenable this bad boy when we switch over to petsc-dev
#        state.q_da.view(viewer)
        if write_p:
            state.gpVec.view(viewer)
        else:
            state.gqVec.view(viewer)
        
        if write_aux:
            state.gauxVec.view(aux_viewer)
    
    viewer.flush()
    viewer.destroy()
    if rank==0:
        pickle_file.close()
    if write_aux:
        aux_viewer.flush()
        aux_viewer.destroy()

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
     - *read_aux* (bool) Whether or not an auxiliary file will try to be read 
       in.  ``default = False``
     - *options* - (dict) Optional argument dictionary, see 
       `PETScIO Option Table`_
    
    .. _`PETScIO Option Table`:
    
    format   : one of 'ascii' or 'binary'
     
    """

    # Option parsing
    option_defaults = {'format':'binary'}

    for k in option_defaults.iterkeys():
        if options.has_key(k):
            pass
        else:
            options[k] = option_defaults[k]
    
    pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    viewer_filename = os.path.join(path, '%s.ptc' % file_prefix) + str(frame).zfill(4)
    aux_viewer_filename1 = os.path.join(path, '%s_aux.ptc' % file_prefix) + str(frame).zfill(4)
    aux_viewer_filename2 = os.path.join(path, '%s_aux.ptc' % file_prefix) + str(0).zfill(4)
    if os.path.exists(aux_viewer_filename1):
         aux_viewer_filename = aux_viewer_filename1
    else:
         aux_viewer_filename = aux_viewer_filename2


    if frame < 0:
        # Don't construct file names with negative frameno values.
        raise IOError("Frame " + str(frame) + " does not exist ***")

    pickle_file = open(pickle_filename,'rb')

    # this dictionary is mostly holding debugging information, only nstates is needed
    # most of this information is explicitly saved in the individual patches
    value_dict = pickle.load(pickle_file)
    nstates    = value_dict['nstates']                    
    num_dim       = value_dict['num_dim']
    num_aux       = value_dict['num_aux']
    num_eqn       = value_dict['num_eqn']

    # now set up the PETSc viewer
    if options['format'] == 'ascii':
        viewer = PETSc.Viewer().createASCII(viewer_filename, PETSc.Viewer.Mode.READ)
        if read_aux:
            aux_viewer = PETSc.Viewer().createASCII(aux_viewer_filename, PETSc.Viewer.Mode.READ)
    elif options['format'] == 'binary':
        if hasattr(PETSc.Viewer,'createMPIIO'):
            viewer = PETSc.Viewer().createMPIIO(viewer_filename, PETSc.Viewer.Mode.READ)
        else:
            viewer = PETSc.Viewer().createBinary(viewer_filename, PETSc.Viewer.Mode.READ)
        if read_aux:
            if os.path.exists(aux_viewer_filename):
                if hasattr(PETSc.Viewer,'createMPIIO'):
                    aux_viewer = PETSc.Viewer().createMPIIO(aux_viewer_filename, PETSc.Viewer.Mode.READ)
                else:
                    aux_viewer = PETSc.Viewer().createBinary(aux_viewer_filename, PETSc.Viewer.Mode.READ)
            else:
                from warnings import warn
                aux_file_path = os.path.join(path,aux_viewer_filename)
                warn('read_aux=True but aux file %s does not exist' % aux_file_path)
                read_aux=False
    else:
        raise IOError('format type %s not supported' % options['format'])

    patches = []
    for m in xrange(nstates):
        patch_dict = pickle.load(pickle_file)

        level   = patch_dict['level']
        names   = patch_dict['names']
        lower   = patch_dict['lower']
        n       = patch_dict['num_cells']
        d       = patch_dict['delta']

        from clawpack import petclaw
        dimensions = []
        for i in xrange(num_dim):
            dimensions.append(
                #pyclaw.solution.Dimension(names[i],lower[i],lower[i] + n[i]*d[i],n[i]))
                petclaw.Dimension(names[i],lower[i],lower[i] + n[i]*d[i],n[i]))
        #patch = pyclaw.solution.Patch(dimensions)
        patch = petclaw.Patch(dimensions)
        patch.level = level 
        #state = pyclaw.state.State(patch)
        state = petclaw.State(patch,num_eqn,num_aux) ##
        state.t = value_dict['t']
        state.problem_data = value_dict['problem_data']

#       DA View/Load is broken in Petsc-3.1.8, we can load/view the DA if needed in petsc-3.2
#       state.q_da.load(viewer)
        state.gqVec.load(viewer)
        
        if read_aux:
            state.gauxVec.load(aux_viewer)
        
        solution.states.append(state)
        patches.append(state.patch)
    solution.domain = petclaw.geometry.Domain(patches)

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
      - *num_eqn* - (int) Number of equations in the frame
      - *npatches* - (int) Number of patches
      - *num_aux* - (int) Auxillary value in the frame
      - *num_dim* - (int) Number of dimensions in q and aux
    
    """

    base_path = os.path.join(path,)
    path = os.path.join(base_path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    try:
        f = open(path,'rb')
        logger.debug("Opening %s file." % path)
        patch_dict = pickle.load(f)

        t      = patch_dict['t']
        num_eqn   = patch_dict['num_eqn']
        nstates = patch_dict['nstates']                    
        num_aux   = patch_dict['num_aux']                    
        num_dim   = patch_dict['num_dim']

        f.close()
    except(IOError):
        raise
    except:
        logger.error("File " + path + " should contain t, num_eqn, npatches, num_aux, num_dim")
        print "File " + path + " should contain t, num_eqn, npatches, num_aux, num_dim"
        raise
        
    return t,num_eqn,nstates,num_aux,num_dim
