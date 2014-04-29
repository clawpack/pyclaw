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
    

def write(solution,frame,path='./',prefix='claw',file_format='binary',clobber=True,
          write_aux=False,write_p=False,options={}, **kwargs):
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
     - *prefix* - (string) Prefix for the file name. ``default = 'claw'``
     - *file_format* - (string) Format to output data, ``default = 'binary'``
     - *clobber* - (bool) Bollean controlling whether to overwrite files
     - *write_aux* - (bool) Boolean controlling whether the associated 
       auxiliary array should be written out. ``default = False``     
     - *write_p* - (bool) Boolean controlling whether the associated 
       p array should be written out. ``default = False``
     - *options* - (dict) Optional argument dictionary     
    
    : Available formats
     - 'ascii', 'binary', 'hdf5', 'netcdf', 'vtk'

    """    
    import os

    if 'format' in kwargs:
        file_format = kwargs['file_format']
    if 'file_prefix' in kwargs:
        prefix = kwargs['file_prefix']

    pickle_filename = os.path.join(path, '%s.pkl' % prefix) + str(frame).zfill(4)
    
    if file_format == 'vtk':
        viewer_filename = os.path.join(path, prefix+str(frame).zfill(4)+'.vtk')
    else:
        viewer_filename = os.path.join(path, '%s.ptc' % prefix) + str(frame).zfill(4)

    if solution.num_aux == 0:
        write_aux = False
    if write_aux:
        aux_filename = os.path.join(path, '%s_aux.ptc' % prefix) + str(frame).zfill(4)

    if not clobber:
        for f in (pickle_filename, viewer_filename, aux_filename):
            if os.path.exists(f):
                raise IOError('Cowardly refusing to clobber %s!' % f)

    rank =  PETSc.Comm.getRank(PETSc.COMM_WORLD)
    if rank==0:
        pickle_file = open(pickle_filename,'wb')
        # explicitly dumping a dictionary here to help out anybody trying to read the pickle file
        sol_dict = {'t':solution.t,'num_eqn':solution.num_eqn,'nstates':len(solution.states),
                         'num_aux':solution.num_aux,'num_dim':solution.domain.num_dim,
                         'write_aux':write_aux,
                         'problem_data' : solution.problem_data,
                         'mapc2p': solution.state.grid.mapc2p}
        if write_p:
            sol_dict['num_eqn'] = solution.mp

        pickle.dump(sol_dict, pickle_file)

    # now set up the PETSc viewers
    if file_format == 'ascii':
        viewer = PETSc.Viewer().createASCII(viewer_filename, PETSc.Viewer.Mode.WRITE)
        if write_aux:
            aux_viewer = PETSc.Viewer().createASCII(aux_filename, PETSc.Viewer.Mode.WRITE)
    elif file_format == 'binary':
        if hasattr(PETSc.Viewer,'createMPIIO'):
            viewer = PETSc.Viewer().createMPIIO(viewer_filename, PETSc.Viewer.Mode.WRITE)
        else:
            viewer = PETSc.Viewer().createBinary(viewer_filename, PETSc.Viewer.Mode.WRITE)
        if write_aux:
            if hasattr(PETSc.Viewer,'createMPIIO'):
                aux_viewer = PETSc.Viewer().createMPIIO(aux_filename, PETSc.Viewer.Mode.WRITE)
            else:
                aux_viewer = PETSc.Viewer().createBinary(aux_filename, PETSc.Viewer.Mode.WRITE)
    elif file_format == 'vtk':
        viewer = PETSc.Viewer().createASCII(viewer_filename, PETSc.Viewer.Mode.WRITE, format=PETSc.Viewer.Format.ASCII_VTK)
        if write_aux:
            aux_viewer = PETSc.Viewer().createASCII(aux_filename, PETSc.Viewer.Mode.WRITE)
    elif file_format=='hdf5':
        if hasattr(PETSc.Viewer, 'createHDF5'):
            viewer = PETSc.Viewer().createHDF5(viewer_filename, PETSc.Viewer.Mode.WRITE)
        else:
            viewer = PETSc.Viewer().createBinary(viewer_filename, PETSc.Viewer.Mode.WRITE)
        if write_aux:
            if hasattr(PETSc.Viewer,'createHDF5'):
                aux_viewer = PETSc.Viewer().createHDF5(aux_filename, PETSc.Viewer.Mode.WRITE)
            else:
                aux_viewer = PETSc.Viewer().createBinary(aux_filename, PETSc.Viewer.Mode.WRITE)
    elif file_format=='netcdf':
        if hasattr(PETSc.Viewer, 'createNetCDF'):
            viewer = PETSc.Viewer().createHDF5(viewer_filename, PETSc.Viewer.Mode.WRITE)
        else:
            viewer = PETSc.Viewer().createBinary(viewer_filename, PETSc.Viewer.Mode.WRITE)
            msg = 'PETSc.Viewer has no attribute createHDF5, falling back to binary output'
            petsc_warning=True
            file_format='binary'
        if write_aux:
            if hasattr(PETSc.Viewer,'createHDF5'):
                aux_viewer = PETSc.Viewer().createHDF5(aux_filename, PETSc.Viewer.Mode.WRITE)
            else:
                aux_viewer = PETSc.Viewer().createBinary(aux_filename, PETSc.Viewer.Mode.WRITE)
    else:
        raise IOError('file format type %s not supported' % file_format)
    

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


def read(solution,frame,path='./',prefix='claw',file_format='binary',read_aux=False,options={},**kwargs):
    r"""
    Read in pickles and PETSc data files representing the solution
    
    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Solution object to 
       read the data into.
     - *frame* - (int) Frame number to be read in
     - *path* - (string) Path to the current directory of the file
     - *prefix* - (string) Prefix of the files to be read in.  
       ``default = 'fort'``
     - *read_aux* (bool) Whether or not an auxiliary file will try to be read 
       in.  ``default = False``
     - *options* - (dict) Optional argument dictionary, see 
       `PETScIO Option Table`_
    
    .. _`PETScIO Option Table`:
    
    format   : one of 'ascii' or 'binary'
     
    """
    import os

    if 'format' in kwargs:
        file_format = kwargs['file_format']
    if 'file_prefix' in kwargs:
        prefix = kwargs['file_prefix']

    pickle_filename = os.path.join(path, '%s.pkl' % prefix) + str(frame).zfill(4)
    viewer_filename = os.path.join(path, '%s.ptc' % prefix) + str(frame).zfill(4)
    aux_viewer_filename1 = os.path.join(path, '%s_aux.ptc' % prefix) + str(frame).zfill(4)
    aux_viewer_filename2 = os.path.join(path, '%s_aux.ptc' % prefix) + str(0).zfill(4)
    if os.path.exists(aux_viewer_filename1):
         aux_viewer_filename = aux_viewer_filename1
    else:
         aux_viewer_filename = aux_viewer_filename2

    try:
        pickle_file = open(pickle_filename,'rb')
    except IOError:
        print "Error: file " + pickle_filename + " does not exist or is unreadable."
        raise

    # this dictionary is mostly holding debugging information, only nstates is needed
    # most of this information is explicitly saved in the individual patches
    value_dict = pickle.load(pickle_file)
    nstates    = value_dict['nstates']                    
    num_dim       = value_dict['num_dim']
    num_aux       = value_dict['num_aux']
    num_eqn       = value_dict['num_eqn']

    # now set up the PETSc viewer
    if file_format == 'ascii':
        viewer = PETSc.Viewer().createASCII(viewer_filename, PETSc.Viewer.Mode.READ)
        if read_aux:
            aux_viewer = PETSc.Viewer().createASCII(aux_viewer_filename, PETSc.Viewer.Mode.READ)
    elif file_format == 'binary':
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
    elif file_format == 'hdf5':
        if hasattr(PETSc.Viewer,'createHDF5'):
            viewer = PETSc.Viewer().createHDF5(viewer_filename, PETSc.Viewer.Mode.READ)
        else:
            viewer = PETSc.Viewer().createBinary(viewer_filename, PETSc.Viewer.Mode.READ)
        if read_aux:
            if os.path.exists(aux_viewer_filename):
                if hasattr(PETSc.Viewer,'createHDF5'):
                    aux_viewer = PETSc.Viewer().createHDF5(aux_viewer_filename, PETSc.Viewer.Mode.READ)
                else:
                    aux_viewer = PETSc.Viewer().createBinary(aux_viewer_filename, PETSc.Viewer.Mode.READ)
            else:
                from warnings import warn
                aux_file_path = os.path.join(path,aux_viewer_filename)
                warn('read_aux=True but aux file %s does not exist' % aux_file_path)
                read_aux=False
    elif file_format == 'netcdf':
        if hasattr(PETSc.Viewer,'createNetCDF'):
            viewer = PETSc.Viewer().createNetCDF(viewer_filename, PETSc.Viewer.Mode.READ)
        else:
            viewer = PETSc.Viewer().createBinary(viewer_filename, PETSc.Viewer.Mode.READ)
        if read_aux:
            if os.path.exists(aux_viewer_filename):
                if hasattr(PETSc.Viewer,'createNetCDF'):
                    aux_viewer = PETSc.Viewer().createNetCDF(aux_viewer_filename, PETSc.Viewer.Mode.READ)
                else:
                    aux_viewer = PETSc.Viewer().createBinary(aux_viewer_filename, PETSc.Viewer.Mode.READ)
            else:
                from warnings import warn
                aux_file_path = os.path.join(path,aux_viewer_filename)
                warn('read_aux=True but aux file %s does not exist' % aux_file_path)
                read_aux=False
    else:
        raise IOError('file format type %s not supported' % file_format)

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
                petclaw.Dimension(names[i],lower[i],lower[i] + n[i]*d[i],n[i]))
        patch = petclaw.Patch(dimensions)
        patch.level = level 
        state = petclaw.State(patch,num_eqn,num_aux)
        state.t = value_dict['t']
        state.problem_data = value_dict.get('problem_data',{})
        if value_dict.has_key('mapc2p'):
            # If no mapc2p is provided, leave the default identity map in grid
            state.grid.mapc2p = value_dict['mapc2p']

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

def read_t(frame,path='./',prefix='claw',**kwargs):
    r"""Read only the petsc.pkl file and return the data
    
    :Input:
     - *frame* - (int) Frame number to be read in
     - *path* - (string) Path to the current directory of the file
     - *prefix* - (string) Prefix of the files to be read in.  
       ``default = 'claw'``
     
    :Output:
     - (list) List of output variables
      - *t* - (int) Time of frame
      - *num_eqn* - (int) Number of equations in the frame
      - *npatches* - (int) Number of patches
      - *num_aux* - (int) Auxillary value in the frame
      - *num_dim* - (int) Number of dimensions in q and aux
    
    """
    import os
    import logging
    logger = logging.getLogger('io')

    if 'file_prefix' in kwargs:
        prefix = kwargs['file_prefix']
    
    base_path = os.path.join(path,)
    path = os.path.join(base_path, '%s.pkl' % prefix) + str(frame).zfill(4)
    try:
        f = open(path,'rb')
    except IOError:
        print "Error: file " + path + " does not exist or is unreadable."
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
