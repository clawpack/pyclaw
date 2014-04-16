#!/usr/bin/env python
# encoding: utf-8
r"""
Routines for reading and writing multifabs.

These routines preserve petclaw/pyclaw syntax for i/o while taking advantage of
PETSc's parallel i/o capabilities to allow for parallel reads and writes of
frame data.
"""

import boxlib
import pickle
import numpy as np

def write(solution,frame,path='./',file_prefix='claw',write_aux=False,
          options={},write_p=False):
    r"""
        Write out pickle and multifab data files representing the
        solution.  Common data is written from process 0 in pickle
        files.  Shared data is written from all processes into BoxLib
        plot files.

    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) pyclaw
       object to be output
     - *frame* - (int) Frame number
     - *path* - (string) Root path
     - *file_prefix* - (string) Prefix for the file name. ``default =
        'claw'``
     - *write_aux* - (bool) Boolean controlling whether the associated
       auxiliary array should be written out. ``default = False``
     - *options* - (dict) Optional argument dictionary.
    """

    import os

    # Option parsing
    # option = {'format':'binary','clobber':True}.update(options)
    # clobber = options['clobber']

    clobber = True

    pickle_filename = os.path.join(path, '%s.pkl' % file_prefix) + str(frame).zfill(4)
    plt_filename = os.path.join(path, '%s.plt' % file_prefix) + str(frame).zfill(4)

    # if solution.num_aux == 0:
    #     write_aux = False
    # if write_aux:
    #     aux_filename = os.path.join(path, '%s_aux.ptc' % file_prefix) + str(frame).zfill(4)

    if not clobber:
        for f in (pickle_filename, viewer_filename, aux_filename):
            if os.path.exists(f):
                raise IOError('Cowardly refusing to clobber %s!' % f)

    rank = boxlib.rank()
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

    for state in solution.states:
        patch = state.patch
        if rank==0:
            pickle.dump({'level':patch.level,
                         'names':patch.name,'lower':patch.lower_global,
                         'num_cells':patch.num_cells_global,'delta':patch.delta}, pickle_file)

        mf = solution.state._create_multifab(solution.num_eqn, 0)
        q  = solution.state._get_array(mf)
        # print q.shape
        # print solution.state.q.shape
        # print np.rollaxis(solution.state.q, 0, q.ndim).shape
        q[...] = np.rollaxis(solution.state.q, 0, q.ndim)

        try:
            if not os.path.exists(plt_filename):
                os.mkdir(plt_filename)
            if not os.path.exists(plt_filename + '/Level_0'):
                os.mkdir(plt_filename + '/Level_0')
        except:
            pass

        # build multifab
        mf = solution.state._create_multifab(solution.state.num_eqn, 0)
        solution.state._copy_into_multifab(mf, 0, solution.state.q)

        # write header
        if rank==0:
            with open(os.path.join(plt_filename, 'Header'), 'w') as f:
                write = lambda s: f.write(str(s) + '\n')

                plo = solution.domain.patch.lower_global
                phi = solution.domain.patch.upper_global
                dx  = solution.domain.patch.delta

                write("HyperCLaw-V1.1")                 # yeah, i have no idea why...
                write(solution.state.num_eqn)           # number of variables
                for i in range(solution.state.num_eqn): # variable names
                    write('q'+str(i))
                write(solution.state.num_dim)           # number of dimensions
                write(solution.state.t)                 # time
                write('0')                              # finest amr level number
                write(' '.join(map(str, plo)))          # lower corner of problem domain
                write(' '.join(map(str, phi)))          # upper corner of problem domain
                write('1')                              # refinement ratio
                write(solution.state.patch._gbox)       # grid domain
                write(frame)                            # time step number
                write(' '.join(map(str, dx)))           # dx
                write('0')                              # cartesian coords
                write('0')                              # boundary data? (0=no)

                # now write info for each amr level (only one in this case)
                write('0 %d %f' % (mf.size(), solution.state.t))
                write(frame)
                for i in range(mf.size()):
                    bx = solution.state.patch._ba.get(i)
                    lo = bx.smallEnd()
                    hi = bx.bigEnd()
                    for d in range(solution.state.num_dim):
                        write("%lf %lf" % (plo[d] + lo[d]*dx[d], plo[d] + hi[d]*dx[d]))
                write("Level_0/Cell")

        # write cell data
        mf.writeOut(plt_filename + '/Level_0/Cell')

    if rank==0:
        pickle_file.close()


