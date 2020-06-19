#!/usr/bin/python

# Author: Xinsheng (Shawn) Qin
# date: 03/2016
# Modified by Katy Barnhart
# date 05/2020

import sys
import os
import numpy as np
from clawpack.pyclaw import Solution
import vtk

def write(
    solution,
    frame,
    path="_output",
    file_prefix='claw',
    write_aux=None,
    options=None,
    write_p=None,
    ):
    """Write out a VTK representation of solution

    For each input frame the following files and directories are created:
      - input_prefixXXXX.vthb. This file provides the metadata to describe
        how AMR patches are represented in the .vti files (including a relative
        path to the .vti files).
      - directory: input_prefixXXXX containing multiple files called
        input_prefixXXXX_<level>_<patch>.vti. <level> represents the AMR level
        and <patch> indicates the AMR patch number at the given AMR level.


https://vtk.org/doc/nightly/html/classvtkOverlappingAMR.html
https://lorensen.github.io/VTKExamples/site/Python/CompositeData/OverlappingAMR/

    To open in paraview, choose the group of .vthb files, not the group of
    folders. This will be read in as cell data. In order to use filters like
    WarpByScalar you must use the CellDataToPointData filter first.

    :Input:
     - *solution* - (:class:`~pyclaw.solution.Solution`) Pyclaw object to be
       output
     - *frame* - (int) Frame number
     - *path* - (string) Root path
     - *file_prefix* - (string) Prefix for the file name. ``default = 'claw'``
     - *write_aux* - (bool) Not implemented.
     - *options* - (dict) Not implemented.
     - *write_p* - (bool) Not implemented.

    Note that some keyword arguments are not used. This is to maintain
    compatibility with the function signature expected by
    :py:class:`~pyclaw.Solution`

    Notes on what is not yet implemented
        - Add options for writing aux files.
        - Consider making an equilvalent vtk.read function.
    """
    assert(isinstance(frame, int))
    assert(isinstance(solution, Solution))
    _set_overlapped_status(solution)

    global_origin = solution.state.patch.lower_global  # base patch
    global_origin.append(0.0)  # append z
    global_origin = np.array(global_origin)
    levels = [state.patch.level-1 for state in solution.states]

    # shift base level to 0, since the base level in clawpack
    # is 1 while the base level in VTK is 0
    level_count = {}
    level_spacing = {}  # spacing of each level
    for i, level in enumerate(levels):
        if level in level_count.keys():
            level_count[level] = level_count[level] + 1
        else:
            level_count[level] = 1
            spacing = solution.states[i].patch.delta
            spacing.append(spacing[0])  # dz = dx
            spacing = np.array(spacing)
            level_spacing[level] = spacing
    numLevels = len(level_count.keys())

    # a list of num of patches at each level
    blocksPerLevel = [item[1] for item in
                     sorted(level_count.items(),
                            key=lambda a: a[0])]

    amr.Initialize(numLevels, blocksPerLevel)
    states_sorted = sorted(solution.states, key=lambda a: a.patch.level)

FIX LEVEL AND BLOCK INDEXING TO SIMPLIFY



    for level in level_count.keys():
        nbox = level_count[level]

        for index in range(box_per_level[level]):
            # ----each vtkAMRBlock can have multiple vtkAMRBox
            local_index = global_index + index
            origin = states_sorted[local_index].patch.lower_global
            origin.append(0.0)  # append z
            origin = np.array(origin)
            ndim = states_sorted[local_index].patch.num_cells_global
            ndim.append(0.0)  # mz
            ndim = np.array(ndim, dtype=np.int)
            ndim = ndim + 1  # ndim should be num of nodes
            amrbox = vtkAMRBox(origin, ndim)

            q = states_sorted[local_index].q
            for i in range(q.shape[0]-1):
                q_i = q[i, ...]
                q_i = q_i.transpose()
                amrbox.set_cell_data(q_i, "q_"+str(i))

            # this is where writing out aux files wouild happen.

            q_ol = q[-1, ...]  # last piece is used to mark overlapped cells
            q_ol = q_ol.transpose()
            amrbox.set_cell_data(q_ol, "vtkGhostType", "UInt8")

            # set vtkGhostType data
            # ghost_q = np.zeros(q[0, ...].shape, dtype=int)
            # ghost_q = ghost_q.transpose()
            # amrbox.set_cell_data(ghost_q, "vtkGhostType", data_type="UInt8")

            # shape = list(q1.shape)
            # shape.append(1)
            # point_data = np.ones( np.array(shape) + 1)
            # amrbox.set_point_data(point_data)

            # create a uniform grid.
            ug = vtk.vtkUniformGrid()
            # Geometry
            ug.SetOrigin(origin)
            ug.SetSpacing(spacing)
            ug.SetDimensions(dims)

            # Data
            scalars = vtk.vtkFloatArray()
            ug.GetPointData().SetScalars(scalars)
            MakeScalars(dims, origin, spacing, scalars)

            # make AMR box
            box1 = vtk.vtkAMRBox()

            # add AMR box and uniform grid to the overlapping AMR object.
            amr.SetAMRBox(level, block, box)
            amr.SetDataSet(level, block, ug)

        # write out.
        amr.write()
        # https://lorensen.github.io/VTKExamples/site/Python/IO/WriteXMLLinearCells/

def _set_overlapped_status(sol):
    """
    return a list, overlapped_states,
    whose entries denote overlapped status for each patch.

    @type sol:  pyclaw.Solution
    @param sol: Solution obejct of pyclaw that contains all information
                of this time step.
    @rtype:     list
    @return:    add a component to the solution q,
                which contains overlapped status of each patch

    """
    levels = [state.patch.level-1 for state in sol.states]
    # shift base level to 0
    level_count = {}
    level_spacing = {}  # spacing of each level
    for i, level in enumerate(levels):
        if level in level_count.keys():
            level_count[level] = level_count[level] + 1
        else:
            level_count[level] = 1
            spacing = sol.states[i].patch.delta
            spacing.append(spacing[0])  # dz = dx
            spacing = np.array(spacing)
            level_spacing[level] = spacing

    # a list of num of patches at each level
    box_per_level = [item[1] for item in
                     sorted(level_count.items(),
                            key=lambda a: a[0])]
    box_per_level = np.array(box_per_level)

    """
    @type level_count:       dictionary
    @variable level_count:   a dictionary that maps levels
                             to number of patches of certain levels
                             e.g. {0:1, 1:2, 2:12}
    @type num_levels:        int
    @variable num_levels:    number of levels in total
    @type box_per_level:     list
    @variable box_per_level: [number of patches on level 0,
                              number of patches on level1, ...]
    """
    for state in sol.states:
        level = state.patch.level-1
        xlower_coarse = state.patch.dimensions[0].lower
        # xupper_coarse = state.patch.dimensions[0].upper
        ylower_coarse = state.patch.dimensions[1].lower
        # yupper_coarse = state.patch.dimensions[1].upper
        dx = state.patch.delta[0]
        dy = state.patch.delta[1]
        nx = state.patch.num_cells_global[0]
        ny = state.patch.num_cells_global[1]
        # in overlapped_status, entry with value 0 denotes
        # that the cell is not overlapped
        # entry with value 8 denotes that the cell is overlapped
        overlapped_status = np.zeros((1, state.q.shape[1], state.q.shape[2]))
        # convert from Fortran-Style to C-style
        # overlapped_status = overlapped_status.transpose(0, 2, 1)
        # In the future, efficiency of this part can be improved
        # by mapping grid levels to
        # a list of states of corresponding levels.
        # Otherwise, we need to scan each states in each outer loop as below
        for state_fine in sol.states:
            # find states with grid level of one higher
            if ((state_fine.patch.level-1) == level + 1):
                xlower_fine = state_fine.patch.dimensions[0].lower
                xupper_fine = state_fine.patch.dimensions[0].upper
                ylower_fine = state_fine.patch.dimensions[1].lower
                yupper_fine = state_fine.patch.dimensions[1].upper
                x_idx_lower = \
                    max(int(round((xlower_fine - xlower_coarse) /
                            float(dx))), 0)
                x_idx_upper = \
                    min(int(round((xupper_fine-xlower_fine)/dx)) +
                        x_idx_lower, nx)
                y_idx_lower = \
                    max(int(round((ylower_fine - ylower_coarse) /
                            float(dy))), 0)
                y_idx_upper = \
                    min(int(round((yupper_fine-ylower_fine)/dy)) +
                        y_idx_lower, ny)
                # set these cells to 8
                overlapped_status[0, x_idx_lower:x_idx_upper,
                                  y_idx_lower:y_idx_upper].fill(8)

            else:
                continue
        # state.q is in fortran style
        # state.q = np.vstack((state.q, overlapped_status.transpose(0, 2, 1)))
        state.q = np.vstack((state.q, overlapped_status))
