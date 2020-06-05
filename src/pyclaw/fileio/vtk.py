#!/usr/bin/python

# Author: Xinsheng (Shawn) Qin
# date: 03/2016
# Modified by Katy Barnhart
# date 05/2020

import sys
import os
import numpy as np
from clawpack.pyclaw import Solution


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
    num_levels = len(level_count.keys())

    # a list of num of patches at each level
    box_per_level = [item[1] for item in
                     sorted(level_count.items(),
                            key=lambda a: a[0])]
    box_per_level = np.array(box_per_level)
    AMRdata = vtkOverlappingAMR(global_origin, num_levels, box_per_level)

    states_sorted = sorted(solution.states, key=lambda a: a.patch.level)
    global_index = 0
#################################################
    for level in level_count.keys():
        nbox = level_count[level]
        block = vtkAMRBlock(level, nbox, level_spacing[level], global_origin)

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
            block.attached_amrbox(amrbox)
            AMRdata.attached_block(level, block)

        global_index = global_index + box_per_level[level]

    filename = file_prefix+str(frame).zfill(4)
    AMRdata.write_ascii(path, filename)


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
    # num_levels = len(level_count.keys())

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


class vtkOverlappingAMR(object):
    r""" A vtkOverlappingAMR object """
    def __init__(self, origin, num_levels, box_per_level):
        r"""
        Input:
        origin - numpy array (x, y, z)
        num_levels - int
        box_per_level - numpy array (int, int, ...)
        """
        self.xml_version = "1.0"
        self.vtk_file_type_version = "1.1"
        # self.vtk_file_type_version = "0.1"
        self.byte_order = "LittleEndian"
        self.header_type = "UInt32"
        # self.header_type = "UInt64"
        self.compressor = "vtkZLibDataCompressor"
        self.grid_description = "XY"  # todo: this is just for 2d grid
        assert(isinstance(origin, np.ndarray)), \
            "origin should be a numpy array"
        assert(isinstance(box_per_level, np.ndarray)), \
            "box_per_level should be a numpy array"
        assert(isinstance(num_levels, int)), \
            "num_levels should be an int"

        self.origin = origin
        self.num_levels = num_levels
        self.box_per_level = box_per_level

        self.blocks = {}
        self.spacing = np.zeros((num_levels, 3))

    def set_spacing(self, level, h):
        assert(isinstance(h, np.ndarray)), \
            "h should be a numpy array"
        assert(isinstance(level, int)), \
            "level should be an int"
        self.spacing[level, :] = h

    def attached_block(self, level, vtk_amr_block):
        assert(isinstance(level, int)), \
            "level should be an int"
        assert(isinstance(vtk_amr_block, vtkAMRBlock)), \
            "vtk_amr_block should be an vtkAMRBlock object"
        self.blocks[level] = vtk_amr_block

    def set_xml_version(self, str1):
        assert (type(str1) == str), "str1 should be a string type!"
        self.xml_version = str1

    # todo: define other setting functions
    def write_ascii(self, path, filename):
        """
        Write the .vthb function

        path
        filename
        """
        import os

        with open(os.path.join(path, filename) + '.vthb', 'w') as op_file:
            op_file.write('<?xml version=\"'+self.xml_version+'\"?>\n')
            op_file.write('<VTKFile type=\"vtkOverlappingAMR" version=\"' +
                          self.vtk_file_type_version+'\" byte_order=\"' +
                          self.byte_order + '\" header_type=\"' +
                          self.header_type + '\" compressor=\"' +
                          self.compressor + '\">\n')
            op_file.write('  <vtkOverlappingAMR origin=\"' +
                          str(self.origin[0]) + ' ' +
                          str(self.origin[1]) + ' ' +
                          str(self.origin[2]) +
                          '\" grid_description=\"' +
                          self.grid_description + '\">\n')

        # write individual block files in folder. This appends information to
        # the .vthb file.
        for i in range(self.num_levels):
            # pdb.set_trace()
            self.blocks[i].write_root_file(path, filename)

        # write end of .vthb file.
        with open(os.path.join(path, filename) + '.vthb', 'a') as op_file:
            op_file.write('  </vtkOverlappingAMR>\n')
            op_file.write('</VTKFile>')


class vtkAMRBlock(object):
    r"""
    A vtkAMRBlock contains all amr_box at the same level.
    """
    def __init__(self, level, nbox, spacing, global_origin):
        if not isinstance(level, int):
            raise TypeError
        if not isinstance(nbox, int):
            raise TypeError
        if not isinstance(spacing, np.ndarray):
            raise TypeError
        if not isinstance(global_origin, np.ndarray):
            raise TypeError
        self.level = level
        self.nbox = nbox
        self.boxes = []
        self.spacing = spacing
        self.global_origin = global_origin

    def attached_amrbox(self, amrbox):
        assert(isinstance(amrbox, vtkAMRBox)), \
            "vtk_amr_block should be an vtkAMRBox object"
        amrbox.set_spacing(self.spacing)
        self.boxes.append(amrbox)

    def write_root_file(self, path, filename):
        # write block info to .vthb file.
        import os

        op_file = open(os.path.join(path, filename) + '.vthb', 'a')
        op_file.write('    <Block level=\"' + str(self.level) +
                      '\" spacing=\"' +
                      str(self.spacing[0]) + ' ' +
                      str(self.spacing[1]) + ' ' +
                      str(self.spacing[2]) + '\">\n')

        if filename not in os.listdir(path):
            os.makedirs(os.path.join(path, filename))

        for i in range(self.nbox):
            # pdb.set_trace()
            boundary_index = self.boxes[i].get_global_boundary_index(
                             self.global_origin)

            fn = filename+'_'+str(self.level)+'_'+str(i)+'.vti'

            child_path = os.path.join(filename, fn)
            op_file.write('      <DataSet index=\"' +
                          str(i) + '\" amr_box=\"' +
                          str(boundary_index[0]) + ' ' +
                          str(boundary_index[1]) + ' ' +
                          str(boundary_index[2]) + ' ' +
                          str(boundary_index[3]) + ' ' +
                          str(boundary_index[4]) + ' ' +
                          str(boundary_index[5]) +
                          '\" file=\"' + child_path + '\">\n')
            op_file.write('      </DataSet>\n')
            # write data in children directory
            self.boxes[i].write_child_ascii(path, child_path)
        op_file.write('    </Block>\n')
        op_file.close()


class vtkAMRBox(object):
    r"""
    An vtk amr_box object.
    """
    def __init__(self, origin, ndim):
        assert(isinstance(origin, np.ndarray)), \
            "origin should be a numpy array"
        assert(isinstance(ndim, np.ndarray)), \
            "ndim should be a numpy array"

        self.origin = origin
        self.ndim = ndim  # (nx, ny, nz) - number of nodes
        self.spacing = None  # (dx, dy, dz)

        self.cell_data = []
        self.cell_data_name = []
        self.cell_data_type = []

        self.point_data = []
        self.point_data_name = []

        self.xml_version = "1.0"
        # self.vtk_file_type_version = "0.1"
        self.vtk_file_type_version = "1.0"
        # self.vtk_file_type_version = "1.1"
        # self.vtk_file_type_version = "2.0"
        self.byte_order = "LittleEndian"
        self.header_type = "UInt32"
        # self.header_type = "UInt64"
        self.compressor = "vtkZLibDataCompressor"

    def set_spacing(self, spacing):
        assert(isinstance(spacing, np.ndarray)), \
            "spacing should be a numpy array"
        self.spacing = spacing

    def get_ndim(self):
        return self.ndim

    def set_point_data(self, data, name):
        if not isinstance(data, np.ndarray):
            raise TypeError
        assert(isinstance(name, str)), "name must be a str."
        self.point_data.append(data)
        self.point_data_name.append(name)

    def set_cell_data(self, data, name, data_type="Float64"):
        if not isinstance(data, np.ndarray):
            raise TypeError
        assert(isinstance(name, str)), "name must be a str."
        assert(isinstance(data_type, str)), "data_type must be a str."
        self.cell_data.append(data)
        self.cell_data_name.append(name)
        self.cell_data_type.append(data_type)

    def get_global_boundary_index(self, global_origin):
        relative_pos = self.origin - global_origin
        # round to nearest int
        i_low = int(round(relative_pos[0]/self.spacing[0]))
        j_low = int(round(relative_pos[1]/self.spacing[1]))
        k_low = int(round(relative_pos[2]/self.spacing[2]))

        i_high = i_low + self.ndim[0] - 2  # todo: figure out why it is -2
        j_high = j_low + self.ndim[1] - 2  # ndim is num of nodes
        k_high = k_low + self.ndim[2] - 2
        # pdb.set_trace()

        return np.array([i_low, i_high, j_low, j_high, k_low, k_high])

    def get_local_boundary_index(self):
        i_low = 0
        j_low = 0
        k_low = 0

        i_high = i_low + self.ndim[0] - 1  # todo: figure out why it is -1
        j_high = j_low + self.ndim[1] - 1
        k_high = k_low + self.ndim[2] - 1

        return np.array([i_low, i_high, j_low, j_high, k_low, k_high])

    def write_child_ascii(self, path, filename):
        import os
        with open(os.path.join(path, filename), 'w') as op_file:
            op_file.write('<?xml version=\"'+self.xml_version+'\"?>\n')
            op_file.write('<VTKFile type=\"ImageData\" version=\"' +
                          self.vtk_file_type_version+'\" byte_order=\"' +
                          self.byte_order + '\" header_type=\"' +
                          self.header_type + '\" compressor=\"' +
                          self.compressor + '\">\n')
            extent = self.get_local_boundary_index()
            op_file.write('  <ImageData WholeExtent=\"' +
                          str(extent[0]) + ' ' +
                          str(extent[1]) + ' ' +
                          str(extent[2]) + ' ' +
                          str(extent[3]) + ' ' +
                          str(extent[4]) + ' ' +
                          str(extent[5]) + '\" ' +
                          'Origin=\"' +
                          str(self.origin[0]) + ' ' +
                          str(self.origin[1]) + ' ' +
                          str(self.origin[2]) +
                          '\" Spacing=\"' +
                          str(self.spacing[0]) + ' ' +
                          str(self.spacing[1]) + ' ' +
                          str(self.spacing[2]) + '\">\n')
            op_file.write('  <Piece Extent=\"' +
                          str(extent[0]) + ' ' +
                          str(extent[1]) + ' ' +
                          str(extent[2]) + ' ' +
                          str(extent[3]) + ' ' +
                          str(extent[4]) + ' ' +
                          str(extent[5]) + '\">\n')
            op_file.write('    <PointData>\n')
            if len(self.point_data) != 0:  # write point data
                for data_name, point_data in zip(self.point_data_name,
                                                 self.point_data):
                    op_file.write('      <DataArray type=\"Float64\" Name=\"' +
                                  data_name + '\" ' +
                                  'format=\"ascii\" RangeMin=\"' +
                                  str(point_data.min()) + '\" ' +
                                  'RangeMax=\"' + str(point_data.max()) + '\">\n')
                    data_flat = point_data.flatten()
                    op_file.write('       ')
                    for i, item in enumerate(data_flat):
                        # write every 6 numbers in a line
                        if (i % 6 == 0) and (i != 0):
                            op_file.write('\n')
                            op_file.write('       ')
                        op_file.write(str(item) + ' ')
                    op_file.write('\n')
                    op_file.write('      </DataArray>\n')
            op_file.write('    </PointData>\n')
            op_file.write('    <CellData>\n')
            if len(self.cell_data) != 0:  # write cell data
                for data_name, cell_data, data_type in \
                        zip(
                            self.cell_data_name,
                            self.cell_data,
                            self.cell_data_type
                           ):
                    op_file.write('      <DataArray type=\"' +
                                  data_type + '\" Name=\"' +
                                  data_name + '\" ' +
                                  'format=\"ascii\" RangeMin=\"' +
                                  str(cell_data.min()) + '\" ' +
                                  'RangeMax=\"' + str(cell_data.max()) + '\">\n')
                    data_flat = cell_data.flatten()
                    op_file.write('       ')
                    for i, item in enumerate(data_flat):
                        # write every 6 numbers in a line
                        if (i % 6 == 0) and (i != 0):
                            op_file.write('\n')
                            op_file.write('       ')
                        # we should write data as integer
                        if ("Int" in data_type):
                            op_file.write(str(int(item)) + ' ')
                        else:
                            op_file.write(str(item) + ' ')
                    op_file.write('\n')
                    op_file.write('      </DataArray>\n')
            op_file.write('    </CellData>\n')
            op_file.write('  </Piece>\n')
            op_file.write('  </ImageData>\n')
            op_file.write('</VTKFile>')
