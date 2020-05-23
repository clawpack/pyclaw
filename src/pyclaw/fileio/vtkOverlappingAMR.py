import numpy as np
"""
This module defines several classes of vtk structures
"""


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
    def write_ascii(self, filename):
        op_file = open(filename + '.vthb', 'w')
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
        op_file.close()
        for i in range(self.num_levels):
            # pdb.set_trace()
            self.blocks[i].write_root_file(filename)
        op_file = open(filename + '.vthb', 'a')
        op_file.write('  </vtkOverlappingAMR>\n')
        op_file.write('</VTKFile>')
        op_file.close()


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

    def write_root_file(self, filename):
        op_file = open(filename + '.vthb', 'a')
        op_file.write('    <Block level=\"' + str(self.level) +
                      '\" spacing=\"' +
                      str(self.spacing[0]) + ' ' +
                      str(self.spacing[1]) + ' ' +
                      str(self.spacing[2]) + '\">\n')
        import os
        if filename not in os.listdir('./'):
            os.mkdir(filename)
        for i in range(self.nbox):
            # pdb.set_trace()
            boundary_index = self.boxes[i].get_global_boundary_index(
                             self.global_origin)
            child_path = filename + '/' + filename + \
                '_' + str(self.level) + '_' + str(i) + '.vti'
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
            self.boxes[i].write_child_ascii(child_path)
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

    def write_child_ascii(self, filename):
        op_file = open(filename, 'w')
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
        op_file.close()
