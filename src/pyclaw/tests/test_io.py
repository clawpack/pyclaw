from __future__ import absolute_import
import os
import numpy as np
from clawpack.pyclaw import Solution
from clawpack.pyclaw.util import check_solutions_are_same
import six

class IOTest():
    @property
    def solution(self):
        return Solution()

    @property
    def file_formats(self):
        return ['hdf5','ascii']

    @property
    def this_dir(self):
        return os.path.dirname(os.path.abspath(__file__))

    @property
    def test_data_dir(self):
        return os.path.join(self.this_dir, './test_data')

    def test_io_from_binary(self):
        # Read regression data
        regression_dir = os.path.join(self.test_data_dir,'./advection_2d_binary')
        self.read_write_and_compare(self.file_formats,regression_dir,'binary',0)

    def test_io_from_hdf5(self):
        regression_dir = os.path.join(self.test_data_dir,'./Sedov_regression_hdf')
        self.read_write_and_compare(self.file_formats,regression_dir,'hdf5',1)

    def test_io_from_hdf5_with_aux(self):
        regression_dir = os.path.join(self.test_data_dir,'./advection_2d_with_aux')
        self.read_write_and_compare(self.file_formats,regression_dir,'hdf5',0,aux=True)

    def read_write_and_compare(self, file_formats,regression_dir,regression_format,frame_num,aux=False):
        r"""Test IO file formats:
            - Reading in an HDF file
            - Writing  files in all formats (for now just ASCII, HDF5)
            - Reading those files back in
            - Checking that all the resulting Solution objects are identical
        """
        ref_sol = self.solution
        ref_sol.read(frame_num,path=regression_dir,file_format=regression_format,read_aux=aux)
        if aux:
            assert (ref_sol.state.aux is not None)

        # Write solution file in each format
        io_test_dir = os.path.join(self.this_dir,'./io_test')
        for fmt in file_formats:
            ref_sol.write(frame_num,path=io_test_dir,file_format=fmt,write_aux=aux)

        # Read solutions back in
        s = {}
        for fmt in file_formats:
            s[fmt] = self.solution
            s[fmt].read(frame_num,path=io_test_dir,file_format=fmt,write_aux=aux)

        # Compare solutions
        # Probably better to do this by defining __eq__ for each class
        for fmt, sol in six.iteritems(s):
            check_solutions_are_same(sol,ref_sol)