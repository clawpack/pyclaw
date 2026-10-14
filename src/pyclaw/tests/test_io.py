import os
import glob
import numpy as np
from clawpack.pyclaw import Solution
from clawpack.pyclaw.util import check_solutions_are_same

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
        for fmt, sol in s.items():
            check_solutions_are_same(sol,ref_sol)

    def test_io_to_vtk(self):
        # since the VTK only has a write and no read, I'm making sure that
        # the output matches sample output.

        # Input format (Binary)
        regression_dir = os.path.join(self.test_data_dir,'./advection_2d_binary')

        # read in from ascii.
        ref_sol = self.solution
        ref_sol.read(0,path=regression_dir,file_format="binary")

        # write out to vtk in io_test_dir
        io_test_dir = os.path.join(self.this_dir,'./io_test')
        ref_sol.write(0, file_format="vtk", path=io_test_dir)

        # Known correct data:
        comparison_dir = os.path.join(self.test_data_dir,'./advection_2d_vtk')

        # use glob to find all vthb and vti files.
        # compare  line by line.
        assert self.compare_vtk(comparison_dir, io_test_dir)==True

    def compare_vtk(self, dir1, dir2):
        # compare vtk files in two directories, return True if files are
        # equivalent.

        # assumption is that there are vthb files at top level and one level
        # down in folders are vthb files.

        # number of vthb and vti files is known. assert them.
        dir1_vthb_files = glob.glob(os.path.join(dir1, "*.vthb"))
        assert len(dir1_vthb_files) == 1

        dir1_vti_files = glob.glob(os.path.join(dir2, *["*", "*.vti"]))
        assert len(dir1_vti_files) == 17

        # for all files do line by line comparison.
        for dir1_file in dir1_vti_files + dir1_vthb_files:
            dir2_file = dir1_file.replace(dir1, dir2)

            # open files
            with open(dir1_file, 'r') as  f1:
                lines1 = f1.readlines()
            with open(dir2_file, 'r') as  f2:
                lines2 = f2.readlines()

            # compare each line.
            for l1, l2 in zip(lines1, lines2):
                assert l1.strip() == l2.strip()

            # delete lines for this file.
            del lines1, lines2

        # Only if all assertions succeed, return True.
        return True
