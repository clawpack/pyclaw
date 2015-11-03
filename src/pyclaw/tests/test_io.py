import os
from clawpack import pyclaw
import numpy as np
from clawpack.pyclaw.util import check_solutions_are_same
from clawpack.pyclaw.util import gen_io_variants

file_formats = ['hdf5','ascii']
thisdir = os.path.dirname(__file__)

def test_io_from_binary():
    # Read regression data
    regression_dir = os.path.join(thisdir,'./test_data/advection_2d_binary')
    for test in gen_io_variants(read_write_and_compare,regression_dir=regression_dir,file_formats=file_formats,regression_format='binary',frame_num=0):
        yield test

def test_io_from_hdf5():
    regression_dir = os.path.join(thisdir,'./test_data/Sedov_regression_hdf')
    for test in gen_io_variants(read_write_and_compare,regression_dir=regression_dir,file_formats=file_formats,regression_format='hdf5',frame_num=1):
        yield test

def test_io_from_hdf5_with_aux():
    regression_dir = os.path.join(thisdir,'./test_data/advection_2d_with_aux')
    for test in gen_io_variants(read_write_and_compare,regression_dir=regression_dir,file_formats=file_formats,regression_format='hdf5',frame_num=0,aux=True):
        yield test

def read_write_and_compare(file_formats,regression_dir,regression_format,frame_num,aux=False):
    r"""Test IO file formats:
        - Reading in an HDF file
        - Writing  files in all formats (for now just ASCII, HDF5)
        - Reading those files back in
        - Checking that all the resulting Solution objects are identical
    """
    ref_sol = pyclaw.Solution()
    ref_sol.read(frame_num,path=regression_dir,file_format=regression_format,read_aux=aux)
    if aux:
        assert (ref_sol.state.aux is not None)

    # Write solution file in each format
    io_test_dir = os.path.join(thisdir,'./io_test')
    for fmt in file_formats:
        ref_sol.write(frame_num,path=io_test_dir,file_format=fmt,write_aux=aux)

    # Read solutions back in
    s = {}
    for fmt in file_formats:
        s[fmt] = pyclaw.Solution()
        s[fmt].read(frame_num,path=io_test_dir,file_format=fmt,write_aux=aux)

    # Compare solutions
    # Probably better to do this by defining __eq__ for each class
    for fmt, sol in s.iteritems():
        check_solutions_are_same(sol,ref_sol)