import os
from clawpack import pyclaw
from clawpack.pyclaw import examples
import numpy as np

thisdir = os.path.dirname(__file__)
file_formats = ['hdf5','ascii']

def test_io_from_binary():
    # Read regression data
    regression_dir = os.path.join(thisdir,'./test_data/advection_2d_binary')
    read_write_and_compare(file_formats,regression_dir,'binary',0)

def test_io_from_hdf5():
    regression_dir = os.path.join(thisdir,'./test_data/Sedov_regression_hdf')
    read_write_and_compare(file_formats,regression_dir,'hdf5',1)

def test_io_from_hdf5_with_aux():
    regression_dir = os.path.join(thisdir,'./test_data/advection_2d_with_aux')
    read_write_and_compare(file_formats,regression_dir,'hdf5',0,aux=True)



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


def check_solutions_are_same(sol_a,sol_b):
    assert len(sol_a.states) == len(sol_b.states)
    assert sol_a.t == sol_b.t
    for state in sol_a.states:
        for ref_state in sol_b.states:
            if ref_state.patch.patch_index == state.patch.patch_index:
                break

        # Required state attributes
        assert np.linalg.norm(state.q - ref_state.q) < 1.e-6 # Not sure why this can be so large
        if ref_state.aux is not None:
            assert np.linalg.norm(state.aux - ref_state.aux) < 1.e-16
        for attr in ['t', 'num_eqn', 'num_aux']:
            assert getattr(state,attr) == getattr(ref_state,attr)
        # Optional state attributes
        for attr in ['patch_index', 'level']:
            if hasattr(ref_state,attr):
                assert getattr(state,attr) == getattr(ref_state,attr)

        patch = state.patch
        ref_patch = ref_state.patch
        # Required patch attributes
        for attr in ['patch_index', 'level']:
            assert getattr(patch,attr) == getattr(ref_patch,attr)

        dims = patch.dimensions
        ref_dims = ref_patch.dimensions
        for dim, ref_dim in zip(dims,ref_dims):
            # Required dim attributes
            for attr in ['num_cells','lower','delta']:
                assert getattr(dim,attr) == getattr(ref_dim,attr)
            # Optional dim attributes
            for attr in ['units','on_lower_boundary','on_upper_boundary']:
                if hasattr(ref_dim,attr):
                    assert getattr(dim,attr) == getattr(ref_dim,attr)
