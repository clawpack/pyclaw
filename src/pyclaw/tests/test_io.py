import os
from clawpack import pyclaw
from clawpack.pyclaw import examples

thisdir = os.path.dirname(__file__)
file_formats = ['hdf5','ascii']

def test_io_from_binary():
    # Read regression data
    regression_dir = os.path.join(thisdir,'./test_data/advection_2d_binary')
    read_write_and_compare(file_formats,regression_dir,'binary',0)

def test_io_from_hdf5():
    regression_dir = os.path.join(thisdir,'./test_data/Sedov_regression_hdf')
    read_write_and_compare(file_formats,regression_dir,'hdf5',1)

def read_write_and_compare(file_formats,regression_dir,regression_format,frame_num):
    r"""Test IO file formats:
        - Reading in an HDF file
        - Writing  files in all formats (for now just ASCII, HDF5)
        - Reading those files back in
        - Checking that all the resulting Solution objects are identical
    """
    ref_sol = pyclaw.Solution()
    ref_sol.read(frame_num,path=regression_dir,file_format=regression_format)

    # Write solution file in each format
    io_test_dir = os.path.join(thisdir,'./io_test')
    for fmt in file_formats:
        ref_sol.write(frame_num,path=io_test_dir,file_format=fmt)

    # Read solutions back in
    s = {}
    for fmt in file_formats:
        s[fmt] = pyclaw.Solution()
        s[fmt].read(frame_num,path=io_test_dir,file_format=fmt)

    # Compare solutions
    # Probably better to do this by defining __eq__ for each class
    for fmt, sol in s.iteritems():
        check_solutions_are_same(sol,ref_sol)


def check_solutions_are_same(sol_a,sol_b):
    assert len(sol_a.states) == len(sol_b.states)
    assert sol_a.t == sol_b.t
    for i,state in enumerate(sol_a.states):
        for j, ref_state in enumerate(sol_b.states):
            if ref_state.patch.patch_index == state.patch.patch_index:
                break

        # Required state attributes
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
