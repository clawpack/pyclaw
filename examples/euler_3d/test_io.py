import os
from clawpack import pyclaw

def test_io():
    r"""Test IO file formats by:
        - Reading in an HDF file
        - Writing  files in all formats (for now just ASCII, HDF5)
        - Reading those files back in
        - Checking that all the resulting Solution objects are identical
    """
    file_formats = ['hdf5','ascii']
    # Read regression data
    thisdir = os.path.dirname(__file__)
    regression_dir = os.path.join(thisdir,'./Sedov_regression')

    ref_sol = pyclaw.Solution()
    ref_sol.read(1,path=regression_dir,file_format='hdf5')

    # Write solution file in each format
    io_test_dir = os.path.join(thisdir,'./io_test')
    for fmt in file_formats:
        ref_sol.write(1,path=io_test_dir,file_format=fmt)

    # Read solutions back in
    s = {}
    for fmt in file_formats:
        s[fmt] = pyclaw.Solution()
        s[fmt].read(1,path=io_test_dir,file_format=fmt)

    # Compare solutions
    # Probably better to do this by defining __eq__ for each class
    for sol in s.itervalues():
        assert sol.t == ref_sol.t
        for i,state in enumerate(sol.states):

            ref_state = ref_sol.states[i]
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

