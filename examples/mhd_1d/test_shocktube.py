"Runs the shocktube test problem with the Fortran Roe solver."

def test_shocktube():
    "Shock tube test (1D MHD with Roe solver)"
    import numpy as np
    import os
    from . import shocktube
    from clawpack.riemann.mhd_1D_constants import B_2, B_3

    claw = shocktube.setup(kernel_language='Fortran')
    claw.run()
    roe = claw.solution.state.get_q_global()

    if roe is not None:
        thisdir = os.path.dirname(__file__)
        expected_B2 = np.loadtxt(os.path.join(thisdir, 'shocktube_regression_B2.txt'))
        expected_B3 = np.loadtxt(os.path.join(thisdir, 'shocktube_regression_B3.txt'))
        test_B2 = roe[B_2, :]
        test_B3 = roe[B_3, :]
        test_err = np.linalg.norm(expected_B2 - test_B2) + np.linalg.norm(expected_B3 - test_B3)
        assert test_err < 1.e-4

