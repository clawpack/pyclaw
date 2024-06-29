"Runs the shocktube test problem with the Python and Fortran HLL solvers."

def test_shocktube():
    "Shock tube test (Euler 1D)"
    import numpy as np
    import os
    from . import shocktube

    thisdir = os.path.dirname(__file__)

    claw = shocktube.setup(kernel_language='Python',outdir=os.path.join(thisdir,'./_output'))
    claw.run()
    hllc = claw.solution.state.get_q_global()

    claw = shocktube.setup(kernel_language='Fortran',outdir=os.path.join(thisdir,'./_output'))
    claw.run()
    hlle = claw.solution.state.get_q_global()

    # Differences of about 0.03 appear near the contact wave
    assert np.max(np.abs(hlle-hllc)) < 0.04

    if hllc is not None:
        expected_density = np.loadtxt(os.path.join(thisdir,'shocktube_regression_density.txt'))
        test_density = hllc[0,:]
        test_err = np.linalg.norm(expected_density-test_density)
        assert test_err < 1.e-4

