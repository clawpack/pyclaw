"Runs the shocktube test problem with the Python and Fortran HLL solvers."
from __future__ import absolute_import

def test_shocktube():
    "Shock tube test (Euler 1D)"
    import numpy as np
    import os
    from . import shocktube

    claw = shocktube.setup(kernel_language='Python')
    claw.run()
    hllc = claw.solution.state.get_q_global()

    claw = shocktube.setup(kernel_language='Fortran')
    claw.run()
    hlle = claw.solution.state.get_q_global()

    # Differences of about 0.03 appear near the contact wave
    assert np.max(np.abs(hlle-hllc)) < 0.04

    if hllc is not None:
        thisdir = os.path.dirname(__file__)
        expected_density = np.loadtxt(os.path.join(thisdir,'shocktube_regression_density.txt'))
        test_density = hllc[0,:]
        test_err = np.linalg.norm(expected_density-test_density)
        assert test_err < 1.e-4


if __name__=="__main__":
    import nose
    nose.main()
