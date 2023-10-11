"Runs the quadrants test problem with the HLLC solver."
from __future__ import absolute_import
import numpy as np
import os
from . import quadrants

def test_quadrants():
    claw = quadrants.setup(riemann_solver='hlle')
    claw.run()
    hlle = claw.solution.state.get_q_global()

    if hlle is not None:
        thisdir = os.path.dirname(__file__)
        expected_density = np.loadtxt(os.path.join(thisdir,'quadrants_regression_density.txt'))
        test_density = hlle[0,:]
        test_err = np.linalg.norm(expected_density-test_density)
        assert test_err < 1.e-4

if __name__=="__main__":
    import nose
    nose.main()
