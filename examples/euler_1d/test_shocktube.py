from __future__ import absolute_import
import numpy as np
import os
from . import shocktube

"Runs the shocktube test problem with the Python HLLC solver."

claw = shocktube.setup(kernel_language='Python')
claw.run()
test_solution = claw.solution.state.get_q_global()

if test_solution is not None:
    thisdir = os.path.dirname(__file__)
    expected_density = np.loadtxt(os.path.join(thisdir,'shocktube_regression_density.txt'))
    test_density = test_solution[0,:]
    test_err = np.linalg.norm(expected_density-test_density)
    assert test_err < 1.e-4
