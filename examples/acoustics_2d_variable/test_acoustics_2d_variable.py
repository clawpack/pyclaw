from . import acoustics_2d_interface
import numpy as np
from clawpack.pyclaw.util import check_diff
import os

def check_error(**kwargs):

    solver_type = kwargs['solver_type']
    claw = acoustics_2d_interface.setup(disable_output=True,**kwargs)
    claw.run()
    test_pressure = claw.frames[-1].q[0,:,:]

    thisdir = os.path.dirname(__file__)
    expected_pressure = np.loadtxt(os.path.join(thisdir, 'pressure_%s.txt' % solver_type))
    test_err = np.max(np.abs(expected_pressure[:].reshape(-1) - 
                             test_pressure[:].reshape(-1)))
    print(test_err)
    return check_diff(0, test_err, abstol=1e-1)


class TestAcoustics2D:
    def test_classic(self):
        assert check_error(solver_type='classic',num_cells=(50,50))==None

    def test_sharpclaw(self):
        assert check_error(solver_type='sharpclaw',num_cells=(50,50))==None

    def test_sharpclaw_multistep(self):
        assert check_error(solver_type='sharpclaw',num_cells=(50,50),
                           time_integrator='SSPLMMk2')==None
