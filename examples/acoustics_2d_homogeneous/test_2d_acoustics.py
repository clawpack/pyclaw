from . import acoustics_2d
import numpy as np
from clawpack.pyclaw.util import check_diff
import os

def check_error(data_filename,**kwargs):

    claw = acoustics_2d.setup(disable_output=True,**kwargs)
    claw.run()

    test_pressure = claw.frames[-1].q[0,:,:]
    thisdir = os.path.dirname(__file__)
    expected_pressure = np.loadtxt(os.path.join(thisdir,data_filename))
    return check_diff(expected_pressure, test_pressure, reltol=1e-3,
                        delta=claw.solution.grid.delta)


class TestAcoustics2D:
    def test_classic(self):
        assert check_error('verify_classic.txt',solver_type='classic')==None

    def test_classic_ptwise(self):
        assert check_error('verify_classic.txt',solver_type='classic',ptwise=True)==None

    def test_sharpclaw(self):
        assert check_error('verify_sharpclaw.txt',solver_type='sharpclaw')==None

    def test_sharpclaw_multistep(self):
        assert check_error('verify_sharpclaw_lmm.txt',solver_type='sharpclaw',
                           time_integrator='SSPLMMk2')==None
