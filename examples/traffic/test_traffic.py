from . import traffic
import numpy as np
from pytest import importorskip
import os

thisdir = os.path.dirname(__file__)
path_expected_sols = os.path.join(thisdir,'expected_sols.npy')
expected_sols_dict = np.load(path_expected_sols,allow_pickle=True).item()

def error(test_name,**kwargs):
    claw = traffic.setup( outdir=None, **kwargs)
    claw.run()
    density_test = claw.frames[claw.num_output_times].state.get_q_global()[0].reshape([-1])
    dx = claw.solution.domain.grid.delta[0]
    assert test_name in expected_sols_dict.keys(), f"Test name {test_name} not found in {path_expected_sols}"
    density_expected = expected_sols_dict[test_name]
    diff_L1 = dx*np.sum(np.abs(density_test-density_expected))
    return diff_L1

class TestTraffic1D:
    def test_classic(self):
        assert error(test_name='classic', solver_type="classic")<1e-6

    def test_sharpclaw(self):
        assert error(test_name='sharpclaw', solver_type="sharpclaw")<1e-6
    