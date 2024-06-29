from . import quadrants, shock_bubble_interaction, shock_forward_step
import numpy as np
from pytest import importorskip
import os

thisdir = os.path.dirname(__file__)
path_expected_sols = os.path.join(thisdir,'expected_sols.npy')
expected_sols_dict = np.load(path_expected_sols,allow_pickle=True).item()

def error(test_name,**kwargs):
    if "quadrants" in test_name:
        claw = quadrants.setup(**kwargs)
    elif "shock_bubble_interaction" in test_name:
        claw = shock_bubble_interaction.setup(**kwargs)
    elif "shock_forward_step" in test_name:
        claw = shock_forward_step.setup(**kwargs)
    else:
        raise ValueError(f"Test name {test_name} not recognized")
    claw.run()
    density_test = claw.frames[claw.num_output_times].state.get_q_global()[0].reshape([-1])
    dx = claw.solution.domain.grid.delta[0]
    dy = claw.solution.domain.grid.delta[1]
    assert test_name in expected_sols_dict.keys(), f"Test name {test_name} not found in {path_expected_sols}"
    density_expected = expected_sols_dict[test_name]
    diff_L1 = dx*dy*np.sum(np.abs(density_test-density_expected))
    return diff_L1

class TestEuler2D:
    def test_quadrants_hlle(self):
        assert error(test_name='quadrants_hlle', riemann_solver='hlle')<1e-6

    def test_quadrants_roe(self):
        assert error(test_name='quadrants_roe', riemann_solver='roe')<1e-6

    def test_shock_bubble_interaction(self):
        _ = importorskip("scipy") # skip test if scipy is not installed
        assert error(test_name='shock_bubble_interaction', mx=160, my=40, tfinal=0.2, num_output_times=1, disable_output=True)<1e-6
    
    # def test_shock_forward_step(self):
    #     assert error(test_name='shock_forward_step',tfinal=0.1)<1e-6