from . import advection_annulus
import numpy as np
import os

thisdir = os.path.dirname(__file__)
path_expected_sols = os.path.join(thisdir,'expected_sols.npy')
expected_sols_dict = np.load(path_expected_sols,allow_pickle=True).item()

def error(test_name,**kwargs):
    """
    Compute difference between test and expected solutions.
    Using L1 norm over the physical domain.
    """
    claw = advection_annulus.setup(outdir=None,**kwargs)
    claw.run()
    qtest = claw.frames[claw.num_output_times].state.get_q_global().reshape([-1])
    dx = claw.solution.domain.grid.delta[0]
    dy = claw.solution.domain.grid.delta[1]
    relative_area = claw.frames[0].state.aux[2].reshape([-1])
    area = relative_area*dx*dy

    assert test_name in expected_sols_dict.keys(), f"Test name {test_name} not found in {path_expected_sols}"
    qexpected = expected_sols_dict[test_name]
    diff_L1 = np.sum(abs(area*(qtest-qexpected)))
    return diff_L1

class TestAdvectionAnnulus2D:

    def test_classic(self):
        assert error(test_name='classic', solver_type='classic')<1e-6
