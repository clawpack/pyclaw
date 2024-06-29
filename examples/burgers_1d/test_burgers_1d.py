from . import burgers_1d
import numpy as np
import os

thisdir = os.path.dirname(__file__)
path_expected_sols = os.path.join(thisdir,'expected_sols.npy')
expected_sols_dict = np.load(path_expected_sols,allow_pickle=True).item()

def error(test_name,**kwargs):
    """
    Compute difference between initial and final solutions.
    This should vanish due to the periodic boundary conditions
    """

    claw = burgers_1d.setup(outdir=None,**kwargs)
    claw.run()

    # tests are done across the entire domain of q normally
    qtest = claw.frames[claw.num_output_times].state.get_q_global()
    qtest = qtest.reshape([-1])
    dx = claw.solution.domain.grid.delta[0]

    assert test_name in expected_sols_dict.keys(), f"Test name {test_name} not found in {path_expected_sols}"
    qexpected = expected_sols_dict[test_name]
    
    diff = dx*np.sum(np.abs(qtest-qexpected))
    return diff


class TestBurgers1D:
    def test_python_classic(self):
        assert error(test_name="python_classic",kernel_language='Python',solver_type='classic')<1e-6

    def test_fortran_classic(self):
        assert error(test_name="fortran_classic",kernel_language='Fortran',solver_type='classic')<1e-6

    def test_python_sharpclaw(self):
        assert error(test_name="python_sharpclaw",kernel_language='Python',solver_type='sharpclaw')<1e-6

    def test_fortran_sharpclaw(self):
        assert error(test_name="fortran_sharpclaw",kernel_language='Fortran',solver_type='sharpclaw')<1e-6



