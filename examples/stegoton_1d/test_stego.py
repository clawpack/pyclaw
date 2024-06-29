from . import stegoton
import numpy as np
from pytest import importorskip
import os

thisdir = os.path.dirname(__file__)
path_expected_sols = os.path.join(thisdir,'expected_sols.npy')
expected_sols_dict = np.load(path_expected_sols,allow_pickle=True).item()

def error(test_name,**kwargs):
    claw = stegoton.setup(tfinal=50.0, outdir=None, **kwargs)
    claw.run()
    strain_test = claw.frames[claw.num_output_times].state.get_q_global()[0].reshape([-1])
    dx = claw.solution.domain.grid.delta[0]
    assert test_name in expected_sols_dict.keys(), f"Test name {test_name} not found in {path_expected_sols}"
    strain_expected = expected_sols_dict[test_name]
    diff_L1 = dx*np.sum(np.abs(strain_test-strain_expected))
    return diff_L1

class TestStegoton1D:
    def test_fortran_classic(self):
        assert error(test_name='fortran_classic', kernel_language='Fortran', solver_type="classic")<1e-6

    def test_fortran_sharpclaw(self):
        assert error(test_name='fortran_sharpclaw', kernel_language='Fortran', solver_type="sharpclaw")<1e-6
    
    def test_python_classic(self):
        assert error(test_name='python_classic', kernel_language='Python', solver_type="classic")<1e-6
    
    def test_python_sharpclaw(self):
        assert error(test_name='python_sharpclaw', kernel_language='Python', solver_type="sharpclaw")<1e-6