from . import variable_coefficient_advection
import numpy as np
from clawpack.pyclaw.util import check_diff
import os

thisdir = os.path.dirname(__file__)
path_expected_sols = os.path.join(thisdir,'expected_sols.npy')
expected_sols_dict = np.load(path_expected_sols,allow_pickle=True).item()

def error(test_name,**kwargs):
    """
    Compute difference between test and expected solutions.
    """

    claw = variable_coefficient_advection.setup(outdir=None,**kwargs)
    claw.run()

    # tests are done across the entire domain of q normally
    qtest = claw.frames[claw.num_output_times].state.get_q_global()
    qtest = qtest.reshape([-1])
    dx = claw.solution.domain.grid.delta[0]

    
    assert test_name in expected_sols_dict.keys(), f"Test name {test_name} not found in {path_expected_sols}"
    qexpected = expected_sols_dict[test_name]
    
    diff = dx*np.sum(np.abs(qtest-qexpected))
    return diff


class TestAdvectionVarCoeff1D:
    def test_python_classic(self):
        assert error(test_name="python_classic",kernel_language='Python',solver_type='classic')<1e-6

    def test_fortran_classic(self):
        assert error(test_name="fortran_classic",kernel_language='Fortran',solver_type='classic')<1e-6

    def test_python_sharpclaw(self):
        assert error(test_name="python_sharpclaw",kernel_language='Python',solver_type='sharpclaw')<1e-6

    def test_fortran_sharpclaw(self):
        assert error(test_name="fortran_sharpclaw",kernel_language='Fortran',solver_type='sharpclaw')<1e-6

    def test_sharpclaw_custom_time_integrator(self):
        #Load Butcher Tableaux for custom time integrators
        rk_methods_dict = np.load(os.path.join(thisdir,'rk_methods.npy'),allow_pickle=True).item()
        rk_names = list(rk_methods_dict.keys())
        for rk_name in rk_names:
            rk_coeffs = rk_methods_dict[rk_name]
            assert error(test_name='sharpclaw_custom_time_integrator_'+rk_name,kernel_language='Fortran',solver_type='sharpclaw',
                             time_integrator='RK',rk_coeffs=rk_coeffs)<1e-6, f"Failed for {rk_name}"


