from . import acoustics_1d
from clawpack.pyclaw.util import check_diff
import numpy as np

def error(**kwargs):
    """
    Compute difference between initial and final solutions.
    This should vanish due to the periodic boundary conditions
    """

    claw = acoustics_1d.setup(disable_output=True,**kwargs)
    claw.run()

    # tests are done across the entire domain of q normally
    q0 = claw.frames[0].state.get_q_global()
    qfinal = claw.frames[claw.num_output_times].state.get_q_global()

    q0 = q0.reshape([-1])
    qfinal = qfinal.reshape([-1])
    dx = claw.solution.domain.grid.delta[0]
    diff = dx*np.sum(np.abs(qfinal-q0))
    return diff

class TestRegression:
    def test_python_classic(self):
        assert abs(error(kernel_language='Python',solver_type='classic')-0.001981)<1e-5

    def test_fortran_classic(self):
        assert abs(error(kernel_language='Fortran',solver_type='classic')-0.001981)<1e-5

    def test_sharpclaw(self):
        assert abs(error(kernel_language='Fortran',solver_type='sharpclaw')-0.001540)<1e-5
        assert abs(error(kernel_language='Fortran',solver_type='sharpclaw',weno_order=11)-0.000521)<1e-5
        assert abs(error(kernel_language='Fortran',solver_type='sharpclaw',time_integrator='SSPLMMk3')-0.001545)<1e-5

class TestAccuracy:

    def test_python_classic(self):
        assert abs(error(num_cells=2000,kernel_language='Python',solver_type='classic'))<1e-5

    def test_fortran_classic(self):
        assert abs(error(num_cells=2000,kernel_language='Fortran',solver_type='classic'))<1e-5
        assert abs(error(num_cells=4000,kernel_language='Fortran',solver_type='classic'))<2e-6

    def test_sharpclaw(self):
        assert abs(error(num_cells=2000,kernel_language='Fortran',solver_type='sharpclaw'))<1e-8
        assert abs(error(num_cells=2000,kernel_language='Fortran',solver_type='sharpclaw',weno_order=11))<1e-8
        assert abs(error(num_cells=2000,kernel_language='Fortran',solver_type='sharpclaw',time_integrator='SSPLMMk3'))<2e-8
