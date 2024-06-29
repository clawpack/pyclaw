from . import advection_1d
import numpy as np
from clawpack.pyclaw.util import check_diff
import os

def error(**kwargs):
    """
    Compute difference between initial and final solutions.
    This should vanish due to the periodic boundary conditions
    """

    claw = advection_1d.setup(outdir=None,**kwargs)
    claw.run()

    # tests are done across the entire domain of q normally
    q0 = claw.frames[0].state.get_q_global()
    qfinal = claw.frames[claw.num_output_times].state.get_q_global()

    q0 = q0.reshape([-1])
    qfinal = qfinal.reshape([-1])
    dx = claw.solution.domain.grid.delta[0]
    diff = dx*np.sum(np.abs(qfinal-q0))
    return diff


class TestAdvection1D:
    def test_python_classic(self):
        assert abs(error(kernel_language='Python',solver_type='classic')-3.203924e-04)<1e-4

    def test_fortran_classic(self):
        assert abs(error(kernel_language='Fortran',solver_type='classic')-3.203924e-04)<1e-4

    def test_python_sharpclaw(self):
        assert abs(error(kernel_language='Python',solver_type='sharpclaw')-1.163605e-05)<1e-4

    def test_fortran_sharpclaw(self):
        assert abs(error(kernel_language='Fortran',solver_type='sharpclaw')-1.163605e-05)<1e-4

    def test_sharpclaw_multistep(self):
        assert abs(error(kernel_language='Fortran',solver_type='sharpclaw',time_integrator='SSPLMMk3')-1.500727e-05)<1e-4

    def test_weno17(self):
        assert abs(error(kernel_language='Fortran',solver_type='sharpclaw',weno_order=17)-7.489618e-06)<1e-4

