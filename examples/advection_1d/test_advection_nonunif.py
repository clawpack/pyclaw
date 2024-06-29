from . import advection_1d_nonunif
import numpy as np
from clawpack.pyclaw.util import check_diff
import os

thisdir = os.path.dirname(__file__)
path_expected_sols = os.path.join(thisdir,'expected_sols.npy')
expected_sols_dict = np.load(path_expected_sols,allow_pickle=True).item()

def error(test_name,**kwargs):
    """
    Compute L1 norm of difference between test and expected
    solutions in the physical domain.
    """

    claw = advection_1d_nonunif.setup(outdir=None,**kwargs)
    claw.run()
    qtest = claw.frames[claw.num_output_times].state.get_q_global().reshape([-1])

    assert test_name in expected_sols_dict.keys(), f"Test name {test_name} not found in {path_expected_sols}"
    qexpected = expected_sols_dict[test_name]

    physical_nodes = claw.frames[0].state.grid.p_nodes
    dx=claw.solution.domain.grid.delta[0]
    comp_nodes = claw.frames[0].state.grid.c_nodes
    jac_mapc2p = np.diff(physical_nodes)/np.diff(comp_nodes)
    return np.sum(abs(dx*jac_mapc2p*(qtest-qexpected)))

class TestAdvectionNonUnif1D:
    def test_python_classic(self):
        assert error(test_name="python_classic",kernel_language='Python',solver_type='classic')<1e-6

    def test_fortran_classic(self):
        assert error(test_name="fortran_classic",kernel_language='Fortran',solver_type='classic')<1e-6
