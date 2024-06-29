from . import dam_break, sill
import numpy as np
import os

thisdir = os.path.dirname(__file__)
path_expected_sols = os.path.join(thisdir,'expected_sols.npy')
expected_sols_dict = np.load(path_expected_sols,allow_pickle=True).item()

def error(test_name,**kwargs):
    if "dambreak" in test_name:
        claw = dam_break.setup(**kwargs)
    elif "sill" in test_name:
        claw = sill.setup(**kwargs)
    else:
        raise ValueError(f"Test name {test_name} not recognized")
    claw.run()
    test_depth = claw.frames[claw.num_output_times].state.get_q_global()[0].reshape([-1])
    dx = claw.solution.domain.grid.delta[0]
    assert test_name in expected_sols_dict.keys(), f"Test name {test_name} not found in {path_expected_sols}"
    expected_depth = expected_sols_dict[test_name]
    diff_L1 = dx*np.sum(np.abs(test_depth-expected_depth))
    return diff_L1

class TestShallowWater1D:
    def test_dambreak(self):
        solver_types = ['classic','sharpclaw']
        riemann_solvers = ['hlle','roe']

        #Tests Fortran
        for solver_type in solver_types:
            for riemann_solver in riemann_solvers:
                test_name = f'dambreak_fortran_{solver_type}_{riemann_solver}'
                assert error(test_name=test_name, kernel_language="Fortran", solver_type=solver_type,
                            riemann_solver=riemann_solver, disable_output=True)<1e-6, f"Test {test_name} failed"
        #Tests Python
        for solver_type in solver_types:
            test_name = f'dambreak_python_{solver_type}'
            assert error(test_name=test_name, kernel_language="Python", solver_type=solver_type,
                         riemann_solver="hlle", disable_output=True)<1e-6, f"Test {test_name} failed"
        
    def test_sill(self):
        for kernel_language in ['Fortran','Python']:
            test_name = f'sill_{kernel_language}'
            assert error(test_name=test_name, kernel_language=kernel_language, outdir=None)<1e-6, f"Test {test_name} failed"
    
    