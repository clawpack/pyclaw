from . import shock_bubble, shocktube, Sedov
import numpy as np
from pytest import importorskip
import os

thisdir = os.path.dirname(__file__)
path_expected_sols = os.path.join(thisdir,'expected_sols.npy')
expected_sols_dict = np.load(path_expected_sols,allow_pickle=True).item()

def error(test_name,**kwargs):
    if "shocktube" in test_name:
        claw = shocktube.shocktube(**kwargs)
    elif "shock_bubble" in test_name:
        claw = shock_bubble.setup(**kwargs)
    else:
        raise ValueError(f"Test name {test_name} not recognized")
    claw.run()
    density_test = claw.frames[claw.num_output_times].state.get_q_global()[0].reshape([-1])
    assert test_name in expected_sols_dict.keys(), f"Test name {test_name} not found in {path_expected_sols}"
    density_expected = expected_sols_dict[test_name]
    dx,dy,dz = claw.solution.domain.grid.delta
    diff_L1 = dx*dy*dz*np.sum(np.abs(density_test-density_expected))
    return diff_L1

class TestEuler3D:
    def test_shocktube(self):
        assert error(test_name='shocktube',mx=10,my=10,mz=128, tfinal=0.2,disable_output=True)<1e-6
    
    def test_shock_bubble(self):
        _ = importorskip("scipy") # skip test if scipy is not installed
        assert error(test_name='shock_bubble',num_cells=(64,16,16),tfinal=0.1,disable_output=True)<1e-6

    def test_sedov_and_hdf(self):
        _ = importorskip("scipy") # skip test if scipy is not installed
        _ = importorskip("h5py") # skip test if h5py is not installed
        from clawpack.pyclaw import Solution

        # Read expected solution
        sol_expected = Solution()
        sol_expected.read(1,path=os.path.join(thisdir,'Sedov_regression'),
                          file_format='hdf',read_aux=False, file_prefix='claw')
        density_expected = sol_expected.state.get_q_global()[0].reshape([-1])

        #Run test
        tempdir = os.path.join(thisdir,'_sedov_test_results')
        claw = Sedov.setup(solver_type='classic',
                                 outdir=tempdir, num_cells=(16, 16, 16),
                                 num_output_times=1, tfinal=0.1)
        claw.run()

        # Read test solution
        sol_test = Solution()
        sol_test.read(1,path=tempdir,file_format='hdf',read_aux=False,file_prefix='claw')
        density_test = sol_test.state.get_q_global()[0].reshape([-1])
        dx,dy,dz = sol_expected.domain.grid.delta
        diff_L1 = dx*dy*dz*np.sum(np.abs(density_test-density_expected))
        assert diff_L1<1e-6
    
