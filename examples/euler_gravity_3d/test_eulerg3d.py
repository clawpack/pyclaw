import numpy as np
from pytest import importorskip
import os
import subprocess

thisdir = os.path.dirname(__file__)

class TestEulerGravity3D:

    def test_rising_hot_sphere(self):
        _ = importorskip("scipy") # skip test if scipy is not installed
        
        # Read expected solution
        expected_density = np.loadtxt(os.path.join(thisdir,'verify_rising_hot_sphere_classic_1.txt'))

        #Run test
        from . import rising_hot_sphere
        claw = rising_hot_sphere.euler3d(kernel_language='Fortran',solver_type='classic', 
                                         disable_output=True,mx=80, my=80, mz=80, tfinal=1.0, num_output_times=1)
        claw.run()
        test_q=claw.solution.state.get_q_global()
        nx = np.size(test_q,1)
        test_density = np.reshape(test_q[0,nx//2,:,:],np.size(test_q[0,nx//2,:,:]),order='F')
        test_err = np.linalg.norm(expected_density-test_density)
        assert test_err<1e-6
    
