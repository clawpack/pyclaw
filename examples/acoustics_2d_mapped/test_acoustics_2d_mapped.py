from . import acoustics_2d_inclusions
import numpy as np
from clawpack.pyclaw.util import check_diff
import os

def check_error(data_filename,**kwargs):

    claw = acoustics_2d_inclusions.setup(disable_output=True,**kwargs)
    claw.run()
    test_pressure = claw.frames[-1].q[0,:,:]

    thisdir = os.path.dirname(__file__)
    expected_pressure = np.load(os.path.join(thisdir, 'pressure.npz'))['arr_0']
    test_err = np.max(np.abs(expected_pressure[:].reshape(-1) - 
                             test_pressure[:].reshape(-1)))
    return check_diff(0, test_err, abstol=1e-7)


def test_classic():
    assert check_error('verify_classic.txt',solver_type='classic',
                       num_cells=100,num_output_times=1)==None
