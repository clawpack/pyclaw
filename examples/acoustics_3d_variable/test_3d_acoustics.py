import os
from itertools import chain
import numpy as np
from clawpack.pyclaw.util import check_diff
from . import acoustics_3d_interface

def check_error(**kwargs):

    claw = acoustics_3d_interface.setup(problem='homogeneous',mx=128,my=4,mz=4,disable_output=True,**kwargs)
    claw.run()

    test_pressure = claw.frames[-1].q[0,:,:,:]
    expected_pressure = claw.frames[0].q[0,:,:,:]
    return check_diff(expected_pressure, test_pressure, reltol=1e-1,
                        delta=claw.solution.grid.delta)


class TestAcoustics3D:
    def test_classic(self):
        assert check_error(solver_type='classic')==None

    def test_sharpclaw(self):
        assert check_error(solver_type='sharpclaw')==None
