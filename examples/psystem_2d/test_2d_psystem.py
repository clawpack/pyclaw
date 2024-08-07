from . import psystem_2d
import numpy as np
import os

thisdir = os.path.dirname(__file__)


class TestPsystem2D:
    def test_psystem2d_gauges(self):
        tempdir = os.path.join(thisdir,'_for_temp_pyclaw_test')
        claw = psystem_2d.setup(kernel_language='Fortran',solver_type='classic', 
                                 outdir=tempdir, cells_per_layer=6, tfinal=1.)
        claw.run()
        test_state = claw.solution.state
        gauge_files = test_state.grid.gauge_files
        test_gauge_data_mem = test_state.gauge_data
        expected_gauges=[]
        for i, gauge in enumerate(gauge_files):
            test_gauge_data_io = np.loadtxt(gauge.name)
            verify_file = os.path.join(thisdir,'verify_' + gauge.name.split('/')[-1])
            expected_gauges.append(np.loadtxt(verify_file))
            diff_mem = np.allclose(expected_gauges[i], test_gauge_data_mem[i])
            diff_io = np.allclose(expected_gauges[i], test_gauge_data_io)
            assert diff_mem, f'Error with data in memory, gauge index: {i}.'
            assert diff_io, f'Error with data in file, gauge index: {i}.'
