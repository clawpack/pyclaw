from __future__ import absolute_import
def test_2d_acoustics():
    """test_2d_acoustics"""

    def verify_data(data_filename):
        def verify(claw):
            """ verifies 2d homogeneous acoustics from a previously verified run """
            import os
            import numpy as np
            from clawpack.pyclaw.util import check_diff


            #grabs parallel results to process 0, None to other processes
            test_q=claw.solution.state.get_q_global()

            if test_q is not None:
                test_pressure = test_q[0,:,:]
                thisdir = os.path.dirname(__file__)
                expected_pressure = np.loadtxt(os.path.join(thisdir,data_filename))
                return check_diff(expected_pressure, test_pressure, reltol=1e-3,
                                    delta=claw.solution.grid.delta)
            else:
                return
        return verify

    from clawpack.pyclaw.util import gen_variants
    from . import acoustics_2d

    classic_tests = gen_variants(acoustics_2d.setup, verify_data('verify_classic.txt'),
                                 kernel_languages=('Fortran',), solver_type='classic', disable_output=True)

    ptwise_tests = gen_variants(acoustics_2d.setup, verify_data('verify_classic.txt'),
                                 kernel_languages=('Fortran',), ptwise=True, solver_type='classic', disable_output=True)

    sharp_tests_rk   = gen_variants(acoustics_2d.setup, verify_data('verify_sharpclaw.txt'),
                                 kernel_languages=('Fortran',), solver_type='sharpclaw', 
                                 time_integrator='SSP104', disable_output=True)

    sharp_tests_lmm   = gen_variants(acoustics_2d.setup, verify_data('verify_sharpclaw_lmm.txt'),
                                 kernel_languages=('Fortran',), solver_type='sharpclaw', 
                                 time_integrator='SSPLMMk2', disable_output=True)

    from itertools import chain
    for test in chain(classic_tests, ptwise_tests, sharp_tests_rk, sharp_tests_lmm):
        yield test


if __name__=="__main__":
    import nose
    nose.main()
