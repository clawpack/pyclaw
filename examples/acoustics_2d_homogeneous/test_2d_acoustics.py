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
                return check_diff(expected_pressure, test_pressure, reltol=1e-3)
            else:
                return
        return verify

    from clawpack.pyclaw.util import gen_variants
    import acoustics_2d

    classic_tests = gen_variants(acoustics_2d.setup, verify_data('verify_classic.txt'),
                                 kernel_languages=('Fortran',), solver_type='classic', disable_output=True)

    sharp_tests   = gen_variants(acoustics_2d.setup, verify_data('verify_sharpclaw.txt'),
                                 kernel_languages=('Fortran',), solver_type='sharpclaw', disable_output=True)

    from itertools import chain
    for test in chain(classic_tests, sharp_tests):
        yield test
