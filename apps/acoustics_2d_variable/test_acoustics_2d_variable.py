def test_acoustics_2d_variable():
    """Test variable-coefficient 2D acoustics"""
    def verify_classic_acoustics(test_state):
        import os
        from clawpack.pyclaw.util import check_diff
        import numpy as np
        """ Verifies 2d variable-coefficient acoustics from a previously verified classic run """

        test_q=test_state.get_q_global()

        if test_q != None:
            thisdir = os.path.dirname(__file__)
            expected_pressure = np.loadtxt(os.path.join(thisdir,'pressure_classic2.txt'))
            test_pressure = test_q[0,:,:]
            test_err = np.linalg.norm(expected_density-test_density)
            expected_err = 0
            return check_diff(expected_err, test_err, abstol=1e-12)

    from acoustics import acoustics2D
    from clawpack.pyclaw.util import gen_variants
    for test in gen_variants(acoustics, verify_classic_acoustics, python_kernel=False, solver_type='classic'):
        yield test
