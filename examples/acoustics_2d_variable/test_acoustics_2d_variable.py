def test_acoustics_2d_variable():
    """Test variable-coefficient 2D acoustics"""

    from acoustics import acoustics2D

    def verify_classic_acoustics(controller):
        import os
        from clawpack.pyclaw.util import check_diff
        import numpy as np
        """ Verifies 2d variable-coefficient acoustics from a previously verified classic run """

        state = controller.frames[controller.num_output_times].state
        dx, dy = controller.solution.domain.grid.delta
        test_q=state.get_q_global()

        if test_q != None:
            thisdir = os.path.dirname(__file__)
            expected_pressure = np.loadtxt(os.path.join(thisdir,'pressure_classic.txt'))
            test_pressure = test_q[0,:,:]
            #test_err = dx*dy*np.linalg.norm(expected_pressure-test_pressure)
            test_err = np.max(np.abs(expected_pressure[:]-test_pressure[:]))
            return check_diff(0, test_err, abstol=1e-1)


    from clawpack.pyclaw.util import gen_variants

    classic_tests = gen_variants(acoustics2D, verify_classic_acoustics,
                                 solver_type='classic', disable_output=True)

    sharp_tests   = gen_variants(acoustics2D, verify_classic_acoustics,
                                 solver_type='sharpclaw', disable_output=True)

    from itertools import chain
    for test in chain(classic_tests, sharp_tests):
        yield test
