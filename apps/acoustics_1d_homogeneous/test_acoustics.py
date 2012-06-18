def test_1d_acoustics():
    """test_1d_acoustics

    tests against known classic, sharpclaw, and high-order weno results """

    from acoustics import acoustics

    def verify_expected(expected):
        """ binds the expected value to the acoustics_verify methods """
        def acoustics_verify(claw):
            from clawpack.pyclaw.util import check_diff
            import numpy as np

            # tests are done across the entire domain of q normally
            q0=claw.frames[0].state.get_q_global()
            qfinal=claw.frames[claw.num_output_times].state.get_q_global()

            # and q_global is only returned on process 0
            if q0 != None and qfinal != None:
                q0 = q0.reshape([-1])
                qfinal = qfinal.reshape([-1])
                dx=claw.solution.domain.grid.delta[0]
                test = dx*np.sum(np.abs(qfinal-q0))
                return check_diff(expected, test, abstol=1e-5)
            else:
                return
        return acoustics_verify

    from clawpack.pyclaw.util import gen_variants

    classic_tests = gen_variants(acoustics, verify_expected(0.00104856594174),
                                 python_kernel=True, solver_type='classic')

    sharp_tests   = gen_variants(acoustics, verify_expected(0.000298879563857),
                                 python_kernel=True, solver_type='sharpclaw')

    weno_tests    = gen_variants(acoustics, verify_expected(0.000153070447918),
                                 python_kernel=False, solver_type='sharpclaw',
                                 weno_order=17)

    from itertools import chain
    for test in chain(classic_tests, sharp_tests, weno_tests):
        yield test
