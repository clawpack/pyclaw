from __future__ import absolute_import


def test_1d_acoustics():
    """test_1d_acoustics

    tests against known classic, sharpclaw, and high-order weno results """

    from . import acoustics_1d

    def verify_expected(expected):
        """ binds the expected value to the acoustics_verify methods """
        def acoustics_verify(claw):
            from clawpack.pyclaw.util import check_diff
            import numpy as np

            # tests are done across the entire domain of q normally
            q0 = claw.frames[0].state.get_q_global()
            qfinal = claw.frames[claw.num_output_times].state.get_q_global()

            # and q_global is only returned on process 0
            if q0 is not None and qfinal is not None:
                q0 = q0.reshape([-1])
                qfinal = qfinal.reshape([-1])
                dx = claw.solution.domain.grid.delta[0]
                test = dx*np.sum(np.abs(qfinal-q0))
                return check_diff(expected, test, abstol=1e-4)
            else:
                return
        return acoustics_verify

    from clawpack.pyclaw.util import gen_variants

    classic_tests = gen_variants(acoustics_1d.setup, verify_expected(0.001049),
                                 kernel_languages=('Python', 'Fortran'),
                                 solver_type='classic', disable_output=True)

    time_step_test = gen_variants(acoustics_1d.setup, verify_expected(0.002020),
                                  kernel_languages=('Python',),
                                  solver_type='classic', disable_output=True,
                                  output_style=(3))

    ptwise_tests = gen_variants(acoustics_1d.setup, verify_expected(0.001049),
                                kernel_languages=('Fortran',), ptwise=True,
                                solver_type='classic', disable_output=True)

    sharp_tests_rk = gen_variants(acoustics_1d.setup, verify_expected(0.000299),
                                  kernel_languages=('Python', 'Fortran'),
                                  solver_type='sharpclaw',
                                  time_integrator='SSP104', disable_output=True)

    sharp_tests_lmm = gen_variants(acoustics_1d.setup,
                                   verify_expected(0.000231),
                                   kernel_languages=('Python', 'Fortran'),
                                   solver_type='sharpclaw',
                                   time_integrator='SSPLMMk3',
                                   disable_output=True)

    weno_tests = gen_variants(acoustics_1d.setup, verify_expected(0.000153),
                              kernel_languages=('Fortran',),
                              solver_type='sharpclaw', time_integrator='SSP104',
                              weno_order=17, disable_output=True)

    from itertools import chain
    for test in chain(classic_tests, time_step_test, ptwise_tests,
                      sharp_tests_rk, sharp_tests_lmm, weno_tests):
        yield test


if __name__ == "__main__":
    import nose
    nose.main()
