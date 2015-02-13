def test_1d_acoustics():
    """test_1d_acoustics

    tests against known classic, sharpclaw, and high-order weno results """

    import acoustics_1d

    def verify_expected(expected):
        """ binds the expected value to the acoustics_verify methods """
        def acoustics_verify(claw):
            from clawpack.pyclaw.util import check_diff
            import numpy as np

            # tests are done across the entire domain of q normally
            q0=claw.frames[0].state.get_q_global()
            qfinal=claw.frames[claw.num_output_times].state.get_q_global()

            # and q_global is only returned on process 0
            if q0 is not None and qfinal is not None:
                q0 = q0.reshape([-1])
                qfinal = qfinal.reshape([-1])
                dx=claw.solution.domain.grid.delta[0]
                test = dx*np.sum(np.abs(qfinal-q0))
                return check_diff(expected, test, abstol=1e-4)
            else:
                return
        return acoustics_verify

    from clawpack.pyclaw.util import gen_variants

    classic_tests = gen_variants(acoustics_1d.setup, verify_expected(0.00104856594174),
                                 kernel_languages=('Python','Fortran'), solver_type='classic', disable_output=True)

    ptwise_tests = gen_variants(acoustics_1d.setup, verify_expected(0.00104856594174),
                                 kernel_languages=('Fortran',), ptwise=True,
                                 solver_type='classic', disable_output=True)


    sharp_tests_rk   = gen_variants(acoustics_1d.setup, verify_expected(0.000298879563857),
                                 kernel_languages=('Python','Fortran'), solver_type='sharpclaw',
                                 time_integrator='SSP104', disable_output=True)

    sharp_tests_lmm   = gen_variants(acoustics_1d.setup, verify_expected(0.00227996627104),
                                 kernel_languages=('Python','Fortran'), solver_type='sharpclaw',
                                 time_integrator='SSPMS32', disable_output=True)

    weno_tests    = gen_variants(acoustics_1d.setup, verify_expected(0.000153070447918),
                                 kernel_languages=('Fortran',), solver_type='sharpclaw',
                                 time_integrator='SSP104', weno_order=17, disable_output=True)

    from itertools import chain
    for test in chain(classic_tests, ptwise_tests, sharp_tests_rk, sharp_tests_lmm, weno_tests):
        yield test


if __name__=='__main__':
    for test in test_1d_acoustics():
        test()

    # monkey patch util.build_variant_arg_dicts and rerun tests
    import clawpack.pyclaw.util
    import clawpack.boxclaw.util
    clawpack.pyclaw.util.build_variant_arg_dicts = clawpack.boxclaw.util.boxlib_build_variant_arg_dicts

    for test in test_1d_acoustics():
        test()
