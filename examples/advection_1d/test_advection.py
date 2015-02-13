def test_1d_advection():
    """test_1d_advection

    tests against expected classic, sharpclaw, and high-order weno results """

    import advection_1d

    def verify_expected(expected):
        """ given an expected value, returns a verification function """
        def advection_verify(claw):
            from clawpack.pyclaw.util import check_diff
            import numpy as np

            q0=claw.frames[0].state.get_q_global()
            qfinal=claw.frames[claw.num_output_times].state.get_q_global()

            if q0 is not None and qfinal is not None:
                dx=claw.solution.domain.grid.delta[0]
                test = dx*np.linalg.norm(qfinal-q0,1)
                return check_diff(expected, test, reltol=1e-4)
            else:
                return
        return advection_verify

    from clawpack.pyclaw.util import gen_variants

    classic_tests = gen_variants(advection_1d.setup, verify_expected(3.203924e-04),
                                 kernel_languages=('Python','Fortran'),
                                 solver_type='classic', outdir=None)

    sharp_tests_rk  = gen_variants(advection_1d.setup, verify_expected(1.163605e-05),
                                 kernel_languages=('Python','Fortran'),
                                 solver_type='sharpclaw',time_integrator='SSP104', outdir=None)

    sharp_tests_lmm = gen_variants(advection_1d.setup, verify_expected(3.39682948116e-05),
                                 kernel_languages=('Python','Fortran'),
                                 solver_type='sharpclaw',time_integrator='SSPMS32', outdir=None)

    weno_tests = gen_variants(advection_1d.setup, verify_expected(7.489618e-06),
                                 kernel_languages=('Fortran',), solver_type='sharpclaw',
                                 time_integrator='SSP104', weno_order=17,
                                 outdir=None)

    from itertools import chain
    for test in chain(classic_tests, sharp_tests_rk, sharp_tests_lmm, weno_tests):
        yield test


if __name__=='__main__':
    for test in test_1d_advection():
        test()

    # monkey patch util.build_variant_arg_dicts and rerun tests
    import clawpack.pyclaw.util
    import clawpack.boxclaw.util
    clawpack.pyclaw.util.build_variant_arg_dicts = clawpack.boxclaw.util.boxlib_build_variant_arg_dicts

    for test in test_1d_advection():
        test()
