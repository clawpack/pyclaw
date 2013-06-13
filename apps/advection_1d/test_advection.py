def test_1d_advection():
    """test_1d_advection

    tests against expected classic, sharpclaw, and high-order weno results """

    from advection import advection

    def verify_expected(expected):
        """ given an expected value, returns a verification function """
        def advection_verify(claw):
            from clawpack.pyclaw.util import check_diff
            import numpy as np

            q0=claw.frames[0].state.get_q_global()
            qfinal=claw.frames[claw.num_output_times].state.get_q_global()

            if q0 != None and qfinal != None:
                dx=claw.solution.domain.grid.delta[0]
                test = dx*np.linalg.norm(qfinal-q0,1)
                return check_diff(expected, test, reltol=1e-4)
            else:
                return
        return advection_verify

    from clawpack.pyclaw.util import gen_variants

    classic_tests = gen_variants(advection, verify_expected(3.203924e-04),
                                 kernel_languages=('Python','Fortran'),
                                 solver_type='classic', outdir=None)

    sharp_tests   = gen_variants(advection, verify_expected(1.163605e-05),
                                 kernel_languages=('Python','Fortran'),
                                 solver_type='sharpclaw', outdir=None)

    weno_tests    = gen_variants(advection, verify_expected(7.489618e-06),
                                 kernel_languages=('Fortran',),
                                 solver_type='sharpclaw', weno_order=17,
                                 outdir=None)

    from itertools import chain
    for test in chain(classic_tests, sharp_tests, weno_tests):
        yield test

if __name__=='__main__':
    test_1d_advection()
