#!/usr/bin/env python
# encoding: utf-8

def test_1d_sill():
    """test_1d_sill

    tests against expected classic solution of shallow water equations over
    a sill."""

    import sill

    def verify_expected(expected):
        """ given an expected value, returns a verification function """
        def sill_verify(claw):
            from clawpack.pyclaw.util import check_diff
            import numpy as np

            q0 = claw.frames[0].state.get_q_global()
            qfinal = claw.frames[claw.num_output_times].state.get_q_global()

            if q0 != None and qfinal != None:
                dx = claw.solution.domain.grid.delta[0]
                test = dx*np.linalg.norm(qfinal-q0,1)
                return check_diff(expected, test, reltol=1e-4)
            else:
                return
        return sill_verify

    from clawpack.pyclaw.util import gen_variants

    import pdb; pdb.set_trace()
    classic_tests = gen_variants(sill.setup, verify_expected(3.203924e-04), 
                                             kernel_languages=("Python"),
                                             solver_type='classic',
                                             outdir=None)

    from itertools import chain
    for test in chain(classic_tests):
        yield test

if __name__=='__main__':
    import pdb; pdb.set_trace()
    test_1d_sill()
