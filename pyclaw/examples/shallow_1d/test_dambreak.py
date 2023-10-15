#!/usr/bin/env python
# encoding: utf-8

from __future__ import absolute_import
def test_1d_dambreak():
    """
    Tests solution against reference solution for shallow water dam break.
    """

    from . import dam_break
    import os
    

    def verify_expected():
        """ given an expected value, returns a verification function """
        def dambreak_verify(claw):
            from clawpack.pyclaw.util import check_diff
            import numpy as np

            thisdir = os.path.dirname(__file__)
            expected_depth = np.loadtxt(os.path.join(thisdir,'./dam_break_ref.txt'))

            test_solution = claw.frames[claw.num_output_times].state.get_q_global()

            if test_solution is not None:
                return check_diff(expected_depth, test_solution[0,:], reltol=1e-3)
            else:
                return
        return dambreak_verify

    from clawpack.pyclaw.util import gen_variants
    classic_tests = gen_variants(dam_break.setup, verify_expected(),
                                             riemann_solver="hlle",
                                             solver_type='classic',
                                             disable_output=True)

    from itertools import chain
    for test in chain(classic_tests):
        yield test

if __name__=='__main__':
    import nose
    nose.main()
