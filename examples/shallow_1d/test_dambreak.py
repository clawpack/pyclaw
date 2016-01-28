#!/usr/bin/env python
# encoding: utf-8

def test_1d_dambreak():
    """test_1d_dambreak

    Tests against expected classic solution of shallow water equations
    for dam-break problem."""

    import dam_break

    def verify_dambreak(controller):
        import numpy as np
        import os

        test_solution = controller.solution.state.get_q_global()

        if test_solution is not None:
            thisdir = os.path.dirname(__file__)
            expected_depth = np.loadtxt(os.path.join(thisdir,'dam_break_depth.txt'))
            test_depth = test_solution[0,:]
            discrepancy = np.linalg.norm(expected_depth - test_depth)
            assert discrepancy < 1.e-12
            return

    from clawpack.pyclaw.util import gen_variants
    classic_tests = gen_variants(dam_break.setup, verify_dambreak, 
                                             kernel_languages=["Python"],
                                             solver_type='classic',
                                             disable_output=True)

    from itertools import chain
    for test in chain(classic_tests):
        yield test

if __name__=='__main__':
    import nose
    nose.main()
