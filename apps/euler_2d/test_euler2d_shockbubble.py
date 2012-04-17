def test_2d_euler_shockbubble():
    def verify_classic_shockbubble(test_density):
        import os
        from clawpack.pyclaw.util import check_diff
        import numpy as np
        """ verifies 2d euler shockbubble from a previously verified classic run """
        thisdir = os.path.dirname(__file__)
        expected_density = np.loadtxt(os.path.join(thisdir,'verify_shockbubble_classic.txt'))
        test_err = np.linalg.norm(expected_density-test_density)
        expected_err = 0
        return check_diff(expected_err, test_err, abstol=1e-12)

    from shockbubble import shockbubble
    from clawpack.pyclaw.util import gen_variants
    for test in gen_variants(shockbubble, verify_classic_shockbubble, python_kernel=False, solver_type='classic'):
        yield test
