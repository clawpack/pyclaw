from __future__ import absolute_import
def test_2d_euler_shockbubble():
    """test_2d_euler_shockbubble"""
    def verify_classic_shockbubble(controller):
        """ verifies 2d euler shockbubble from a previously verified classic run """

        import os
        from clawpack.pyclaw.util import check_diff
        import numpy as np

        test_q=controller.solution.state.get_q_global()


        if test_q is not None:
            thisdir = os.path.dirname(__file__)
            expected_density = np.loadtxt(os.path.join(thisdir,'verify_shockbubble_classic.txt'))
            test_density = test_q[0,:,:]
            test_err = np.linalg.norm(expected_density-test_density)
            expected_err = 0
            return check_diff(expected_err, test_err, abstol=1e-12)

    try:
        import scipy
    except ImportError:
        from nose import SkipTest
        raise SkipTest("Unable to import scipy, is it installed?")

    from . import shock_bubble_interaction
    
    from clawpack.pyclaw.util import gen_variants
    for test in gen_variants(shock_bubble_interaction.setup, verify_classic_shockbubble, 
                             kernel_languages=('Fortran',),
                             solver_type='classic', disable_output=True,
                             mx=160, my=40, tfinal=0.2, num_output_times=1):
        yield test


if __name__=="__main__":
    import nose
    nose.main()
