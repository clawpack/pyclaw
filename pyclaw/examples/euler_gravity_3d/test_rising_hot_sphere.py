from __future__ import absolute_import
def test_euler_3d_rising_hot_sphere():
    """test_euler_3d_rising_hot_sphere"""
    def verify_classic_rising_hot_sphere(controller):
        """ verifies 3d euler rising hot sphere from a previously verified classic run """

        import os
        from clawpack.pyclaw.util import check_diff
        import numpy as np

        test_q=controller.solution.state.get_q_global()

        if test_q != None:
            thisdir = os.path.dirname(__file__)
            expected_density = np.loadtxt(os.path.join(thisdir,'verify_rising_hot_sphere_classic_1.txt'))
            nx = np.size(test_q,1)
            test_density = np.reshape(test_q[0,nx/2,:,:],np.size(test_q[0,nx/2,:,:]),order='F')
            test_err = np.linalg.norm(expected_density-test_density)
            expected_err = 0
            return check_diff(expected_err, test_err, abstol=1e-12)

    #main code
    try:
        import scipy
    except ImportError:
        from nose import SkipTest
        raise SkipTest("Unable to import scipy, is it installed?")

    import rising_hot_sphere
    
    from clawpack.pyclaw.util import gen_variants
    for test in gen_variants(rising_hot_sphere.euler3d, verify_classic_rising_hot_sphere,kernel_languages=('Fortran',),solver_type='classic', disable_output=True,mx=80, my=80, mz=80, tfinal=1.0, num_output_times=1):
        yield test
