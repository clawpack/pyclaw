from __future__ import absolute_import
def test_shocksine():
    """ tests against expected sharpclaw results """
    from . import shocksine
    from clawpack.pyclaw.util import test_app, check_diff

    def verify_shocksine(controller):
        """ given an expected value, returns a verification function """
        import numpy as np
        import os

        test_solution = controller.solution.state.get_q_global()

        if test_solution is not None:
            thisdir = os.path.dirname(__file__)
            expected_density = np.loadtxt(os.path.join(thisdir,'shocksine_regression_density.txt'))
            test_density = test_solution[0,:]
            test_err = np.linalg.norm(expected_density-test_density)
            return check_diff(0, test_err, abstol=1.e-4)

    return test_app(shocksine.setup, verify_shocksine, {})


if __name__=="__main__":
    import nose
    nose.main()
