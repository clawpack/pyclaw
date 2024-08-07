def test_shocksine():
    """ Test shock-sine wave interaction (Euler 1D)"""
    from . import shocksine
    from clawpack.pyclaw.util import test_app, check_diff
    import os

    thisdir = os.path.dirname(__file__)
    
    def verify_shocksine(controller):
        """ given an expected value, returns a verification function """
        import numpy as np

        test_solution = controller.solution.state.get_q_global()

        if test_solution is not None:
            expected_density = np.loadtxt(os.path.join(thisdir,'shocksine_regression_density.txt'))
            test_density = test_solution[0,:]
            test_err = np.linalg.norm(expected_density-test_density)
            return check_diff(0, test_err, abstol=1.e-4)

    return test_app(shocksine.setup, verify_shocksine, {'outdir': os.path.join(thisdir,'./_output')})


