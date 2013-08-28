def test_shocksine():
    """ tests against expected sharpclaw results """
    import shocksine

    def verify_shocksine(controller):
        """ given an expected value, returns a verification function """
        import numpy as np
        import os

        test_solution = controller.solution.state.get_q_global()

        if test_solution != None:
            thisdir = os.path.dirname(__file__)
            expected_density = np.loadtxt(os.path.join(thisdir,'shocksine_regression_density.txt'))
            test_density = test_solution[0,:]
            test_err = np.linalg.norm(expected_density-test_density)
            return check_diff(0, test_err, abstol=1.e-12)

    from clawpack.pyclaw.util import gen_variants

    for test in gen_variants(shocksine.setup, verify_shocksine,
                             solver_type='sharpclaw', outdir=None):
        yield test

if __name__=='__main__':
    test_shocksine()
