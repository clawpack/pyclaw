def test_woodward_colella_blast():
    """ tests against expected sharpclaw results """
    import woodward_colella_blast
    from clawpack.pyclaw.util import test_app, check_diff

    def verify_woodward_colella_blast(controller):
        """ given an expected value, returns a verification function """
        import numpy as np
        import os

        test_solution = controller.solution.state.get_q_global()

        if test_solution != None:
            thisdir = os.path.dirname(__file__)
            expected_density = np.loadtxt(os.path.join(thisdir,'blast_regression_density.txt'))
            test_density = test_solution[0,:]
            test_err = np.linalg.norm(expected_density-test_density)
            return check_diff(0, test_err, abstol=1.e-4)


    return test_app(woodward_colella_blast.setup, verify_woodward_colella_blast, {})

if __name__=='__main__':
    test_woodward_colella_blast()
