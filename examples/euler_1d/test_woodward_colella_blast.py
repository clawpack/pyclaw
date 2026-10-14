def test_woodward_colella_blast():
    """ Woodward-Colella blast wave test (Euler 1D)"""
    from . import woodward_colella_blast
    from clawpack.pyclaw.util import test_app, check_diff
    import os

    thisdir = os.path.dirname(__file__)

    def verify_woodward_colella_blast(controller):
        """ given an expected value, returns a verification function """
        import numpy as np

        test_solution = controller.solution.state.get_q_global()

        if test_solution is not None:
            expected_density = np.loadtxt(os.path.join(thisdir,'blast_regression_density.txt'))
            test_density = test_solution[0,:]
            return check_diff(expected_density, test_density, reltol=1.e-5,delta=controller.solution.grid.delta)

    return test_app(woodward_colella_blast.setup, verify_woodward_colella_blast, {'outdir': os.path.join(thisdir,'./_output')})


