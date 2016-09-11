def test_1d_euler():
    """test_1d_euler

    Tests against expected Classic and SharpClaw results.

    Runs the shocktube test problem with
    a) Classic solver and the Python HLLC solver;
	b) SharpClaw solver and given downwind Runge-Kutta method.

    Runs the shocksine test problem with a total fluctuation solver, and
    WENO characteristic decomposition in SharpClaw by using a given Runge-Kutta method.

    Runs the Woodward-Colella blast wave problem with a total fluctuation solver, and
    TVD characteristic decomposition in SharpClaw.

    """

    import shocktube, shocksine, woodward_colella_blast

    def verify_expected(expected):
        """ binds the expected results to the verify_shocktube methods """
        def euler_1d_verify(controller):
            """ given expected results, returns a verification function """
            from clawpack.pyclaw.util import check_diff
            import numpy as np
            import os

            test_solution = controller.solution.state.get_q_global()

            if test_solution is not None:
                thisdir = os.path.dirname(__file__)
                expected_density = np.loadtxt(os.path.join(thisdir,expected))
                test_density = test_solution[0,:]
                test_err = np.linalg.norm(expected_density-test_density)
                return check_diff(expected_density,test_density,reltol=1.e-5,delta=controller.solution.grid.delta)
            else:
                return
        return euler_1d_verify

    from clawpack.pyclaw.util import gen_variants

    shocktube_classic_tests = gen_variants(shocktube.setup, verify_expected('shocktube_regression_density_classic.txt'),
        kernel_languages=('Python',), solver_type='classic', disable_output=True)

    shocktube_sharp_tests   = gen_variants(shocktube.setup, verify_expected('shocktube_regression_density_sharpclaw.txt'),
        kernel_languages=('Fortran',), solver_type='sharpclaw', time_integrator='DWRK', disable_output=True)

    shocksine_test = gen_variants(shocksine.setup, verify_expected('shocksine_regression_density.txt'),
        kernel_languages=('Fortran',), solver_type='sharpclaw', use_char_decomp=True, tfluct_solver=True)

    woodward_colella_blast_test = gen_variants(woodward_colella_blast.setup, verify_expected('blast_regression_density.txt'),
        kernel_languages=('Fortran',), solver_type='sharpclaw', tfluct_solver=True)

    from itertools import chain
    for test in chain(shocktube_classic_tests,shocktube_sharp_tests,shocksine_test,woodward_colella_blast_test):
        yield test


if __name__=="__main__":
    import nose
    nose.main()
