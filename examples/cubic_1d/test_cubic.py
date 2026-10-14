"Runs the test problem."

def test_cubic():
    import numpy as np
    from . import cubic

    claw = cubic.setup(solver_type="sharpclaw")
    # For the classic solver with solver.limiters = pyclaw.limiters.tvd.vanleer,
    # nonclassical solutions occur...
    claw.run()
    sol = claw.solution.state.get_q_global()

    if sol is not None:
        xc = claw.solution.state.grid.x.centers
        qL = 4.0
        qR = -2.0
        x0 = -0.5 + (qL**3 - qR**3) / (qL - qR) * 0.2
        expected_sol = (xc < x0) * qL + (xc >= x0) * qR
        test_err = np.linalg.norm(expected_sol - sol, ord=1) / len(xc)
        print(test_err)
        assert test_err < 1.e-3


if __name__=="__main__":
    import nose
    nose.main()
