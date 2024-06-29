def test_advection_reaction():
    """ tests against expected sharpclaw results """
    from clawpack.pyclaw.examples.advection_reaction_2d import advection_reaction
    from clawpack.pyclaw.util import test_app, check_diff

    def verify_advection_reaction(controller):
        """ given an expected value, returns a verification function """
        import numpy as np
        import os

        test_solution = controller.solution.state.get_q_global()

        if test_solution is not None:
            thisdir = os.path.dirname(__file__)
            expected_density = np.loadtxt(os.path.join(thisdir,'advection_reaction.txt'))
            test_density = test_solution[0,:,:]
            return check_diff(expected_density, test_density, reltol=1.e-5,delta=controller.solution.grid.delta)

    return test_app(advection_reaction.setup, verify_advection_reaction, {"outdir":''})


# if __name__=="__main__":
#     import nose
#     nose.main()
