from __future__ import absolute_import
def test_acoustics_2d_variable():
    """Test variable-coefficient 2D acoustics on mapped grids"""

    import acoustics_2d_inclusions

    def verify_acoustics_inclusions(controller, solver_type='classic'):
        """ Regression test against data from a previous run."""
        import os
        from clawpack.pyclaw.util import check_diff
        import numpy as np

        state = controller.frames[controller.num_output_times].state
        dx, dy = controller.solution.domain.grid.delta
        test_q = state.get_q_global()

        if test_q is not None:
            thisdir = os.path.dirname(__file__)
            expected_pressure = np.load(os.path.join(thisdir,
                                        'pressure.npz'))['arr_0']
            test_pressure = test_q[0,:,:]
            test_err = np.max(np.abs(expected_pressure[:].reshape(-1) - 
                                     test_pressure[:].reshape(-1)))
            return check_diff(0, test_err, abstol=1e-7)


    from clawpack.pyclaw.util import gen_variants

    verify_func = lambda controller: verify_acoustics_inclusions(controller, solver_type='classic')
    classic_tests = gen_variants(acoustics_2d_inclusions.setup, verify_func,
                                 solver_type='classic', disable_output=True, 
                                 num_cells=100,num_output_times=1)

    from itertools import chain
    for test in chain(classic_tests):
        yield test


if __name__=="__main__":
    import nose
    nose.main()
