from __future__ import absolute_import
def test_acoustics_2d_variable():
    """Test variable-coefficient 2D acoustics"""

    from . import acoustics_2d_interface

    def verify_acoustics(controller, solver_type='classic'):
        """ Verifies 2d variable-coefficient acoustics from a previously verified classic run """
        import os
        from clawpack.pyclaw.util import check_diff
        import numpy as np

        state = controller.frames[controller.num_output_times].state
        dx, dy = controller.solution.domain.grid.delta
        test_q = state.get_q_global()

        if test_q is not None:
            thisdir = os.path.dirname(__file__)
            expected_pressure = np.loadtxt(os.path.join(thisdir,
                                               'pressure_%s.txt' % solver_type))
            test_pressure = test_q[0,:,:]
            #test_err = dx*dy*np.linalg.norm(expected_pressure-test_pressure)
            test_err = np.max(np.abs(expected_pressure[:].reshape(-1) - 
                                     test_pressure[:].reshape(-1)))
            return check_diff(0, test_err, abstol=1e-1)


    from clawpack.pyclaw.util import gen_variants

    verify_func = lambda controller: verify_acoustics(controller, solver_type='classic')
    classic_tests = gen_variants(acoustics_2d_interface.setup, verify_func,
                                 solver_type='classic', disable_output=True, 
                                 num_cells=(50, 50))

    verify_func = lambda controller: verify_acoustics(controller, solver_type='sharpclaw')
    sharp_tests_rk   = gen_variants(acoustics_2d_interface.setup, 
                                    verify_func,
                                    solver_type='sharpclaw', 
                                    time_integrator='SSP104', 
                                    disable_output=True, num_cells=(50, 50))

    sharp_tests_lmm   = gen_variants(acoustics_2d_interface.setup, 
                                     verify_func, lim_type=1,
                                     solver_type='sharpclaw', 
                                     time_integrator='SSPLMMk2', 
                                     disable_output=True,
                                     num_cells=(50, 50))

    from itertools import chain
    for test in chain(classic_tests, sharp_tests_rk, sharp_tests_lmm):
        yield test


if __name__=="__main__":
    import nose
    nose.main()
