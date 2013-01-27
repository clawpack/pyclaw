def test_acoustics_2d_variable_io():
    """Test I/O on variable-coefficient 2D acoustics application"""

    from acoustics import acoustics2D

    def verify_acoustics_io(controller):
        import os
        from clawpack.pyclaw.util import check_diff
        import numpy as np
        from  pyclaw import Solution
        """ Verifies I/O on 2d variable-coefficient acoustics application"""
        thisdir = os.path.dirname(__file__)
        verify_dir = os.path.join(thisdir,'./io_test_verification')
        
        sol_0 = Solution()
        sol_0.read(0,path=verify_dir,file_format='ascii',
                               file_prefix=None,read_aux=True)
        expected_aux = sol_0.state.aux

        sol_20 = Solution()
        sol_20.read(20,path=verify_dir,file_format='ascii',
                               file_prefix=None,read_aux=False)
        expected_q = sol_20.state.q


        state = controller.frames[controller.num_output_times].state
        test_q=state.get_q_global()
        test_aux = state.get_aux_global()
        
        test_passed = True
        if test_q is not None:
            q_err = check_diff(expected_q, test_q, reltol=1e-4)
            if q_err is not None:
                test_passed = False
        else:
            return

        if test_aux is not None:
            aux_err = check_diff(expected_aux, test_aux, reltol=1e-4)
            if aux_err is not None:
                test_passed = False
        else:
            return

        if test_passed:
            return None
        else:
            return ([q_err[0], aux_err[0]], [q_err[1], aux_err[1]]  ,q_err[2] )


    from clawpack.pyclaw.util import gen_variants

    classic_tests = gen_variants(acoustics2D, verify_acoustics_io,
                                 solver_type='classic')


    from itertools import chain
    for test in chain(classic_tests):
        yield test
