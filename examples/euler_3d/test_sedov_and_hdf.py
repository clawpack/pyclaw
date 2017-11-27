from __future__ import absolute_import
from __future__ import print_function
def test_sedov_and_hdf():
    """Test HDF I/O on Sedov 3D Euler application"""

    try:
        import h5py
        import scipy
    except ImportError:
        import nose
        raise nose.SkipTest

    from . import Sedov

    def verify_sedov(controller):
        import os
        from clawpack.pyclaw.util import check_diff
        from clawpack.pyclaw import Solution
        
        thisdir = os.path.dirname(__file__)
        verify_dir = os.path.join(thisdir,'./Sedov_regression')

        # Expected solution
        sol_expected = Solution()
        sol_expected.read(1,path=verify_dir,file_format='hdf',read_aux=False)
        assert sol_expected.t == 0.1
        expected_q = sol_expected.state.q

        # Test solution
        sol_test = Solution()
        sol_test.read(1,path=controller.outdir,
                        file_format=controller.output_format,
                        read_aux=False,
                        options=controller.output_options)
        test_q = sol_test.state.get_q_global()


        if test_q is not None:
            q_err = check_diff(expected_q, test_q, reltol=1e-6)
            if q_err is not None:
                return q_err
        else:
            return

    from clawpack.pyclaw.util import gen_variants
    tempdir = './_sedov_test_results'
    classic_tests = gen_variants(Sedov.setup, 
                                 verify_sedov, solver_type='classic',
                                 disable_petsc=True,
                                 outdir=tempdir, num_cells=(16, 16, 16),
                                 num_output_times=1)

    import shutil
    from itertools import chain
    try:
        for test in chain(classic_tests):
            yield test
    finally:
        ERROR_STR= """Error removing %(path)s, %(error)s """
        try:
            shutil.rmtree(tempdir )
        except OSError as xxx_todo_changeme:
            (errno, strerror) = xxx_todo_changeme.args
            print(ERROR_STR % {'path' : tempdir, 'error': strerror })


if __name__=="__main__":
    import nose
    nose.main()
