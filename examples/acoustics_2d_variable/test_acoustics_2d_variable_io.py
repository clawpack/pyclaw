def test_acoustics_2d_variable_io():
    """Test I/O on variable-coefficient 2D acoustics application"""
    import os
    from . import acoustics_2d_interface
    thisdir = os.path.dirname(__file__)

    def verify_acoustics_io(controller):
        """ Verifies I/O on 2d variable-coefficient acoustics application"""
        from clawpack.pyclaw.util import check_diff
        from clawpack.pyclaw import Solution
        
        verify_dir = os.path.join(thisdir,'./io_test_verification')
        
        # Expected solution
        sol_0_expected = Solution()
        sol_0_expected.read(0,path=verify_dir,file_format='ascii',
                               file_prefix=None,read_aux=True)
        expected_aux = sol_0_expected.state.aux

        sol_20_expected = Solution()
        sol_20_expected.read(20,path=verify_dir,file_format='ascii',
                               file_prefix=None,read_aux=False)
        expected_q = sol_20_expected.state.q

        # Test solution
        sol_0_test = Solution()
        sol_0_test.read(0,path=controller.outdir,
                        file_format=controller.output_format,
                        file_prefix=None,read_aux=True,
                        options=controller.output_options)
        test_aux = sol_0_test.state.get_aux_global()

        sol_20_test = Solution()
        sol_20_test.read(20,path=controller.outdir,
                        file_format=controller.output_format,
                        file_prefix=None,read_aux=False,
                        options=controller.output_options)
        test_q = sol_20_test.state.get_q_global()


        if test_q is not None:
            q_err = check_diff(expected_q, test_q, reltol=1e-6)
            if q_err is not None:
                return q_err
        else:
            return

        if test_aux is not None:
            aux_err = check_diff(expected_aux, test_aux, reltol=1e-6)
            if aux_err is not None:
                return aux_err
        else:
            return
    
    tempdir = os.path.join(thisdir,'_io_test_results')
    claw = acoustics_2d_interface.setup(outdir=tempdir,num_cells=(50,50))
    claw.run()
    verify_acoustics_io(claw)
