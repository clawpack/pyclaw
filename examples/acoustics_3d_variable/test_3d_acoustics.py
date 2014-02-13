def test_3d_acoustics():
    """ test_3d_acoustics """

    def acoustics_verify_homogeneous(claw):
        from clawpack.pyclaw.util import check_diff
        import numpy as np

        pinitial = claw.frames[0].state.get_q_global()
        pfinal   = claw.frames[claw.num_output_times].state.get_q_global()

        pinitial = pinitial[0,:,:,:].reshape(-1)
        pfinal   = pfinal[0,:,:,:].reshape(-1)
        grid = claw.solution.state.grid
        final_difference =np.prod(grid.delta)*np.linalg.norm(pfinal-pinitial,ord=1)

        return check_diff(0., final_difference, abstol=1e-1)

    def acoustics_verify_heterogeneous(claw):
        import os
        import numpy as np
        from clawpack.pyclaw.util import check_diff

        pinitial = claw.frames[0].state.get_q_global()
        pfinal   = claw.frames[claw.num_output_times].state.get_q_global()

        pinitial = pinitial[0,:,:,:].reshape(-1)
        pfinal   = pfinal[0,:,:,:].reshape(-1)
        grid = claw.solution.state.grid
        final_difference =np.prod(grid.delta)*np.linalg.norm(pfinal-pinitial,ord=1)

        thisdir = os.path.dirname(__file__)
        verify_pfinal = np.loadtxt(os.path.join(thisdir,'verify_classic_heterogeneous.txt'))
        norm_err = np.linalg.norm(pfinal-verify_pfinal)
        return check_diff(0, norm_err, abstol=10.)



    from clawpack.pyclaw.util import gen_variants
    import acoustics_3d_interface

    classichomogeneous_tests   = gen_variants(acoustics_3d_interface.setup, acoustics_verify_homogeneous,
                                       kernel_languages=('Fortran',), 
                                       solver_type='classic', test='homogeneous',
                                       disable_output=True)

    classicheterogeneous_tests = gen_variants(acoustics_3d_interface.setup, acoustics_verify_heterogeneous,
                                       kernel_languages=('Fortran',), 
                                       solver_type='classic', test='heterogeneous',
                                       disable_output=True)

    sharphomogeneous_tests   = gen_variants(acoustics_3d_interface.setup, acoustics_verify_homogeneous,
                                       kernel_languages=('Fortran',), 
                                       solver_type='sharpclaw', test='homogeneous',
                                       disable_output=True)

    sharpheterogeneous_tests = gen_variants(acoustics_3d_interface.setup, acoustics_verify_heterogeneous,
                                       kernel_languages=('Fortran',), 
                                       solver_type='sharpclaw', test='heterogeneous',
                                       disable_output=True)
    
    from itertools import chain
    for test in chain(classichomogeneous_tests, classicheterogeneous_tests, sharphomogeneous_tests, sharpheterogeneous_tests):
        yield test
