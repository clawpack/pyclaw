def test_shallow_sphere():
    """ test shallow sphere """

    try:
        import problem
        import classic2
    except ImportError:
        import warnings
        warnings.warn("missing extension modules, running python setup.py build_ext -i")
        import subprocess
        import os
        thisdir = os.path.dirname(__file__)
        subprocess.check_call('python setup.py build_ext -i', shell=True, cwd=thisdir)
        
    from shallow_4_Rossby_Haurwitz_wave import shallow_4_Rossby_Haurwitz

    def verify_shallow_sphere(claw):
        import os
        import numpy as np
        from clawpack.pyclaw.util import check_diff

        test_q = claw.frames[-1].state.get_q_global()
        test_height = test_q[0,:,:]

        thisdir = os.path.dirname(__file__)
        data_filename='swsphere_height.txt'
        expected_height = np.loadtxt(os.path.join(thisdir,data_filename))
        test_err = np.linalg.norm(expected_height-test_height)
        expected_err = 0
        return check_diff(expected_err, test_err, abstol=1e-4)

    from clawpack.pyclaw.util import test_app
    kwargs = {}   
    kwargs['disable_output']= True
    return test_app(shallow_4_Rossby_Haurwitz,
                    verify_shallow_sphere,
                    kwargs)

