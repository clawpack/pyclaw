from __future__ import absolute_import
def test_1d_advection():
    """test_1d_advection

    tests against expected classic results """

    from . import advection_1d_nonunif

    def verify_expected(expected):
        """ given an expected value, returns a verification function """
        def advection_verify(claw):
            from clawpack.pyclaw.util import check_diff
            import numpy as np

            q0=claw.frames[0].state.get_q_global()
            qfinal=claw.frames[claw.num_output_times].state.get_q_global()

            if q0 is not None and qfinal is not None:
                dx=claw.solution.domain.grid.delta[0]
                test = dx*np.linalg.norm(qfinal-q0,1)
                return check_diff(expected, test, reltol=1e-4)
            else:
                return
        return advection_verify

    def verify_expected_nonunif(expected):
        """ given an expected value, returns a verification function for the nonuniform advection ex"""
        from clawpack.pyclaw.util import check_diff
        import numpy as np

        def mapc2p_nonunif(xc):
            neg = -1*(xc < 0) + (xc > 0)
            xp = xc**2
            xp = neg*xp
            return xp

        def advection_nu_verify(claw):
            q0=claw.frames[0].state.get_q_global()
            qfinal=claw.frames[claw.num_output_times].state.get_q_global()
            if q0 is not None and qfinal is not None:
                dx=claw.solution.domain.grid.delta[0]
                grid1d = claw.frames[0].state.grid
                grid1d.mapc2p = mapc2p_nonunif
                nx = 100
                aux = np.zeros((1,nx))
                aux[0,:] = np.diff(grid1d.p_nodes)/np.diff(grid1d.x.nodes)
                test = abs(np.sum(dx*aux[0,:]*(qfinal-q0)))
                return check_diff(expected, test, reltol=1e-4)
            else:
                return
        return advection_nu_verify

    from clawpack.pyclaw.util import gen_variants

    classic_tests = gen_variants(advection_1d_nonunif.setup, verify_expected_nonunif(7.817097663620487e-17),
                                 kernel_languages=('Python','Fortran'),
                                 solver_type='classic', outdir=None)

    
    from itertools import chain
    for test in chain(classic_tests):
        yield test


if __name__=="__main__":
    import nose
    nose.main()
