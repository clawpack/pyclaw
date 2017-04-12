#!/usr/bin/env python
# encoding: utf-8

from __future__ import absolute_import
def test_2d_sill():
    """test_2d_sill

    Tests against expected classic solution of shallow water equations over
    a sill."""

    from . import sill

    def verify_expected(expected):
        def sill_verify(claw):
            from clawpack.pyclaw.util import check_diff
            import numpy as np

            q0 = claw.frames[0].state.get_q_global()
            qfinal = claw.frames[claw.num_output_times].state.get_q_global()

            if q0 is not None and qfinal is not None:
                dx, dy = claw.solution.domain.grid.delta
                total_mass = dx * dy * np.linalg.norm(qfinal[0,:,:].reshape(-1), 1)
                return check_diff(expected, total_mass, reltol=1e-3)
            else:
                return
        return sill_verify

    from clawpack.pyclaw.util import gen_variants
    classic_tests = gen_variants(sill.setup, verify_expected(3.7439),
                                 kernel_languages=["Fortran"],
                                 solver_type='classic', outdir=None)

    from itertools import chain
    for test in chain(classic_tests):
        yield test


if __name__=='__main__':
    import nose
    nose.main()
