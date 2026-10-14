"""
Tests for the SharpClaw polynomial reconstructions (lim_type 4, 5 and 6),
the lim_type names and the deprecated weno_order alias.
"""
import warnings

import numpy as np
import pytest

from clawpack import pyclaw, riemann

sharpclaw1 = pytest.importorskip("clawpack.pyclaw.sharpclaw.sharpclaw1")

POLY = {4: 'poly', 5: 'poly_dissipative', 6: 'poly_cell_centered'}
ORDERS = (4, 6, 8, 10)


def formal_order(lim_type, order):
    """Formal order of accuracy of the interface values."""
    return order + 1 if lim_type == 6 else order


def sine_averages(xc, h, kwave=1):
    """Exact cell averages of sin(2*pi*kwave*x)."""
    w = 2*np.pi*kwave
    return (np.cos(w*(xc - h/2)) - np.cos(w*(xc + h/2)))/(w*h)


def reconstruct(num_cells, lim_type, order):
    """Reconstruct interface values of sin(2*pi*x) on a periodic grid;
    return the max interface error and the max jump ql(i+1) - qr(i)."""
    num_ghost = order//2 + 1
    h = 1./num_cells
    xc = (np.arange(-num_ghost, num_cells + num_ghost) + 0.5)*h
    q = sine_averages(xc, h)[None, :]
    ql, qr = sharpclaw1.poly.poly_recon(q, num_cells, num_ghost, lim_type, order)
    real = slice(num_ghost, num_cells + num_ghost)
    err_l = np.abs(ql[0, real] - np.sin(2*np.pi*(xc[real] - h/2))).max()
    err_r = np.abs(qr[0, real] - np.sin(2*np.pi*(xc[real] + h/2))).max()
    jump = np.abs(ql[0, num_ghost+1:num_cells+num_ghost+1]
                  - qr[0, num_ghost:num_cells+num_ghost]).max()
    return max(err_l, err_r), jump


def advect(num_cells, lim_type, order, kwave=1, cfl=0.4, kernel_language='Fortran'):
    """Advect sin(2*pi*kwave*x) once around a periodic unit domain;
    return the L1 error and the ratio of final to initial L2 norm."""
    solver = pyclaw.SharpClawSolver1D(riemann.advection_1D)
    solver.kernel_language = kernel_language
    solver.lim_type = lim_type
    solver.reconstruction_order = order
    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic
    solver.cfl_desired = cfl
    solver.cfl_max = 1.05*cfl

    x = pyclaw.Dimension(0., 1., num_cells, name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain, 1)
    state.problem_data['u'] = 1.
    state.q[0, :] = sine_averages(x.centers, x.delta, kwave)

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.tfinal = 1.
    claw.num_output_times = 1
    claw.output_format = None
    claw.keep_copy = True
    claw.verbosity = 0
    claw.run()

    q0 = claw.frames[0].q[0]
    qf = claw.frames[-1].q[0]
    error = x.delta*np.abs(qf - q0).sum()
    energy_ratio = np.linalg.norm(qf)/np.linalg.norm(q0)
    return error, energy_ratio


class TestReconstruction:

    @pytest.mark.parametrize("lim_type", sorted(POLY))
    @pytest.mark.parametrize("order", ORDERS)
    def test_order_of_accuracy(self, lim_type, order):
        # Coarse grids for the high orders keep the errors above round-off
        grids = (16, 32) if order < 8 else (8, 16)
        coarse, fine = (reconstruct(n, lim_type, order)[0] for n in grids)
        observed = np.log2(coarse/fine)
        assert observed > formal_order(lim_type, order) - 0.4

    @pytest.mark.parametrize("order", ORDERS)
    def test_poly_is_single_valued_at_faces(self, order):
        """lim_type 4 gives the same value on both sides of each face,
        so the Riemann solver adds no dissipation."""
        assert reconstruct(32, 4, order)[1] < 1e-14

    @pytest.mark.parametrize("lim_type", (5, 6))
    def test_dissipative_variants_jump_at_faces(self, lim_type):
        assert reconstruct(32, lim_type, 4)[1] > 1e-6


class TestSolver:

    @pytest.mark.parametrize("lim_type", sorted(POLY))
    def test_convergence(self, lim_type):
        coarse, fine = (advect(n, lim_type, 4)[0] for n in (40, 80))
        assert np.log2(coarse/fine) > formal_order(lim_type, 4) - 0.3

    @pytest.mark.parametrize("lim_type", sorted(POLY))
    @pytest.mark.parametrize("order", (8, 10))
    def test_high_orders_are_stable(self, lim_type, order):
        error, energy_ratio = advect(40, lim_type, order)
        assert error < 1e-5
        assert energy_ratio <= 1. + 1e-12

    def test_poly_conserves_energy(self):
        """An unresolved mode (4 cells per wavelength) keeps its energy
        with 'poly' but is damped by the dissipative variants."""
        ratios = {name: advect(16, name, 6, kwave=4)[1] for name in POLY.values()}
        assert ratios['poly'] > 0.999
        assert ratios['poly_dissipative'] < ratios['poly_cell_centered'] < 0.5

    @pytest.mark.parametrize("order", (5, 7, 12))
    def test_rejects_invalid_order(self, order):
        with pytest.raises(ValueError):
            advect(20, 'poly', order)

    def test_rejects_python_kernels(self):
        with pytest.raises(NotImplementedError):
            advect(20, 'poly', 4, kernel_language='Python')


class TestLimTypeNames:

    @pytest.mark.parametrize("name, value", [
        ('godunov', -1), ('tvd', 1), ('weno', 2), ('weno5', 3),
        ('poly', 4), ('poly_dissipative', 5), ('poly_cell_centered', 6)])
    def test_name_maps_to_integer(self, name, value):
        solver = pyclaw.SharpClawSolver1D(riemann.advection_1D)
        solver.lim_type = name
        assert solver.lim_type == value
        solver.lim_type = name.upper()
        assert solver.lim_type == value
        solver.lim_type = value
        assert solver.lim_type == value

    def test_unknown_name_raises(self):
        solver = pyclaw.SharpClawSolver1D(riemann.advection_1D)
        with pytest.raises(ValueError):
            solver.lim_type = 'spline'

    @pytest.mark.parametrize("name, value", [('weno', 2), ('weno5', 3), ('poly', 4)])
    def test_name_and_integer_give_same_result(self, name, value):
        order = 4 if value == 4 else 5
        assert advect(20, name, order) == advect(20, value, order)


class TestWenoOrderAlias:

    def test_weno_order_maps_to_reconstruction_order(self):
        solver = pyclaw.SharpClawSolver1D(riemann.advection_1D)
        with pytest.warns(FutureWarning, match="reconstruction_order"):
            solver.weno_order = 9
        assert solver.reconstruction_order == 9
        with pytest.warns(FutureWarning):
            assert solver.weno_order == 9

    def test_reconstruction_order_does_not_warn(self):
        solver = pyclaw.SharpClawSolver1D(riemann.advection_1D)
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            solver.reconstruction_order = 7


def test_poly_2d():
    """The 2D solver reuses the 1D reconstruction along each direction."""
    solver = pyclaw.SharpClawSolver2D(riemann.advection_2D)
    solver.lim_type = 'poly'
    solver.reconstruction_order = 6
    solver.bc_lower = [pyclaw.BC.periodic]*2
    solver.bc_upper = [pyclaw.BC.periodic]*2
    solver.cfl_desired = 0.4
    solver.cfl_max = 0.42

    x = pyclaw.Dimension(0., 1., 24, name='x')
    y = pyclaw.Dimension(0., 1., 24, name='y')
    domain = pyclaw.Domain([x, y])
    state = pyclaw.State(domain, 1)
    state.problem_data['u'] = 1.
    state.problem_data['v'] = 1.
    X, Y = state.grid.p_centers
    h = x.delta
    state.q[0, ...] = sine_averages(X, h)*sine_averages(Y, h)

    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.tfinal = 1.
    claw.num_output_times = 1
    claw.output_format = None
    claw.keep_copy = True
    claw.verbosity = 0
    claw.run()

    assert solver.num_ghost == 4
    q0 = claw.frames[0].q[0]
    qf = claw.frames[-1].q[0]
    assert np.abs(qf - q0).max() < 1e-4
