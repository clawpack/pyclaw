"""
Tests for the SharpClaw polynomial reconstructions (lim_type 4, 5 and 6),
the lim_type names, the weno_order and poly_order properties and the
example setup() keywords.
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


@pytest.mark.filterwarnings("ignore:.*may deprecate:FutureWarning")
class TestOrderProperties:

    def make_solver(self, **kwargs):
        solver = pyclaw.SharpClawSolver1D(riemann.advection_1D)
        for key, value in kwargs.items():
            setattr(solver, key, value)
        return solver

    def test_weno_order_reads_by_lim_type(self):
        assert self.make_solver(reconstruction_order=9).weno_order == 9
        assert self.make_solver(lim_type='weno5').weno_order == 5
        assert self.make_solver(lim_type='poly', reconstruction_order=6).weno_order is None

    def test_weno_order_selects_weno(self):
        solver = self.make_solver(lim_type='poly', reconstruction_order=6)
        solver.weno_order = 7
        assert (solver.lim_type, solver.reconstruction_order) == (2, 7)
        solver = self.make_solver(lim_type='weno5')
        solver.weno_order = 5
        assert solver.lim_type == 3

    def test_poly_order_reads_by_lim_type(self):
        assert self.make_solver().poly_order is None
        for lim_type in POLY:
            assert self.make_solver(lim_type=lim_type, reconstruction_order=8).poly_order == 8

    def test_poly_order_selects_poly(self):
        solver = self.make_solver()
        solver.poly_order = 6
        assert (solver.lim_type, solver.reconstruction_order) == (4, 6)
        solver = self.make_solver(lim_type='poly_cell_centered')
        solver.poly_order = 10
        assert (solver.lim_type, solver.reconstruction_order) == (6, 10)

    @pytest.mark.parametrize("order", (5, 12))
    def test_poly_order_rejects_invalid_order(self, order):
        with pytest.raises(ValueError):
            self.make_solver().poly_order = order

@pytest.mark.filterwarnings("ignore:.*may deprecate:FutureWarning")
class TestSetReconstruction:

    def configure(self, **kwargs):
        solver = pyclaw.SharpClawSolver1D(riemann.advection_1D)
        solver.set_reconstruction(**kwargs)
        return solver.lim_type, solver.reconstruction_order

    @pytest.mark.parametrize("kwargs, expected", [
        ({}, (2, 5)),
        ({'reconstruction_order': 7}, (2, 7)),
        ({'weno_order': 11}, (2, 11)),
        ({'weno_order': 9, 'reconstruction_order': 9}, (2, 9)),
        ({'lim_type': 'weno5', 'weno_order': 5}, (3, 5)),
        ({'poly_order': 6}, (4, 6)),
        ({'lim_type': 'poly_dissipative', 'poly_order': 8}, (5, 8)),
        ({'lim_type': 6, 'reconstruction_order': 4}, (6, 4)),
        ])
    def test_valid(self, kwargs, expected):
        assert self.configure(**kwargs) == expected

    @pytest.mark.parametrize("kwargs", [
        {'weno_order': 5, 'poly_order': 4},
        {'weno_order': 7, 'reconstruction_order': 5},
        {'poly_order': 6, 'reconstruction_order': 4},
        {'lim_type': 'poly', 'weno_order': 5},
        {'lim_type': 'weno', 'poly_order': 4},
        {'lim_type': 'tvd', 'poly_order': 4},
        ])
    def test_contradictions_raise(self, kwargs):
        with pytest.raises(ValueError):
            self.configure(**kwargs)


@pytest.mark.filterwarnings("ignore:.*may deprecate:FutureWarning")
class TestExampleSetup:
    """Example setup() functions accept reconstruction_order, and still
    accept weno_order and poly_order."""

    def solver(self, **kwargs):
        from examples.acoustics_1d_homogeneous import acoustics_1d
        claw = acoustics_1d.setup(solver_type='sharpclaw', disable_output=True, **kwargs)
        return claw.solver

    def test_weno_order(self):
        solver = self.solver(weno_order=11)
        assert (solver.lim_type, solver.reconstruction_order) == (2, 11)

    def test_poly_order(self):
        solver = self.solver(poly_order=6)
        assert (solver.lim_type, solver.reconstruction_order) == (4, 6)
        solver = self.solver(poly_order=6, lim_type='poly_dissipative')
        assert (solver.lim_type, solver.reconstruction_order) == (5, 6)

    def test_defaults_unchanged(self):
        solver = self.solver()
        assert (solver.lim_type, solver.reconstruction_order) == (2, 5)

    def test_weno_order_matches_upstream_result(self):
        """setup(weno_order=11) reproduces the WENO11 regression value."""
        from examples.acoustics_1d_homogeneous.test_acoustics import error
        assert abs(error(kernel_language='Fortran', solver_type='sharpclaw',
                         weno_order=11) - 0.000521) < 1e-5

    @pytest.mark.parametrize("module, kwargs", [
        ('examples.advection_1d.advection_1d', {'kernel_language': 'Fortran'}),
        ('examples.advection_1d_variable.variable_coefficient_advection', {'kernel_language': 'Fortran'}),
        ('examples.cubic_1d.cubic', {'N': 50}),
        ])
    @pytest.mark.parametrize("order_kwarg", [{'weno_order': 7}, {'poly_order': 6}])
    def test_other_examples(self, module, kwargs, order_kwarg):
        import importlib
        setup = importlib.import_module(module).setup
        claw = setup(solver_type='sharpclaw', **kwargs, **order_kwarg)
        expected_lim_type = 2 if 'weno_order' in order_kwarg else 4
        assert claw.solver.lim_type == expected_lim_type
        assert claw.solver.reconstruction_order == list(order_kwarg.values())[0]

    def test_nonuniform_example_accepts_keywords(self):
        from examples.advection_1d import advection_1d_nonunif
        advection_1d_nonunif.setup(weno_order=5)
        advection_1d_nonunif.setup(poly_order=4)


class TestFutureWarnings:
    """weno_order and poly_order still work but warn that future releases
    may deprecate them; reconstruction_order never warns."""

    MATCH = "future releases may deprecate"

    def solver(self):
        return pyclaw.SharpClawSolver1D(riemann.advection_1D)

    @pytest.mark.parametrize("name, value", [('weno_order', 7), ('poly_order', 6)])
    def test_setting_warns(self, name, value):
        solver = self.solver()
        with pytest.warns(FutureWarning, match=self.MATCH):
            setattr(solver, name, value)

    @pytest.mark.parametrize("kwargs", [{'weno_order': 7}, {'poly_order': 6}])
    def test_set_reconstruction_warns(self, kwargs):
        with pytest.warns(FutureWarning, match=self.MATCH):
            self.solver().set_reconstruction(**kwargs)

    @pytest.mark.parametrize("kwargs", [{'weno_order': 7}, {'poly_order': 6}])
    def test_example_setup_warns_at_the_call(self, kwargs):
        from examples.acoustics_1d_homogeneous import acoustics_1d
        with pytest.warns(FutureWarning, match=self.MATCH) as record:
            acoustics_1d.setup(solver_type='sharpclaw', disable_output=True, **kwargs)
        assert len(record) == 1
        assert record[0].filename == __file__

    def test_reading_and_reconstruction_order_do_not_warn(self):
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            solver = self.solver()
            solver.reconstruction_order = 7
            solver.set_reconstruction(lim_type='poly', reconstruction_order=6)
            solver.weno_order, solver.poly_order

    @pytest.mark.parametrize("module, kwargs, order", [
        ('examples.acoustics_1d_homogeneous.acoustics_1d', {'disable_output': True}, 11),
        ('examples.advection_1d.advection_1d', {'kernel_language': 'Fortran'}, 7),
        ('examples.advection_1d_variable.variable_coefficient_advection',
         {'kernel_language': 'Fortran'}, 9),
        ('examples.cubic_1d.cubic', {'N': 50}, 7),
        ])
    def test_examples_take_reconstruction_order(self, module, kwargs, order):
        import importlib
        setup = importlib.import_module(module).setup
        with warnings.catch_warnings():
            warnings.simplefilter("error", FutureWarning)
            claw = setup(solver_type='sharpclaw', reconstruction_order=order, **kwargs)
            assert (claw.solver.lim_type, claw.solver.reconstruction_order) == (2, order)
            claw = setup(solver_type='sharpclaw', lim_type='poly', reconstruction_order=6, **kwargs)
            assert (claw.solver.lim_type, claw.solver.reconstruction_order) == (4, 6)

    def test_nonuniform_example_takes_reconstruction_order(self):
        from examples.advection_1d import advection_1d_nonunif
        with warnings.catch_warnings():
            warnings.simplefilter("error", FutureWarning)
            advection_1d_nonunif.setup(reconstruction_order=5)


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


@pytest.fixture(scope='module')
def generator():
    """development/poly_reconstruction/poly_weights.py, loaded by path."""
    import importlib.util
    import os
    path = os.path.join(os.path.dirname(__file__), '..', '..', '..', 'development',
                        'poly_reconstruction', 'poly_weights.py')
    if not os.path.exists(path):
        pytest.skip('development/ is not available')
    spec = importlib.util.spec_from_file_location('poly_weights', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TestGenerator:
    """poly.f90 is generated by development/poly_reconstruction/poly_weights.py;
    these tests stop the two from drifting apart."""

    def test_poly_f90_matches_generator(self, generator):
        with open(generator.POLY_F90) as f:
            assert f.read() == generator.fortran_source(), \
                'poly.f90 is out of date: run development/poly_reconstruction/poly_weights.py --write'

    @pytest.mark.parametrize("lim_type", sorted(POLY))
    @pytest.mark.parametrize("order", ORDERS)
    def test_compiled_weights_match_exact_weights(self, generator, lim_type, order):
        lo, hi, w, ok = sharpclaw1.poly.poly_weights(lim_type, order)
        cells, exact = generator.ql_weights(lim_type, order)
        assert ok and (lo, hi) == (cells[0], cells[-1])
        np.testing.assert_allclose(w[lo+5:hi+6], [float(v) for v in exact], rtol=0, atol=1e-15)
