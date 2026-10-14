"""
Toy tests for the polynomial reconstructions, using the exact weights from
poly_weights.py. Run from this directory with ``python -m pytest``.

CI does not collect development/; the Fortran and the generator are
covered there by src/pyclaw/tests/test_sharpclaw_poly.py.
"""
from fractions import Fraction

import numpy as np
import pytest

import poly_weights as pw
from sine_reconstruction import interface_error

SCHEMES = [(lim_type, order) for lim_type in pw.LIM_TYPES for order in pw.ORDERS]


@pytest.mark.parametrize("lim_type, order", SCHEMES)
def test_sine_order_of_accuracy(lim_type, order):
    grids = (16, 32) if order < 8 else (8, 16)
    coarse, fine = (interface_error(n, lim_type, order) for n in grids)
    assert np.log2(coarse / fine) > pw.formal_order(lim_type, order) - 0.4


@pytest.mark.parametrize("lim_type, order", SCHEMES)
def test_weights_are_consistent(lim_type, order):
    """Weights sum to one and reproduce averages of polynomials of the
    degree the stencil supports exactly."""
    cells, w = pw.ql_weights(lim_type, order)
    assert sum(w) == 1
    degree = len(cells) - 1
    for p in range(degree + 1):
        # average of x^p over cell j is ((j+1/2)^(p+1) - (j-1/2)^(p+1)) / (p+1)
        averages = [(Fraction(2*j + 1, 2)**(p + 1) - Fraction(2*j - 1, 2)**(p + 1)) / (p + 1)
                    for j in cells]
        assert sum(wj * a for wj, a in zip(w, averages)) == Fraction(-1, 2)**p


@pytest.mark.parametrize("order", pw.ORDERS)
def test_poly_is_non_dissipative(order):
    theta = np.linspace(0, np.pi, 201)
    assert np.abs(pw.fourier_symbol(4, order, theta).real).max() < 1e-13


@pytest.mark.parametrize("lim_type", (5, 6))
@pytest.mark.parametrize("order", pw.ORDERS)
def test_dissipative_variants_are_stable_and_dissipative(lim_type, order):
    theta = np.linspace(1e-3, np.pi, 201)
    real = pw.fourier_symbol(lim_type, order, theta).real
    assert real.max() < 1e-13
    assert real[-1] < -0.5


def test_point_value_weights_are_second_order():
    """The old weights treated averages as point values: second order."""
    cells = list(pw.stencil(6, 4))
    w = [float(v) for v in pw.point_weights(cells, Fraction(-1, 2))]
    errors = []
    for n in (16, 32):
        h = 1. / n
        xc = (np.arange(n) + 0.5) * h
        avg = (np.cos(2*np.pi*(xc - h/2)) - np.cos(2*np.pi*(xc + h/2))) / (2*np.pi*h)
        ql = sum(wj * np.roll(avg, -j) for j, wj in zip(cells, w))
        errors.append(np.abs(ql - np.sin(2*np.pi*(xc - h/2))).max())
    assert abs(np.log2(errors[0] / errors[1]) - 2) < 0.1


def test_generated_fortran_is_current():
    with open(pw.POLY_F90) as f:
        assert f.read() == pw.fortran_source()
