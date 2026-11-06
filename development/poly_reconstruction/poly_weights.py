#!/usr/bin/env python
r"""
Exact weights for the SharpClaw polynomial reconstructions, and the
generator of ``src/pyclaw/sharpclaw/poly.f90``.

PyClaw stores cell averages. A reconstruction computes, for each cell i,
the value ql(i) at its left face x_{i-1/2} and qr(i) at its right face
x_{i+1/2}, as a linear combination of the averages of nearby cells:

    ql(i) = sum_j w_j q(i+j),        qr(i) = sum_j w_j q(i-j).

The weights are found by interpolating the primitive function
Q(x) = int_{x_lo}^{x} q dx, whose values at the faces are exact sums of
cell averages, and differentiating the interpolant at the face. They are
exact rationals; ``poly.f90`` stores each set as integers over a common
denominator.

lim_type selects the cells (offsets j) used for ql(i); with k = N/2 and
N = reconstruction_order in {4, 6, 8, 10}:

    4  'poly'                j = -k .. k-1  (N cells centred on the face)
    5  'poly_dissipative'    j = -k+1 .. k  (central stencil shifted
                                             one cell upwind)
    6  'poly_cell_centered'  j = -k .. k    (N+1 cells centred on the cell)

Usage::

    python poly_weights.py --table   # print the weights as Markdown tables
    python poly_weights.py --check   # exit 1 if poly.f90 is out of date
    python poly_weights.py --write   # regenerate poly.f90
"""
import argparse
import os
import sys
from fractions import Fraction
from math import lcm

import numpy as np

LIM_TYPES = {4: 'poly', 5: 'poly_dissipative', 6: 'poly_cell_centered'}
ORDERS = (4, 6, 8, 10)
DESCRIPTIONS = {4: 'central, non-dissipative',
                5: 'upwind-biased, dissipative',
                6: 'cell-centred, dissipative'}

HERE = os.path.dirname(os.path.abspath(__file__))
POLY_F90 = os.path.join(HERE, '..', '..', 'src', 'pyclaw', 'sharpclaw', 'poly.f90')


# ----------------------------------------------------------------------------
#   Weights
# ----------------------------------------------------------------------------
def formal_order(lim_type, order):
    """Formal order of accuracy of the interface values."""
    return order + 1 if lim_type == 6 else order


def stencil(lim_type, order):
    """Offsets j of the cells used for ql(i) (value at x_{i-1/2})."""
    if lim_type not in LIM_TYPES or order not in ORDERS:
        raise ValueError('lim_type must be in %s and order in %s'
                         % (tuple(LIM_TYPES), ORDERS))
    k = order // 2
    return {4: range(-k, k), 5: range(-k + 1, k + 1), 6: range(-k, k + 1)}[lim_type]


def face_weights(cells, x):
    """Weights giving the value at x from the averages of the unit cells
    centred at the integers in ``cells`` (cell j spans [j-1/2, j+1/2]).

    Interpolates the primitive function at the cell faces with a Lagrange
    polynomial and returns the weights of its derivative at x.
    """
    cells = list(cells)
    x = Fraction(x)
    faces = [Fraction(2*cells[0] - 1, 2) + m for m in range(len(cells) + 1)]
    # d/dx of the Lagrange basis polynomial of each face, at x
    dbasis = []
    for k, fk in enumerate(faces):
        total = Fraction(0)
        for a, fa in enumerate(faces):
            if a == k:
                continue
            term = Fraction(1) / (fk - fa)
            for b, fb in enumerate(faces):
                if b not in (a, k):
                    term *= (x - fb) / (fk - fb)
            total += term
        dbasis.append(total)
    # Q(face_k) = sum of the averages of the cells left of face k
    return [sum(dbasis[k] for k in range(len(faces)) if k > j)
            for j in range(len(cells))]


def point_weights(cells, x):
    """Lagrange weights that treat q as point values at the cell centres
    (the incorrect, second-order form used before; for comparison only)."""
    cells = list(cells)
    x = Fraction(x)
    weights = []
    for j in cells:
        w = Fraction(1)
        for m in cells:
            if m != j:
                w *= (x - m) / Fraction(j - m)
        weights.append(w)
    return weights


def ql_weights(lim_type, order):
    """(offsets, exact weights) for ql(i)."""
    cells = list(stencil(lim_type, order))
    return cells, face_weights(cells, Fraction(-1, 2))


def integer_weights(lim_type, order):
    """(lo, hi, numerators, denominator) for ql(i)."""
    cells, w = ql_weights(lim_type, order)
    denominator = lcm(*[v.denominator for v in w])
    numerators = [int(v * denominator) for v in w]
    return cells[0], cells[-1], numerators, denominator


# ----------------------------------------------------------------------------
#   Numerical helpers (used by the example, the tests and the figures)
# ----------------------------------------------------------------------------
def reconstruct(q, lim_type, order):
    """ql, qr for a periodic array of cell averages q (any length)."""
    cells, w = ql_weights(lim_type, order)
    w = [float(v) for v in w]
    ql = sum(wj * np.roll(q, -j) for j, wj in zip(cells, w))
    qr = sum(wj * np.roll(q, j) for j, wj in zip(cells, w))
    return ql, qr


def fourier_symbol(lim_type, order, theta):
    r"""Semi-discrete symbol lambda(theta) for u_t + u_x = 0 with h = 1 and
    the upwind flux, so that dq/dt = lambda q for q = exp(i j theta).
    The exact value is -i theta; Re(lambda) < 0 is dissipation."""
    cells, w = ql_weights(lim_type, order)
    theta = np.asarray(theta, dtype=float)
    # left state at face i+1/2 is qr(i) = sum_j w_j q(i-j)
    qr = sum(float(wj) * np.exp(-1j * j * theta) for j, wj in zip(cells, w))
    return -(1 - np.exp(-1j * theta)) * qr


# ----------------------------------------------------------------------------
#   Fortran generator
# ----------------------------------------------------------------------------
HEADER = """\
! ===================================================================
module poly
! ===================================================================
!   Linear (unlimited) polynomial reconstruction of interface values
!   from cell averages, for the SharpClaw finite-volume solvers.
!
!   q(:,i) holds the average of q over cell i. For each cell, ql(:,i) is
!   the reconstructed value at the left face x_(i-1/2) and qr(:,i) the
!   value at the right face x_(i+1/2). The weights are exact: they are
!   obtained by differentiating the polynomial that interpolates the
!   primitive function of q at the cell faces, so they are correct for
!   cell averages (not point values). They are applied component-wise.
!
!   lim_type selects the stencil; order = solver.reconstruction_order,
!   an even number N in {4, 6, 8, 10}. With k = N/2:
!
!     4 ('poly')               ql(i) uses cells i-k..i+k-1 (N cells,
!                              centred on the face). Order N. ql(i+1)
!                              equals qr(i), so the scheme adds no
!                              numerical dissipation.
!     5 ('poly_dissipative')   ql(i) uses cells i-k+1..i+k (the central
!                              stencil shifted one cell upwind for each
!                              face state). Order N, dissipative.
!     6 ('poly_cell_centered') ql(i) uses cells i-k..i+k (N+1 cells,
!                              centred on the cell, polynomial degree N).
!                              Order N+1, dissipative; equal to WENO of
!                              order N+1 with its ideal linear weights.
!
!   qr uses the mirror image of the ql stencil. Every variant reaches k
!   cells from the cell being reconstructed, so num_ghost = k + 1.
!
!   Coefficients generated with exact rational arithmetic.
! ===================================================================

    implicit none

contains

    ! ===================================================================
    subroutine poly_weights(lim_type, order, lo, hi, w, ok)
    ! ===================================================================
    !   Weights for the left-face value ql(i) = sum_j w(j) * q(i+j),
    !   j = lo..hi.

        integer, intent(in) :: lim_type, order
        integer, intent(out) :: lo, hi
        double precision, intent(out) :: w(-5:5)
        logical, intent(out) :: ok

        ok = .true.
        lo = 0; hi = 0
        w = 0.d0

        select case (lim_type)
"""

FOOTER = """\
        case default
            ok = .false.
        end select

    end subroutine poly_weights

    ! ===================================================================
    subroutine poly_recon(q,ql,qr,num_eqn,maxnx,num_ghost,lim_type,order)
    ! ===================================================================

        integer,          intent(in) :: num_eqn, maxnx, num_ghost, lim_type, order
        double precision, intent(in) :: q(num_eqn,maxnx+2*num_ghost)
        double precision, intent(out) :: ql(num_eqn,maxnx+2*num_ghost),qr(num_eqn,maxnx+2*num_ghost)

        double precision :: w(-5:5)
        integer :: lo, hi, reach, i, j, mx2
        logical :: ok

        call poly_weights(lim_type, order, lo, hi, w, ok)
        if (.not. ok) then
            print *, 'ERROR: polynomial reconstruction needs lim_type 4, 5 or 6'
            print *, '       and an even reconstruction_order between 4 and 10.'
            error stop
        endif

        mx2 = size(q,2)
        reach = max(-lo, hi)
        if (reach >= num_ghost) then
            print *, 'ERROR: polynomial reconstruction needs num_ghost >= ', reach+1
            error stop
        endif

        ql = 0.d0
        qr = 0.d0
        do i = 1+reach, mx2-reach
            do j = lo, hi
                ql(:,i) = ql(:,i) + w(j)*q(:,i+j)
                ! qr mirrors the ql stencil about cell i
                qr(:,i) = qr(:,i) + w(j)*q(:,i-j)
            end do
        end do

    end subroutine poly_recon

end module poly
"""


def _fortran_array(numerators, denominator):
    values = ['%d.d0' % v for v in numerators]
    if len(values) <= 8:
        return '[' + ', '.join(values) + '] / %d.d0' % denominator
    half = len(values) // 2
    return ('[' + ', '.join(values[:half]) + ', &\n' + ' '*22
            + ', '.join(values[half:]) + '] / %d.d0' % denominator)


def fortran_source():
    """The full text of poly.f90."""
    cases = []
    for lim_type in LIM_TYPES:
        orders = []
        for order in ORDERS:
            lo, hi, numerators, denominator = integer_weights(lim_type, order)
            orders.append('            case (%d)\n'
                          '                lo = %d; hi = %d\n'
                          '                w(%d:%d) = %s'
                          % (order, lo, hi, lo, hi, _fortran_array(numerators, denominator)))
        cases.append('        case (%d) ! %s\n'
                     '            select case (order)\n'
                     '%s\n'
                     '            case default\n'
                     '                ok = .false.\n'
                     '            end select'
                     % (lim_type, DESCRIPTIONS[lim_type], '\n'.join(orders)))
    return HEADER + '\n'.join(cases) + '\n' + FOOTER


def markdown_tables():
    """The weights of every scheme as Markdown tables."""
    out = []
    for lim_type, name in LIM_TYPES.items():
        out.append("**lim_type = %d, `'%s'`** (weights for ql(i), the value at "
                   "x<sub>i-1/2</sub>; qr(i) uses the mirror image)\n" % (lim_type, name))
        out.append('| N | order | cells | weights |')
        out.append('|---|---|---|---|')
        for order in ORDERS:
            lo, hi, numerators, denominator = integer_weights(lim_type, order)
            cells = 'i%+d … i%+d' % (lo, hi) if lo else 'i … i%+d' % hi
            out.append('| %d | %d | %s | (%s) / %d |'
                       % (order, formal_order(lim_type, order), cells,
                          ', '.join(str(v) for v in numerators), denominator))
        out.append('')
    return '\n'.join(out)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--write', action='store_true', help='regenerate poly.f90')
    group.add_argument('--check', action='store_true',
                       help='exit with status 1 if poly.f90 is out of date')
    group.add_argument('--table', action='store_true', help='print Markdown tables')
    args = parser.parse_args(argv)

    if args.table:
        print(markdown_tables())
        return 0
    if args.write:
        with open(POLY_F90, 'w') as f:
            f.write(fortran_source())
        print('Wrote %s' % os.path.normpath(POLY_F90))
        return 0
    with open(POLY_F90) as f:
        current = f.read()
    if current == fortran_source():
        print('poly.f90 is up to date.')
        return 0
    print('poly.f90 differs from the generator; run with --write.')
    return 1


if __name__ == '__main__':
    sys.exit(main())
