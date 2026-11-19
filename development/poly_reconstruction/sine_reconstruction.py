#!/usr/bin/env python
r"""
Toy example: reconstruct sin(2 pi x) on a periodic grid of N cells with
every polynomial reconstruction and order, and report the interface error
and the observed order of accuracy.

The input is the exact cell averages of sin(2 pi x); the output is the
value at every cell face, compared with the exact sin(2 pi x_{i -/+ 1/2}).

Usage::

    python sine_reconstruction.py                 # N = 16, NumPy weights
    python sine_reconstruction.py --cells 32
    python sine_reconstruction.py --fortran       # also check the compiled
                                                  # Fortran (clawpack built)
"""
import argparse

import numpy as np

from poly_weights import LIM_TYPES, ORDERS, formal_order, reconstruct


def sine_averages(num_cells, kwave=1):
    """Exact averages of sin(2 pi kwave x) over the cells of [0, 1]."""
    h = 1. / num_cells
    xc = (np.arange(num_cells) + 0.5) * h
    w = 2 * np.pi * kwave
    return (np.cos(w * (xc - h/2)) - np.cos(w * (xc + h/2))) / (w * h)


def interface_error(num_cells, lim_type, order):
    """Max error of ql and qr against the exact face values."""
    h = 1. / num_cells
    xc = (np.arange(num_cells) + 0.5) * h
    ql, qr = reconstruct(sine_averages(num_cells), lim_type, order)
    return max(np.abs(ql - np.sin(2*np.pi*(xc - h/2))).max(),
               np.abs(qr - np.sin(2*np.pi*(xc + h/2))).max())


def fortran_difference(num_cells, lim_type, order):
    """Max difference between the NumPy and the compiled Fortran result."""
    from clawpack.pyclaw.sharpclaw import sharpclaw1
    ng = order // 2 + 1
    q = sine_averages(num_cells)
    padded = np.concatenate([q[-ng:], q, q[:ng]])[None, :]
    ql_f, qr_f = sharpclaw1.poly.poly_recon(padded, num_cells, ng, lim_type, order)
    ql, qr = reconstruct(q, lim_type, order)
    real = slice(ng, num_cells + ng)
    return max(np.abs(ql_f[0, real] - ql).max(), np.abs(qr_f[0, real] - qr).max())


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--cells', type=int, default=16, help='number of cells N (default 16)')
    parser.add_argument('--fortran', action='store_true',
                        help='compare with the compiled Fortran reconstruction')
    args = parser.parse_args(argv)
    n = args.cells

    print('Interface error for sin(2 pi x) from %d cell averages; the observed order'
          ' compares %d and %d cells.\n' % (n, n, 2*n))
    header = '%-20s %3s %7s %12s %9s' % ('lim_type', 'N', 'formal', 'max error', 'observed')
    if args.fortran:
        header += ' %14s' % '|NumPy-Fortran|'
    print(header)
    print('-' * len(header))
    for lim_type, name in LIM_TYPES.items():
        for order in ORDERS:
            e1 = interface_error(n, lim_type, order)
            e2 = interface_error(2*n, lim_type, order)
            line = '%-20s %3d %7d %12.3e %9.2f' % ('%d %s' % (lim_type, name), order,
                                                    formal_order(lim_type, order), e1,
                                                    np.log2(e1 / e2))
            if args.fortran:
                line += ' %14.1e' % fortran_difference(n, lim_type, order)
            print(line)
    print('\nObserved orders fall below the formal order once the error nears '
          'round-off (about 1e-14); use fewer cells to see the high orders.')


if __name__ == '__main__':
    main()
