#!/usr/bin/env python
r"""
Regenerate the figures in figures/ used by README.md.

Run manually after changing the reconstructions; the PNG files are
committed so that readers need nothing installed. Needs only NumPy and
matplotlib::

    python make_figures.py

Colour roles, the same in every figure: blue is the left-face value ql(i),
orange the right-face value qr(i); the orders N = 4, 6, 8, 10 use a
light-to-dark blue ramp; the three reconstruction types get one panel
each rather than a colour.
"""
import os
from fractions import Fraction

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

import poly_weights as pw  # noqa: E402
from sine_reconstruction import interface_error, sine_averages  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'figures')

# Palette: categorical slots 1-2 and an ordinal blue ramp, validated
# (CVD and contrast) on the light surface below.
DPI = 150

SURFACE = '#fcfcfb'
INK = '#0b0b0b'
INK_2 = '#52514e'
MUTED = '#a3a29c'
GRID = '#e6e5e0'
QL = '#2a78d6'          # left-face value ql(i)
QR = '#eb6834'          # right-face value qr(i)
QL_FILL = '#cde2fb'
QR_FILL = '#fbd9c9'
ORDER_RAMP = {4: '#86b6ef', 6: '#3987e5', 8: '#1c5cab', 10: '#0d366b'}

TITLES = {4: "lim_type 4, 'poly'",
          5: "lim_type 5, 'poly_dissipative'",
          6: "lim_type 6, 'poly_cell_centered'"}

plt.rcParams.update({
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE, 'savefig.facecolor': SURFACE,
    'axes.edgecolor': MUTED, 'axes.labelcolor': INK_2, 'axes.titlecolor': INK,
    'xtick.color': INK_2, 'ytick.color': INK_2, 'text.color': INK,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': 0.8,
    'axes.spines.top': False, 'axes.spines.right': False,
    'font.size': 10, 'axes.titlesize': 10.5, 'axes.titleweight': 'bold',
    'lines.linewidth': 2, 'legend.frameon': False,
})


def save(fig, name):
    os.makedirs(OUT, exist_ok=True)
    path = os.path.join(OUT, name)
    # no Software tag, so the files do not change with the matplotlib version
    fig.savefig(path, dpi=DPI, metadata={'Software': None}, bbox_inches='tight')
    plt.close(fig)
    print('wrote', os.path.relpath(path, HERE))


# ----------------------------------------------------------------------------
def fig_averages_vs_points():
    """Cell averages are not centre values, and treating them as such
    caps the old weights at second order."""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 3.8), gridspec_kw={'width_ratios': [1, 1]})

    n = 4
    h = 1. / n
    x = np.linspace(0, 1, 400)
    avg = sine_averages(n)
    xc = (np.arange(n) + 0.5) * h
    ax1.plot(x, np.sin(2*np.pi*x), color=MUTED, lw=1.5, label='q(x) = sin 2πx')
    for j in range(n):
        ax1.add_patch(Rectangle((j*h, 0), h, avg[j], facecolor=QL_FILL, edgecolor=SURFACE, lw=2))
    ax1.plot(xc, avg, 's', color=QL, ms=6, label='cell average (what PyClaw stores)')
    ax1.plot(xc, np.sin(2*np.pi*xc), 'o', color=INK, ms=4, label='value at the cell centre')
    ax1.set_xlabel('x')
    ax1.set_xlim(0, 1)
    ax1.set_ylim(-1.15, 1.15)
    for j in range(n):
        ax1.annotate('', (xc[j], avg[j]), (xc[j], np.sin(2*np.pi*xc[j])),
                     arrowprops=dict(arrowstyle='<->', color=INK_2, lw=1, shrinkA=4, shrinkB=4))
    ax1.legend(loc='lower left', fontsize=8.5)
    ax1.set_title('Average ≠ centre value: gap ≈ h²q″/24', loc='left')

    cells_n = [8, 16, 32, 64, 128]
    for order in pw.ORDERS:
        cells = list(pw.stencil(6, order))
        w_pt = [float(v) for v in pw.point_weights(cells, Fraction(-1, 2))]
        errors = []
        for m in cells_n:
            hh = 1. / m
            xx = (np.arange(m) + 0.5) * hh
            ql = sum(wj * np.roll(sine_averages(m), -j) for j, wj in zip(cells, w_pt))
            errors.append(np.abs(ql - np.sin(2*np.pi*(xx - hh/2))).max())
        ax2.loglog(cells_n, errors, color=QR, lw=1.5, alpha=0.9)
        fv = np.array([interface_error(m, 6, order) for m in cells_n])
        keep = fv > 1e-13
        ax2.loglog(np.array(cells_n)[keep], fv[keep], color=ORDER_RAMP[order], marker='o', ms=4)
    ax2.text(140, 2e-5, 'point-value weights\n(old), all N:\nslope 2', color=QR, fontsize=8.5, va='center')
    for order in pw.ORDERS:
        ax2.plot([], [], color=ORDER_RAMP[order], marker='o', ms=4,
                 label='cell-average, N=%d: slope %d' % (order, order + 1))
    ax2.legend(loc='center left', bbox_to_anchor=(1.0, 0.45), fontsize=8)
    ax2.set_xlim(7, 600)
    ax2.set_ylim(1e-14, 1)
    ax2.set_xlabel('number of cells')
    ax2.set_ylabel('max interface error')
    ax2.set_title('Old vs new weights (2k+1-cell stencils)', loc='left')
    save(fig, 'averages_vs_points.png')


# ----------------------------------------------------------------------------
def fig_primitive_construction():
    """Interpolate the primitive at the faces, differentiate at the face."""
    cells = list(pw.stencil(6, 4))           # i-2 .. i+2
    q = lambda s: 1.0 + 0.5*np.sin(0.9*s + 0.4) + 0.08*s  # noqa: E731
    big_q = lambda s: 1.0*s - 0.5/0.9*np.cos(0.9*s + 0.4) + 0.04*s**2  # noqa: E731
    faces = np.arange(cells[0] - 0.5, cells[-1] + 1.0)
    avg = [big_q(j + 0.5) - big_q(j - 0.5) for j in cells]
    prim = np.concatenate([[0], np.cumsum(avg)])
    target = -0.5
    xs = np.linspace(faces[0] - 0.3, faces[-1] + 0.3, 400)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 6.2), sharex=True,
                                   gridspec_kw={'height_ratios': [1, 1.15]})
    for j, a in zip(cells, avg):
        ax1.add_patch(Rectangle((j - 0.5, 0), 1, a, facecolor=QL_FILL, edgecolor=SURFACE, lw=2))
        ax1.text(j, a + 0.05, r'$\bar{q}_{i%s}$' % ('' if j == 0 else '%+d' % j), ha='center', fontsize=10, color=INK_2)
    ax1.plot(xs, q(xs), color=MUTED, lw=1.5)
    ax1.plot([target], [q(target)], 'o', color=QL, ms=8, zorder=5)
    ax1.annotate('wanted: ql(i) = q(x$_{i-1/2}$)', (target, q(target)), (target - 2.3, 1.85),
                 color=QL, fontsize=9, arrowprops=dict(arrowstyle='->', color=QL, lw=1))
    ax1.set_ylim(0, 2.1)
    ax1.set_ylabel('q')
    ax1.set_title('1. Known: the cell averages over the stencil', loc='left')

    coeffs = np.polyfit(faces, prim, len(faces) - 1)
    ax2.plot(xs, np.polyval(coeffs, xs), color=INK_2, lw=1.5, label='interpolating polynomial P(x)')
    ax2.plot(faces, prim, 'o', color=INK, ms=6, label=r'Q at the faces = running sum of $h\,\bar{q}$')
    slope = np.polyval(np.polyder(coeffs), target)
    t = np.linspace(target - 1.1, target + 1.1, 2)
    ax2.plot(t, np.polyval(coeffs, target) + slope*(t - target), color=QL, lw=2.5,
             label="tangent at x$_{i-1/2}$: slope P′ = ql(i)")
    ax2.set_ylabel('Q(x) = ∫ q dx')
    ax2.set_title('2. Interpolate the primitive at the faces; its slope at the face is ql(i)', loc='left')
    ax2.legend(loc='upper left', fontsize=8.5)
    ax2.set_xticks(faces)
    ax2.set_xticklabels(['x$_{i%s}$' % ('-1/2' if f == -0.5 else '%+g' % f).replace('+0.5', '+1/2')
                         .replace('-1.5', '-3/2').replace('+1.5', '+3/2').replace('-2.5', '-5/2')
                         .replace('+2.5', '+5/2') for f in faces])
    for ax in (ax1, ax2):
        ax.axvline(target, color=QL, lw=0.8, ls=':')
    fig.text(0.01, -0.02, 'Shown for the 5-cell stencil of lim_type 6, N = 4. The weights do not depend '
             'on q: they follow from this construction once, exactly, in poly_weights.py.',
             fontsize=8.5, color=INK_2)
    save(fig, 'primitive_construction.png')


# ----------------------------------------------------------------------------
def draw_strip(ax, y, used, colour, fill, label, face=None, span=5, centre_cell=True):
    """Row of cells i-span..i+span; cells in `used` are filled."""
    for j in range(-span, span + 1):
        on = j in used
        ax.add_patch(Rectangle((j - 0.5, y - 0.35), 1, 0.7, facecolor=fill if on else SURFACE,
                               edgecolor=colour if on else GRID, lw=1.6 if on else 1))
        if on:
            ax.text(j, y, 'i' if j == 0 else 'i%+d' % j, ha='center', va='center', fontsize=7.5, color=INK)
    if centre_cell:
        ax.add_patch(Rectangle((-0.5, y - 0.35), 1, 0.7, fill=False, edgecolor=INK, lw=1.6))
    if face is not None:
        ax.plot([face, face], [y - 0.5, y + 0.5], color=colour, lw=3, solid_capstyle='butt')
    ax.text(-span - 0.8, y, label, ha='right', va='center', fontsize=9, color=colour)


def fig_stencils():
    """Which cells each type uses for ql(i) and qr(i)."""
    fig, axes = plt.subplots(3, 2, figsize=(11, 5.6))
    for row, lim_type in enumerate(pw.LIM_TYPES):
        for col, order in enumerate((4, 6)):
            ax = axes[row, col]
            cells = list(pw.stencil(lim_type, order))
            draw_strip(ax, 1, cells, QL, QL_FILL, 'ql(i)', face=-0.5)
            draw_strip(ax, 0, [-j for j in cells], QR, QR_FILL, 'qr(i)', face=0.5)
            ax.set_xlim(-8.6, 5.7)
            ax.set_ylim(-0.6, 1.6)
            ax.axis('off')
            ax.set_title('%s, N = %d: %d cells, order %d'
                         % (TITLES[lim_type], order, len(cells), pw.formal_order(lim_type, order)),
                         loc='left', fontsize=9.5)
    fig.text(0.01, -0.02, 'Bold outline: cell i. Thick bar: the face whose value is computed '
             '(x$_{i-1/2}$ for ql, x$_{i+1/2}$ for qr). qr(i) uses the mirror image of the ql(i) cells.',
             fontsize=8.5, color=INK_2)
    fig.tight_layout()
    save(fig, 'stencils.png')


def fig_face_states():
    """The two states the Riemann solver sees at one face."""
    fig, axes = plt.subplots(3, 1, figsize=(7.5, 5.6))
    notes = {4: 'same cells, same weights: one value, no jump, no dissipation',
             5: 'each side leans upwind: the states differ, the jump dissipates',
             6: 'each side centred on its own cell: the states differ, the jump dissipates'}
    for ax, lim_type in zip(axes, pw.LIM_TYPES):
        cells = list(pw.stencil(lim_type, 4))
        # face x_{i+1/2}: left state qr(i), right state ql(i+1)
        draw_strip(ax, 1, [-j for j in cells], QR, QR_FILL, 'left state qr(i)', face=0.5)
        draw_strip(ax, 0, [1 + j for j in cells], QL, QL_FILL, 'right state ql(i+1)', face=0.5,
                   centre_cell=False)
        ax.set_xlim(-11.5, 5.7)
        ax.set_ylim(-0.6, 1.6)
        ax.axis('off')
        ax.set_title('%s (N = 4): %s' % (TITLES[lim_type], notes[lim_type]), loc='left', fontsize=9)
    fig.text(0.01, -0.02, 'Both rows reconstruct the value at the same face x$_{i+1/2}$ (thick bar), '
             'from the left and from the right.', fontsize=8.5, color=INK_2)
    fig.tight_layout()
    save(fig, 'face_states.png')


# ----------------------------------------------------------------------------
def fig_sine_reconstruction():
    """ql and qr of sin(4 pi x) from 8 cell averages (4 per wavelength), N = 4."""
    n = 8
    h = 1. / n
    x = np.linspace(0, 1, 400)
    avg = sine_averages(n, kwave=2)
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), sharey=True)
    for ax, lim_type in zip(axes, pw.LIM_TYPES):
        ql, qr = pw.reconstruct(avg, lim_type, 4)
        for j in range(n):
            ax.add_patch(Rectangle((j*h, min(0, avg[j])), h, abs(avg[j]), facecolor='#f0efec',
                                   edgecolor=SURFACE, lw=2))
        ax.plot(x, np.sin(4*np.pi*x), color=MUTED, lw=1.5)
        faces = np.arange(n) * h
        ax.plot(faces, ql, '>', color=QL, ms=7, label='ql(i) at x$_{i-1/2}$')
        ax.plot(faces + h, qr, '<', color=QR, ms=7, label='qr(i) at x$_{i+1/2}$')
        jump = np.abs(np.roll(ql, -1) - qr).max()
        ax.set_title('%s, N = 4' % TITLES[lim_type], loc='left', fontsize=9.5)
        ax.text(0.03, -1.32, 'max |ql(i+1) − qr(i)| = %.1e' % jump, fontsize=8.5, color=INK_2)
        ax.set_xlabel('x')
        ax.set_xlim(0, 1)
        ax.set_ylim(-1.45, 1.3)
    axes[0].legend(loc='upper right', fontsize=8.5)
    fig.text(0.01, -0.04, 'Grey: exact sin 4πx and its 8 cell averages, a coarse 4 cells per wavelength. '
             'Each face gets a value from each neighbouring cell; with \'poly\' they coincide.',
             fontsize=8.5, color=INK_2)
    fig.tight_layout()
    save(fig, 'sine_reconstruction.png')


def fig_convergence():
    """Interface error against the number of cells."""
    cells_n = np.array([4, 6, 8, 12, 16, 24, 32, 48, 64, 96, 128])
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.8), sharey=True)
    for ax, lim_type in zip(axes, pw.LIM_TYPES):
        for order in pw.ORDERS:
            errors = np.array([interface_error(m, lim_type, order) for m in cells_n])
            keep = errors > 1e-13
            ax.loglog(cells_n[keep], errors[keep], color=ORDER_RAMP[order], marker='o', ms=3.5)
            last = np.flatnonzero(keep)[-1]
            ax.text(cells_n[last] * 1.12, errors[last], 'N=%d' % order, fontsize=8,
                    color=INK_2, va='center')
        ax.set_title(TITLES[lim_type], loc='left', fontsize=9.5)
        ax.set_xlabel('number of cells')
        ax.set_xlim(3.5, 250)
        ax.set_ylim(1e-14, 2)
    axes[0].set_ylabel('max interface error, sin 2πx')
    fig.text(0.01, -0.04, 'Slopes: N for \'poly\' and \'poly_dissipative\', N+1 for '
             '\'poly_cell_centered\'. Curves stop where round-off takes over.',
             fontsize=8.5, color=INK_2)
    fig.tight_layout()
    save(fig, 'convergence.png')


def fig_fourier():
    """Dissipation and dispersion for u_t + u_x = 0 (semi-discrete)."""
    theta = np.linspace(1e-4, np.pi, 600)
    fig, axes = plt.subplots(2, 3, figsize=(11, 6.2), sharex=True, sharey='row')
    for col, lim_type in enumerate(pw.LIM_TYPES):
        top, bottom = axes[0, col], axes[1, col]
        for order in pw.ORDERS:
            lam = pw.fourier_symbol(lim_type, order, theta)
            top.plot(theta, lam.real, color=ORDER_RAMP[order], label='N=%d' % order)
            bottom.plot(theta, -lam.imag, color=ORDER_RAMP[order], label='N=%d' % order)
        bottom.plot(theta, theta, color=MUTED, lw=1.2, ls='--', label='exact')
        top.set_title(TITLES[lim_type], loc='left', fontsize=9.5)
        if lim_type == 4:
            top.text(0.1, -0.35, 'Re λ = 0 for every N:\nno numerical dissipation', fontsize=9, color=INK_2)
        bottom.set_xlabel('θ = kh   (points per wavelength = 2π/θ)')
        bottom.set_xticks([0, np.pi/4, np.pi/2, 3*np.pi/4, np.pi])
        bottom.set_xticklabels(['0', 'π/4\n(8)', 'π/2\n(4)', '3π/4\n(2.7)', 'π\n(2)'])
    axes[0, 0].set_ylabel('Re λ  (damping rate)')
    axes[1, 0].set_ylabel('−Im λ  (modified wavenumber)')
    axes[0, 1].legend(loc='lower left', fontsize=8.5)
    axes[1, 0].legend(loc='upper left', fontsize=8.5)
    axes[0, 0].set_ylim(-2.8, 0.15)
    fig.text(0.01, -0.03, 'Mode exp(ijθ) with h = 1 and upwind flux evolves as exp(λt); exact: '
             'λ = −iθ. Re λ = 0 means no damping; −Im λ below the diagonal means waves lag.',
             fontsize=8.5, color=INK_2)
    fig.tight_layout()
    save(fig, 'fourier.png')


if __name__ == '__main__':
    fig_averages_vs_points()
    fig_primitive_construction()
    fig_stencils()
    fig_face_states()
    fig_sine_reconstruction()
    fig_convergence()
    fig_fourier()
