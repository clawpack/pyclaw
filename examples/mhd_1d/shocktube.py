#!/usr/bin/env python
# encoding: utf-8
r"""
Shock-tube problem
===================================

Solve the one-dimensional ideal MHD equations for inviscid, compressible flow:

.. math::
    \rho_t + (\rho u_1)_x &= 0 \\
    (\rho u_i)_t + (\rho u_1 u_i + (p + B^2 / 2) \delta_{i1} - B_1 B_i)_x &= 0 \\
    E_t + (u_1 (E + p + B^2 / 2) - B_1 (u \cdot B))_x &= 0,
    (B_i)_t + (u_1 B_i - u_i B_1)_x &= 0.

The fluid is an ideal gas, with pressure given by :math:`p=\rho (\gamma-1)e` where
e is internal energy.

This script runs a shock-tube problem, cf. https://doi.org/10.1016/j.jcp.2016.10.034.
"""
from clawpack import riemann
from clawpack.riemann.mhd_1D_constants import num_eqn, density, momentum_1, momentum_2, momentum_3, energy, B_1, B_2, B_3

gamma = 1.4 # Ratio of specific heats

def setup(use_petsc=False, outdir='./_output', solver_type='classic',
          kernel_language='Fortran', disable_output=False):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if kernel_language =='Python':
        raise Exception('Not implemented.')
    elif kernel_language =='Fortran':
        rs = riemann.mhd_roe_1D

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)
    elif solver_type=='classic':
        solver = pyclaw.ClawSolver1D(rs)

    solver.kernel_language = kernel_language

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap

    mx = 800
    x = pyclaw.Dimension(-4.0, 4.0, mx, name='x')
    domain = pyclaw.Domain([x])
    state = pyclaw.State(domain, num_eqn)

    state.problem_data['gamma'] = gamma
    state.problem_data['gamma1'] = gamma - 1.

    x = state.grid.x.centers

    pressure = 1.0
    B1_l = 1.5; B1_r = 1.5
    B2_l = 0.5; B2_r = 1.6
    B3_l = 0.6; B3_r = 0.2
    state.q[density,    :] = 1.0
    state.q[momentum_1, :] = 0.0
    state.q[momentum_2, :] = 0.0
    state.q[momentum_3, :] = 0.0
    state.q[B_1,        :] = (x<0.) * B1_l + (x>=0.) * B1_r
    state.q[B_2,        :] = (x<0.) * B2_l + (x>=0.) * B2_r
    state.q[B_3,        :] = (x<0.) * B3_l + (x>=0.) * B3_r
    state.q[energy,     :] = pressure / (gamma - 1.) + 0.5 * (state.q[momentum_1,:]**2 + state.q[momentum_2,:]**2 + state.q[momentum_3,:]**2) / state.q[density,:] + 0.5 * (state.q[B_1,:]**2 + state.q[B_2,:]**2 + state.q[B_3,:]**2)

    claw = pyclaw.Controller()
    claw.tfinal = 1.0
    claw.solution = pyclaw.Solution(state, domain)
    claw.solver = solver
    claw.num_output_times = 10
    claw.outdir = outdir
    claw.setplot = setplot
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None

    return claw

#--------------------------
def setplot(plotdata):
#--------------------------
    """
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """
    plotdata.clearfigures()  # clear any old figures,axes,items data

    plotfigure = plotdata.new_plotfigure(name='', figno=0)

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.title = '$B_2$'

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = B_2
    plotitem.kwargs = {'linewidth': 3}

    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.title = '$B_3$'

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = B_3
    plotitem.kwargs = {'linewidth': 3}

    return plotdata

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup, setplot)
