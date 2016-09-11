#!/usr/bin/env python
# encoding: utf-8
r"""
Shock-tube problem
===================================

Solve the one-dimensional Euler equations for inviscid, compressible flow:

.. math::
    \rho_t + (\rho u)_x & = 0 \\
    (\rho u)_t + (\rho u^2 + p)_x & = 0 \\
    E_t + (u (E + p) )_x & = 0.

The fluid is an ideal gas, with pressure given by :math:`p=\rho (\gamma-1)e` where
e is internal energy.

This script runs a shock-tube problem.

This example also demonstrates:

 - how to use an arbitrary downwind Runge-Kutta method by providing the
   Shu-Osher array coefficients and the SSP coefficient of the method.
 - How to use a 2nd-order, two-stage downwind Runge-Kutta method.
"""

from __future__ import absolute_import
import numpy as np
from clawpack import riemann
from clawpack.riemann.euler_with_efix_1D_constants import *

gamma = 1.4 # Ratio of specific heats

# Coefficients of Downwind Runge-Kutta(4,4) method
sspcoeff = 0.6850160627358832
v = np.array([1., 0., 0., 0., 0.])
a = np.array([[0., 0., 0., 0., 0.],
              [0.671254015683971, 0., 0., 0., 0.],
              [0.270090108540258, 0.342508031367942, 0., 0., 0.],
              [2.71491162884274e-13, 0., 0.685016062735883, 0., 0.],
              [0.306177880587856, 0.176917664597762, 0.150130853210932, 0.114169343789314, 0.]])
a_tilde = np.array([[ 0., 0., 0., 0., 0.],
                    [0.328745984316029, 0., 0., 0., 0.],
                    [0.387401860091801, 0., 0., 0., 0.],
                    [0.0803604341607595, 0.234623503103086 , 0., 0., 0.],
                    [0.252604257814136, 0., 0., 0., 0.]])
c = np.array([0., 0.5, 0.5, 1.])

def setup(use_petsc=False, outdir='./_output', solver_type='sharpclaw',
          kernel_language='Fortran',time_integrator='DWRK',disable_output=False):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if kernel_language =='Python':
        rs = riemann.euler_1D_py.euler_hllc_1D
    elif kernel_language =='Fortran':
        rs = riemann.euler_with_efix_1D

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)
        solver.time_integrator = time_integrator
        solver.cfl_desired = 0.6
        solver.cfl_max = 0.7
        if solver.time_integrator == 'DWRK':
            solver.v, solver.a, solver.a_tilde, solver.c = v, a, a_tilde, c
            solver.sspcoeff = sspcoeff
    elif solver_type=='classic':
        solver = pyclaw.ClawSolver1D(rs)

    solver.kernel_language = kernel_language

    solver.bc_lower[0]=pyclaw.BC.extrap
    solver.bc_upper[0]=pyclaw.BC.extrap

    mx = 800;
    x = pyclaw.Dimension(-1.0,1.0,mx,name='x')
    domain = pyclaw.Domain([x])
    state = pyclaw.State(domain,num_eqn)

    state.problem_data['gamma'] = gamma

    x = state.grid.x.centers

    rho_l = 1.; rho_r = 1./8
    p_l = 1.; p_r = 0.1
    state.q[density ,:] = (x<0.)*rho_l + (x>=0.)*rho_r
    state.q[momentum,:] = 0.
    velocity = state.q[momentum,:]/state.q[density,:]
    pressure = (x<0.)*p_l + (x>=0.)*p_r
    state.q[energy  ,:] = pressure/(gamma - 1.) + 0.5 * state.q[density,:] * velocity**2

    claw = pyclaw.Controller()
    claw.tfinal = 0.4
    claw.solution = pyclaw.Solution(state,domain)
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
    plotaxes.title = 'Density'

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = density
    plotitem.kwargs = {'linewidth':3}
    
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    plotaxes.title = 'Energy'

    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = energy
    plotitem.kwargs = {'linewidth':3}
    
    return plotdata

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
