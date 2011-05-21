#!/usr/bin/env python
# encoding: utf-8
"""
py_example.py

Created by Kyle Mandli on 2008-08-25.
Copyright (c) 2008 University of Washington. All rights reserved.

Example python script for solving the 1d shallow water equations.

"""

import os

import numpy as np

from pyclaw.controller import Controller
#from pyclaw.solution import Grid, Dimension
from petclaw.grid import Grid, Dimension
from pyclaw.solution import Solution
#from pyclaw.evolve.sharpclaw import SharpClawSolver1D as ClawSolver
from petclaw.evolve.sharpclaw import SharpPetClawSolver1D as ClawSolver
from petclaw import plot

def qinit(grid):
    """Simple shallow water Riemann problem initial condition"""
    q=np.zeros([grid.meqn,grid.n[0]],order='F')
    q[0,:] = hl * (grid.p_center[0] <= sloc) + hr * (grid.p_center[0] > sloc)
    q[1,:] = hl*ul * (grid.p_center[0] <= sloc) + hr*ur * (grid.p_center[0] > sloc)
    grid.empty_q()
    grid.q=q

    
# Problem specific data
sloc = 0.0
hl = 3.0
ul = 0.0
hr = 1.0
ur = 0.0

# Initialize grids and solutions
x = Dimension('x',-5.0,5.0,500,mthbc_lower=1,mthbc_upper=1)
grid = Grid(x)
grid.aux_global['g'] = 1.0
grid.aux_global['efix'] = False
grid.meqn = 2
grid.t = 0.0
grid.init_q_petsc_structures()
qinit(grid)
init_solution = Solution(grid)

# Solver setup
solver = ClawSolver()
solver.kernel_language='Fortran'

solver.time_integrator='SSP33'
solver.lim_type=2
solver.char_decomp=0
mbc=3

solver.dt = 0.1
solver.max_steps = 500
solver.mwaves = 2
solver.set_riemann_solver('shallow_roe')
solver.order = 2
solver.mthlim = [3]*solver.mwaves
solver.cfl_max=0.5
solver.cfl_desired=0.45

# Controller instantiation
claw = Controller()
claw.outdir = './_output'
claw.keep_copy = False
claw.output_format = 'petsc'
#claw.output_format = 'ascii'
claw.nout = 5
claw.outstyle = 1
claw.tfinal = 2.0
claw.solutions['n'] = init_solution
claw.solver = solver

# Solve
status = claw.run()

iplot=True

if htmlplot:  plot.plotHTML(outdir=outdir,format=output_format)
if iplot:     plot.plotInteractive(outdir=outdir,format=output_format)

# Plot
if claw.keep_copy:
    import matplotlib.pyplot as plt
    for n in xrange(0,claw.nout+1):
        sol = claw.frames[n]
        plt.figure(1)
        plt.subplot(3,2,n+1)
        plt.plot(x.center,sol.q[0,:])
        plt.axis([x.lower,x.upper,0.0,3.2])
        plt.title('Height t = %s' % sol.t)
        plt.figure(2)
        plt.subplot(3,2,n+1)
        plt.plot(x.center,sol.q[1,:])
        plt.axis([x.lower,x.upper,0.0,4.6])
        plt.title('Momentum t = %s' % sol.t)
    plt.show()
