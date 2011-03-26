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
from pyclaw.solution import Solution, Grid, Dimension
from pyclaw.evolve.clawpack import ClawSolver1D

def qinit(grid):
    """Simple shallow water Riemann problem initial condition"""
    
    # Initial condition paramters
    # hl = grid.aux_global['hl']
    # hr = grid.aux_global['hr']
    # ul = grid.aux_global['ul']
    # ur = grid.aux_global['ur']
    # sloc = grid.aux_global['sloc']
    
    # Grid
    grid.empty_q()
    grid.q[:,0] = hl * (grid.p_center[0] <= sloc) + hr * (grid.p_center[0] > sloc)
    grid.q[:,1] = hl*ul * (grid.p_center[0] <= sloc) + hr*ur * (grid.p_center[0] > sloc)

    
# Problem specific data
sloc = 0.0
hl = 3.0
ul = 0.0
hr = 1.0
ur = 0.0

# Initialize grids and solutions
x = Dimension('x',-5.0,5.0,500)
grid = Grid(x)
grid.aux_global['g'] = 1.0
grid.aux_global['efix'] = False
grid.meqn = 2
grid.t = 0.0
qinit(grid)
init_solution = Solution(grid)

# Solver setup
solver = ClawSolver1D()
solver.dt = 0.1
solver.max_steps = 500
solver.set_riemann_solver('shallow_roe')
solver.order = 2
solver.mthlim = [3,3]

# Controller instantiation
claw = Controller()
claw.outdir = './output'
claw.keep_copy = True
claw.output_format = None
claw.nout = 4
claw.outstyle = 1
claw.tfinal = 2.0
claw.solutions['n'] = init_solution
claw.solver = solver

# Solve
status = claw.run()

# Plot
if claw.keep_copy:
    import matplotlib.pyplot as plt
    for n in xrange(0,4+1):
        sol = claw.frames[n]
        plt.figure(1)
        plt.subplot(3,2,n+1)
        plt.plot(x.center,sol.q[:,0])
        plt.axis([x.lower,x.upper,0.0,3.2])
        plt.title('Height t = %s' % sol.t)
        plt.figure(2)
        plt.subplot(3,2,n+1)
        plt.plot(x.center,sol.q[:,1])
        plt.axis([x.lower,x.upper,0.0,4.6])
        plt.title('Momentum t = %s' % sol.t)
    plt.show()
