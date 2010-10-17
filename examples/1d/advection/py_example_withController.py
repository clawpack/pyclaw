#!/usr/bin/env python
# encoding: utf-8
"""
py_example.py

Created by Kyle Mandli on 2008-08-25.
Copyright (c) 2008 University of Washington. All rights reserved.

Example python script for solving the 1d advection equation.

To see the equivalent example using the fortran clawpack, see fort_example.py.
"""

import os

import numpy as np

from pyclaw.controller import Controller
from pyclaw.solution import Solution, Dimension, Grid
from pyclaw.evolve.clawpack import ClawSolver1D
from petsc4py import PETSc

def qinit(grid):
    
    # Initial Data parameters
    ic = grid.aux_global['ic']
    beta = grid.aux_global['beta']
    gamma = grid.aux_global['gamma']
    x0 = grid.aux_global['x0']
    x1 = grid.aux_global['x1']
    x2 = grid.aux_global['x2']
    
    # Create an array with fortran native ordering
    x =grid.center(grid.x)
   
    grid.empty_q()
    
    # Grid
    # 
    # x = grid.x.center Amal: removed because dimention object does not know about the mbc
    #x = np.empty(grid.q.size)
    #for i in xrange(0,grid.q.size):
        #x[i] = grid.x.lower + (i+0.5)*grid.x.d
    
    # Gaussian
    qg = np.exp(-beta * (x-x0)**2) * np.cos(gamma * (x - x0))

    # Step Function
    qs = (x > x1) * 1.0 - (x > x2) * 1.0
    
    if ic == 1:
        grid.q[:,0] = qg
    elif ic == 2:
        grid.q[:,0] = qs
    elif ic == 3:
        grid.q[:,0] = qg + qs



    
   
           
    # why this is here not in the homogenous step? because it needs to be set once
    grid.gqVec.setArray(grid.q)
    
    
    

    
# Data paths and objects
example_path = './'
setprob_path = os.path.join(example_path,'setprob.data')

# Initialize grids and solutions
x = Dimension('x',0.0,1.0,100,mthbc_lower=2,mthbc_upper=2)
grid = Grid(x)
grid.set_aux_global(setprob_path)
grid.meqn = 1
grid.t = 0.0
qinit(grid)
init_solution = Solution(grid)

# Solver setup
solver = ClawSolver1D(kernelsType = 'P')

solver.dt = 0.0004
solver.max_steps = 5000
solver.set_riemann_solver('advection')
solver.order = 2
solver.mthlim = 4

# Controller instantiation
claw = Controller()
claw.outdir = './output/py_withController'
claw.keep_copy = True
claw.nout = 10
claw.outstyle = 1
claw.tfinal = 1.0
claw.solutions['n'] = init_solution
claw.solver = solver

# Solve
status = claw.run()


