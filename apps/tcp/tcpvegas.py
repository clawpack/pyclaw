#!/usr/bin/env python
# encoding: utf-8
"""
tcpvegas.py

Python script for solving the 1d advection-reaction system modeling 
TCP vegas cwnd evolution.
We use a variable-coefficient formulation just to enforce the no-outflow
boundary conditions.
"""

import os

import numpy as np


def qinit(grid):

    x =grid.x.center
    
    grid.zeros_q()
    grid.q[0,0]=1./grid.x.d

def auxinit(grid):
    # Initilize petsc Structures for aux
    maux = 2
    xghost=grid.x.centerghost
    grid.zeros_aux(maux)
    grid.aux[:,0] =  1.
    grid.aux[:,1] = -1.
    grid.aux[-grid.mbc-1:,1] = 0. #No outflow at right boundary
    grid.aux[:grid.mbc+1,1]  = 0. #No outflow at left  boundary
    

import pyclaw
solver = pyclaw.ClawSolver1D()
solver.mwaves = 2
solver.kernel_language = 'Python'
solver.mthbc_lower[0] = pyclaw.BC.outflow
solver.mthbc_upper[0] = pyclaw.BC.outflow
import sys
sys.path.append('.')
import rp_tcp_advection
solver.rp = rp_tcp_advection.rp_tcp_advection_1d

# Initialize grids and solutions
alpha = 600.  #Low end of cwnd range
beta =  800.  #High end of cwnd range
xlower=alpha; xupper=beta; mx=200
x = pyclaw.Dimension('x',xlower,xupper,mx,mthbc_lower=1,mthbc_upper=1,mbc=2)
grid = pyclaw.Grid(x)
grid.meqn = 2
qinit(grid)
auxinit(grid)

#Transition probabilities:
grid.aux_global['lamda00'] = 0.5
grid.aux_global['lamda01'] = grid.aux_global['lamda00']
grid.aux_global['lamda10'] = 0.1
grid.aux_global['lamda11'] = grid.aux_global['lamda10']

claw = pyclaw.Controller()
claw.outstyle = 1
claw.nout = 100
claw.tfinal = 100.0
claw.solution = pyclaw.Solution(grid)
claw.solver = solver

# Solve
status = claw.run()

iplot=True
if iplot:     pyclaw.plot.interactive_plot()
if htmlplot:  pyclaw.plot.html_plot()
