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


from petclaw.grid import PCDimension as Dimension
from petclaw.grid import PCGrid as Grid
from pyclaw.solution import Solution
from pyclaw.controller import Controller
from petclaw.evolve.petclaw import PetClawSolver1D
from petsc4py import PETSc

def qinit(grid):

    # Initilize petsc Structures for q
    grid.init_q_petsc_structures()
    
    # Initial Data parameters
    
    x =grid.x.center
    
    q=np.zeros([len(x),grid.meqn],order='F')
    q[0,0]=1./grid.x.d
    
    grid.q=q


def auxinit(grid):
    # Initilize petsc Structures for aux
    maux = 2
    xghost=grid.x.centerghost
    grid.aux=np.empty([len(xghost),maux],order='F')
    grid.aux[:,0] =  1.
    grid.aux[:,1] = -1.
    grid.aux[-grid.mbc-1:,1] = 0. #No outflow at right boundary
    grid.aux[:grid.mbc+1,1]  = 0. #No outflow at left  boundary
    

# Initialize grids and solutions
alpha = 600.  #Low end of cwnd range
beta =  800.  #High end of cwnd range
xlower=alpha; xupper=beta; mx=200
x = Dimension('x',xlower,xupper,mx,mthbc_lower=1,mthbc_upper=1,mbc=2)
grid = Grid(x)
grid.meqn = 2
grid.t = 0.0
qinit(grid)
auxinit(grid)

#Transition probabilities:
grid.aux_global['lamda00'] = 0.5
grid.aux_global['lamda01'] = grid.aux_global['lamda00']
grid.aux_global['lamda10'] = 0.1
grid.aux_global['lamda11'] = grid.aux_global['lamda10']
init_solution = Solution(grid)

# Solver setup
solver = PetClawSolver1D(kernelsType = 'P')
solver.dt = 0.9*grid.x.d
print grid.x.d, solver.dt
solver.max_steps = 5000
solver.set_riemann_solver('tcp_advection')
solver.order = 2
solver.mthlim = [4,4]
solver.dt_variable = True

use_controller = True

if(use_controller):

# Controller instantiation
    claw = Controller()
    claw.outdir = './_output'
    claw.keep_copy = True
    claw.nout = 100
    claw.outstyle = 1
    claw.output_format = 'petsc'
    claw.tfinal = 100.0
    claw.solutions['n'] = init_solution
    claw.solver = solver

    # Solve
    status = claw.run()


    if claw.keep_copy:
    
        for n in xrange(0,11):
            sol = claw.frames[n]
            plotTitle="time: {0}".format(sol.t)
            viewer = PETSc.Viewer()
            viewer.createDraw(  title = plotTitle,  comm=sol.grid.gqVec.comm)


        
            OptDB = PETSc.Options()
            OptDB['draw_pause'] = -1
            sol.grid.gqVec.view(viewer)

else:
    sol = {"n":init_solution}
    
    solver.evolve_to_time(sol,.4)
    sol = sol["n"]

    #viewer = PETSc.Viewer.DRAW(grid.gqVec.comm)
    #OptDB = PETSc.Options()
    #OptDB['draw_pause'] = -1
    #viewer(grid.gqVec)
