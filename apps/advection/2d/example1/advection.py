#!/usr/bin/python
#!/usr/bin/env python
# encoding: utf-8
"""
advection.py

Example python script for solving the 2d advection equation.
"""

import os, sys
    
import numpy as np
from petsc4py import PETSc

from petclaw.grid import PCDimension as Dimension
from petclaw.grid import PCGrid as Grid
from pyclaw.solution import Solution
from petclaw.evolve.petclaw import PetClawSolver2D
from pyclaw.controller import Controller



def qinit(grid):

    # Set initial conditions for q.
    # Sample scalar equation with data that is piecewise constant with
    # q = 1.0  if  0.1 < x < 0.6   and   0.1 < y < 0.6
    #     0.1  otherwise

    # Initialize petsc Structures for q
    grid.init_q_petsc_structures()

    
    # Create an array with fortran native ordering
    x =grid.x.center
    y =grid.y.center
    q=np.empty([len(x),len(y),grid.meqn])
    for i in range(len(x)):
        for j in range(len(y)):
            if x[i] > 0.1 and x[i] < 0.6 and y[j]>0.1 and y[j] < 0.6:
                q[i,j,:] = 1.0
            else:
                q[i,j,:] = 0.1
                
    grid.q=q


# Initialize grids and solutions
x = Dimension('x',0.0,1.0,50,mthbc_lower=2,mthbc_upper=2)
y = Dimension('y',0.0,1.0,50,mthbc_lower=2,mthbc_upper=2)
grid = Grid([x,y])
grid.aux_global['u']=0.5
grid.aux_global['v']=1.0

grid.meqn = 1
grid.mbc = 2
grid.t = 0.0
qinit(grid)
init_solution = Solution(grid)

# Solver setup
solver = PetClawSolver2D(kernelsType = 'F')

solver.dt = 0.016
solver.dt_variable=True
solver.max_steps = 500
solver.set_riemann_solver('advection')
solver.order = 2
solver.order_trans = 2
solver.mthlim = [3]
solver.src_split = 0


useController = False # controller does not work in case of 2D yet
makePlot = False


if useController:

    # Controller instantiation
    claw = Controller()
    claw.outdir = './_output/'
    claw.keep_copy = True
    claw.nout = 1
    claw.outstyle = 1
    claw.output_format = 'petsc'
    claw.tfinal =0.5
    claw.solutions['n'] = init_solution
    claw.solver = solver

    # Solve
    status = claw.run()

    if makePlot:
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
    
    solver.evolve_to_time(sol,2)
    
    sol = sol["n"]

    if makePlot:
        viewer = PETSc.Viewer.DRAW(grid.gqVec.comm)
        OptDB = PETSc.Options()
        OptDB['draw_pause'] = -1
        viewer(grid.gqVec)
    
        


