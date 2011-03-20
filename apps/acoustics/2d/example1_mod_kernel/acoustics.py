

#!/usr/bin/python
#!/usr/bin/env python
# encoding: utf-8
"""
advection.py

Example python script for solving the 2d advection equation.
"""

import os, sys
    
import numpy as np
import math
from petsc4py import PETSc

from petclaw.grid import Dimension
from petclaw.grid import Grid
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
    q=np.empty([grid.meqn,len(x),len(y)], order = 'F')

    pi = 4.0*math.atan(1.0)
    width = 0.20
    
    
    for i in range(len(x)):
        for j in range(len(y)):
            r = math.sqrt(x[i]**2 + y[j]**2)
            if abs(r-0.50) <= width:
                pressure = 1.0 + math.cos(pi*(r - 0.50)/width)
            else:
                pressure = 0.0

            q[0,i,j] = pressure
            q[1,i,j] = 0.0
            q[2,i,j] = 0.0
                
                
    grid.q=q


# Initialize grids and solutions
from dimsp2 import cparam
x = Dimension('x',-1.0,1.0,100,mthbc_lower=1,mthbc_upper=1)
y = Dimension('y',-1.0,1.0,100,mthbc_lower=1,mthbc_upper=1)
grid = Grid([x,y])

rho = 1.0
bulk = 4.0
cc = math.sqrt(bulk/rho)
zz = rho*cc
grid.aux_global['rho']= rho
grid.aux_global['bulk']=bulk
grid.aux_global['zz']= zz
grid.aux_global['cc']=cc
cparam.rho = grid.aux_global['rho']
cparam.bulk = grid.aux_global['bulk']
cparam.zz = grid.aux_global['zz']
cparam.cc = grid.aux_global['cc']

grid.meqn = 3
grid.mbc = 2
grid.t = 0.0
qinit(grid)
init_solution = Solution(grid)

# Solver setup
solver = PetClawSolver2D(kernelsType = 'F')

solver.dt_initial = 0.001
solver.dt_variable=True
solver.dt_max = 1e+99
solver.cfl_max = 0.1
solver.cfl_desired = 0.1
solver.max_steps = 50000
solver.mwaves = 2
solver.order = 2
solver.order_trans = 2
solver.mthlim = [0,0]
solver.src_split = 0


useController = True# controller does not work in case of 2D yet
makePlot = True


if useController:

    # Controller instantiation
    claw = Controller()
    claw.outdir = './_output/'
    claw.keep_copy = True
    claw.nout = 10
    claw.outstyle = 1
    claw.output_format = 'petsc'
    claw.tfinal = .27
    claw.solutions['n'] = init_solution
    claw.solver = solver

    # Solve
    status = claw.run()

    if makePlot:
        if claw.keep_copy:
            for n in xrange(0,claw.nout+1):
                sol = claw.frames[n]
                plotTitle="time: {0}".format(sol.t)
                viewer = PETSc.Viewer.DRAW(sol.grid.gqVec.comm)
                OptDB = PETSc.Options()
                OptDB['draw_pause'] = 1
                viewer(sol.grid.gqVec)


else:
    sol = {"n":init_solution}
    
    solver.evolve_to_time(sol,.27)
    
    sol = sol["n"]

    if makePlot:
        viewer = PETSc.Viewer.DRAW(grid.gqVec.comm)
        OptDB = PETSc.Options()
        OptDB['draw_pause'] = 1
        viewer(grid.gqVec)
    
        


