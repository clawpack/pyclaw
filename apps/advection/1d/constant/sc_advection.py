#!/usr/bin/python
#!/usr/bin/env python
# encoding: utf-8
"""
sc_advection.py

Example python script for solving the 1d advection equation.
"""

import os, sys



import numpy as np
from petsc4py import PETSc


from petclaw.grid import PCDimension as Dimension
from petclaw.grid import PCGrid as Grid
from pyclaw.solution import Solution
from petclaw.evolve.sharpclaw import SharpClawSolver1D
from pyclaw.controller import Controller



def qinit(grid):

    # Initialize petsc Structures for q
    grid.init_q_petsc_structures()
    
    # Initial Data parameters
    ic = 3
    beta = 100.
    gamma = 0.
    x0 = 0.3
    x1 = 0.7
    x2 = 0.9
    
    # Create an array with fortran native ordering
    x =grid.x.center
   
    q=np.zeros([len(x),grid.meqn])
    
    # Gaussian
    qg = np.exp(-beta * (x-x0)**2) * np.cos(gamma * (x - x0))

    # Step Function
    qs = (x > x1) * 1.0 - (x > x2) * 1.0
    
    if ic == 1:
        q[:,0] = qg
    elif ic == 2:
        q[:,0] = qs
    elif ic == 3:
        q[:,0] = qg + qs
    grid.q=q


# Initialize grids and solutions
x = Dimension('x',0.0,1.0,100,mthbc_lower=2,mthbc_upper=2)
grid = Grid(x)
grid.aux_global['u']=-1.
grid.meqn = 1
grid.t = 0.0
qinit(grid)
init_solution = Solution(grid)

# Solver setup
solver = SharpClawSolver1D(kernelsType = 'P')

solver.dt = 0.004
solver.dt_variable=False
solver.max_steps = 50000

solver.set_riemann_solver('advection')
solver.order = 2
solver.mthlim = 4
solver.lim_type = 2
solver.time_integrator='Euler'
solver.char_decomp=0

useController = True
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
    solver.evolve_to_time(sol,0.25)

    sol = sol["n"]

    if makePlot:
        viewer = PETSc.Viewer.DRAW(grid.gqVec.comm)
        OptDB = PETSc.Options()
        OptDB['draw_pause'] = -1
        viewer(grid.gqVec)
        


