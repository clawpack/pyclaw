#!/usr/bin/env python
# encoding: utf-8
"""
advection.py

Example python script for solving the 1d advection equation.
"""

import os, sys



try:
    import numpy as np
    from petsc4py import PETSc
    
except:
    sys.path.append("/opt/share/ksl/petsc4py/dev-aug29/ppc450d/lib/python/")
    sys.path.append("/opt/share/ksl/numpy/dev-aug29/ppc450d/lib/python/")
    
    import numpy as np
    from petsc4py import PETSc



from petclaw.solution import Solution, Dimension, Grid
from petclaw.evolve.clawpack import ClawSolver1D
from petsc4py import PETSc
from petclaw.controller import Controller

def qinit(grid):

    # Initialize petsc Structures for q
    grid.init_q_petsc_structures()
    
    # Initial Data parameters
    ic = grid.aux_global['ic']
    beta = grid.aux_global['beta']
    gamma = grid.aux_global['gamma']
    x0 = grid.aux_global['x0']
    x1 = grid.aux_global['x1']
    x2 = grid.aux_global['x2']
    
    # Create an array with fortran native ordering
    x =grid.center(grid.x)
   
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

useController = True
makePlot = True


if useController:

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
        


