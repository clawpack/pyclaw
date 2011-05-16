#!/usr/bin/env python
# encoding: utf-8
"""
Example python script for solving the 1d variable-coefficient advection 
equation: q_t + u(x)q_x = 0.
"""

import numpy as np

from petclaw.grid import Grid, Dimension
from pyclaw.solution import Solution
from pyclaw.controller import Controller
from petclaw.evolve.clawpack import PetClawSolver1D
from petsc4py import PETSc

def qinit(grid):

    # Initilize petsc Structures for q
    grid.init_q_petsc_structures()
    
    # Initial Data parameters
    ic = 3
    beta = 100.
    gamma = 0.
    x0 = 0.3
    x1 = 0.7
    x2 = 0.9

    x =grid.x.center
    
    q=np.zeros([grid.meqn,len(x)],order='F')
    
    # Gaussian
    qg = np.exp(-beta * (x-x0)**2) * np.cos(gamma * (x - x0))
    # Step Function
    qs = (x > x1) * 1.0 - (x > x2) * 1.0
    
    if ic == 1: q[0,:] = qg
    elif ic == 2: q[0,:] = qs
    elif ic == 3: q[0,:] = qg + qs

    grid.q=q


def auxinit(grid):
    # Initilize petsc Structures for aux
    maux = 1
    xghost=grid.x.centerghost
    grid.aux=np.empty([maux,len(xghost)])
    grid.aux[0,:] = np.sin(2.*np.pi*xghost)+2
    

def vc_advection(kernel_language='Python',iplot=False):
# Initialize grids and solutions
    xlower=0.0; xupper=1.0; mx=100
    x = Dimension('x',xlower,xupper,mx,mthbc_lower=2,mthbc_upper=2,mbc=2)
    grid = Grid(x)
    grid.meqn = 1
    grid.t = 0.0
    qinit(grid)
    auxinit(grid)
    init_solution = Solution(grid)

# Solver setup
    solver = PetClawSolver1D()
    solver.kernel_language='Python'
    solver.dt = 0.9*grid.x.d/np.max(np.abs(grid.aux[0,:]))
    print grid.x.d, solver.dt
    solver.max_steps = 5000
    solver.set_riemann_solver('vc_advection')
    solver.order = 2
    solver.mthlim = 4
    solver.dt_variable = True
    solver.mwaves = 1

    use_controller = True

# Controller instantiation
    claw = Controller()
    claw.outdir = './_output'
    claw.nout = 10
    claw.outstyle = 1
    claw.output_format = 'petsc'
    claw.tfinal = 1.0
    claw.solutions['n'] = init_solution
    claw.solver = solver

# Solve
    status = claw.run()

    iplot=True

    if iplot:
        from petclaw import plot
        plot.plotInteractive(format=claw.output_format)

if __name__=="__main__":
    import sys
    from pyclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    output=vc_advection(*args,**kwargs)
