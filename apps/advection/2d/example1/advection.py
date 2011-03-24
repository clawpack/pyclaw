#!/usr/bin/python
#!/usr/bin/env python
# encoding: utf-8

import numpy as np
from petsc4py import PETSc



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
    #q2=np.arange(len(q.reshape([-1])))
    for i in range(len(x)):
        for j in range(len(y)):
            if x[i] > 0.0 and x[i] < 0.5 and y[j]>0.0 and y[j] < 0.5:
                q[:,i,j] = 1.0
            else:
                q[:,i,j] = 0.1
                
    grid.q=q

def advection2D(iplot=False,petscPlot=False,useController=True,htmlplot=False):
    """
    Example python script for solving the 2d advection equation.
    """

    from petclaw.grid import Dimension
    from petclaw.grid import Grid
    from pyclaw.solution import Solution
    from petclaw.evolve.solver import PetClawSolver2D
    from pyclaw.controller import Controller
    from petclaw import plot

    mx=80; my=80
    # Initialize grids and solutions
    from dimsp2 import comrp
    x = Dimension('x',0.0,1.0,mx,mthbc_lower=1,mthbc_upper=1)
    y = Dimension('y',0.0,1.0,my,mthbc_lower=1,mthbc_upper=1)
    grid = Grid([x,y])

    grid.aux_global['u']=-0.6
    grid.aux_global['v']=0.4
    comrp.ubar = grid.aux_global['u']
    comrp.vbar = grid.aux_global['v']

    grid.meqn = 1
    qinit(grid)
    initial_solution = Solution(grid)

    solver = PetClawSolver2D()

    solver.dt = 0.016
    solver.dt_variable=False
    solver.max_steps = 5000
    solver.set_riemann_solver('advection')
    solver.order = 2
    solver.order_trans = 2
    solver.mwaves=1
    solver.mthlim = [3]

    claw = Controller()
    claw.keep_copy = True
    claw.nout = 10
    claw.output_format = 'petsc'
    claw.tfinal =solver.dt * 200
    claw.solutions['n'] = initial_solution
    claw.solver = solver

    #Solve
    status = claw.run()

    if htmlplot:  plot.plotHTML()
    if petscPlot: plot.plotPetsc(output_object)
    if iplot:     plot.plotInteractive()

if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        advection2D(*args,**kwargs)
    else: advection2D()
