#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(grid,xlower,ylower,dx,dy,hin,hout):
    # Initialize petsc Structures for q
    grid.init_q_petsc_structures()
    
    # Create an array with fortran native ordering
    x =grid.x.center
    y =grid.y.center
    Y,X = np.meshgrid(y,x)
 
    q=np.empty([grid.meqn,len(x),len(y)], order = 'F')


 #   print dx, dy, hin, hout

 #   for i in range(len(x)):
 #       for j in range(len(y)):
 #           r = np.sqrt(grid.x.center[i]**2 + grid.y.center[j]**2)
 #           if grid.x.center[i]<0.0: 
 #               q[0,i,j] = 2
 #           else:
 #               q[0,i,j] = 1

    q[0,:,:] = hin * (X <= 0.0) + hout * (X > 0.0)
    q[1,:,:] = 0.0
    q[2,:,:] = 0.0
    grid.q=q


#def cellave(xlower,ylower,dx,dy,wl)
#    """
#    This function is a translation from fortran 77 to python  of the Clawpack routine
#    written by Prof. LeVeque. See cellave.f in clawpack/2d/lib/
#    """
    
#    xx = array([xlow,xlow,xlow+dx,xlow+dx,xlow], dtype=float)
#    yy - array([ylow,ylow+dy,ylow+dy,ylow,ylow], dtype=float)


def shallow2D(iplot=False,petscPlot=False,useController=True,htmlplot=True):
    """
    Example python script for solving the 2d shallow water equations.
    """

    from petclaw.grid import Dimension
    from petclaw.grid import Grid
    from pyclaw.solution import Solution
    from petclaw.evolve.clawpack import PetClawSolver2D
    from pyclaw.controller import Controller
    from petclaw import plot

    # Create grid
    mx=100; my=100
    xlower = -2.5
    xupper = 2.5
    ylower = -2.5
    yupper = 2.5
    x = Dimension('x',xlower,xupper,mx,mthbc_lower=1,mthbc_upper=1)
    y = Dimension('y',ylower,yupper,my,mthbc_lower=1,mthbc_upper=1)
    grid = Grid([x,y])
   
    # Get cell's sizes needed to initialize the solution 
    dx = grid.d[0]
    dy = grid.d[1]


    # Define number of equations and BC
    grid.meqn = 3
    grid.mbc = 3

    # Define parameters for simulation
    #g = 1.0 # gravity
    #grid.aux_global['g'] = g
    #from dimsp2 import cparam
    #for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)

    # Set initial condition (initial solution)
    hin = 2.0
    hout = 1.0
    qinit(grid,xlower,ylower,dx,dy,hin,hout)
    initial_solution = Solution(grid)


    # Define solver and solver's parameters
    solver = PetClawSolver2D()
    #solver.order = 2
    #solver.order_trans = 2
    solver.cfl_max = 0.15
    solver.cfl_desired = 0.1
    solver.mwaves = 3
    solver.mathlim = [4]*solver.mwaves

    # Define controller and controller's parameters
    claw = Controller()
    claw.keep_copy = True
    claw.output_format = 'petsc' # The output format MUST be set to petsc!
    tfinal = 1.5
    claw.tfinal = tfinal
    claw.solutions['n'] = initial_solution
    claw.solver = solver

    # Solve problem
    status = claw.run()

    if htmlplot:  plot.plotHTML()
    if petscPlot: plot.plotPetsc(claw)
    if iplot:     plot.plotInteractive()


if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        shallow2D(*args,**kwargs)
    else: shallow2D()



    






