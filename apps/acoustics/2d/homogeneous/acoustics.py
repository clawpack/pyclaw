#!/usr/bin/python
#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(grid,width=0.2):
    # Initialize petsc Structures for q
    grid.init_q_petsc_structures()
    
    # Create an array with fortran native ordering
    x =grid.x.center
    y =grid.y.center
    X,Y = np.meshgrid(x,y)
    r = np.sqrt(X**2 + Y**2)

    q=np.empty([grid.meqn,len(x),len(y)], order = 'F')
    q[0,:,:] = (np.abs(r-0.5)<=width)*(1.+np.cos(np.pi*(r-0.5)/width))
    q[1,:,:] = 0.
    q[2,:,:] = 0.
    grid.q=q


def acoustics2D(iplot=False,petscPlot=False,useController=True):
    """
    Example python script for solving the 2d acoustics equations.
    """

    from petclaw.grid import Dimension
    from petclaw.grid import Grid
    from pyclaw.solution import Solution
    from petclaw.evolve.petclaw import PetClawSolver2D
    from pyclaw.controller import Controller
    from petclaw import plot

    # Initialize grid
    mx=100; my=100
    x = Dimension('x',-1.0,1.0,mx,mthbc_lower=1,mthbc_upper=1)
    y = Dimension('y',-1.0,1.0,my,mthbc_lower=1,mthbc_upper=1)
    grid = Grid([x,y])

    #Set global variables
    rho = 1.0
    bulk = 4.0
    cc = np.sqrt(bulk/rho)
    zz = rho*cc
    grid.aux_global['rho']= rho
    grid.aux_global['bulk']=bulk
    grid.aux_global['zz']= zz
    grid.aux_global['cc']=cc
    from dimsp2 import cparam
    for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)

    grid.meqn = 3
    grid.mbc = 2
    grid.t = 0.0
    tfinal = 0.27
    qinit(grid)
    inital_solution = Solution(grid)

    solver = PetClawSolver2D()
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.mwaves = 2
    solver.mthlim = [4]*solver.mwaves

    claw = Controller()
    claw.keep_copy = True
    # The output format MUST be set to petsc!
    claw.output_format = 'petsc'
    claw.tfinal = tfinal
    claw.solutions['n'] = inital_solution
    claw.solver = solver

    # Solve
    status = claw.run()

    if petscPlot:
        plot.plotPetsc(claw)

    if iplot:
        plot.plotInteractive()

    pressure=claw.frames[claw.nout].grid.gqVec.getArray().reshape([mx,my,grid.meqn])[:,:,0]
    return pressure


if __name__=="__main__":
    acoustics2D()
