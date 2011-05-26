#!/usr/bin/env python
# encoding: utf-8

"""
2D shallow water equations.
"""
#===========================================================================
# Import libraries
#===========================================================================

import numpy as np
from petclaw import plot
#import pdb


def qinit(grid,hl,ul,vl,hr,ur,vr,radDam):
    # Create an array with fortran native ordering
    xCenter = grid.x.center
    yCenter = grid.y.center
    Y,X = np.meshgrid(yCenter,xCenter)

    q=np.empty([grid.meqn,len(xCenter),len(yCenter)], order = 'F')
    for i,x in enumerate(xCenter):
        for j,y in enumerate(yCenter):
            #r = np.sqrt(xCenter[i]**2 + yCenter[j]**2)
            if np.sqrt(x**2+y**2)<=radDam: 
                q[0,i,j] = hl
                q[1,i,j] = hl*ul
                q[2,i,j] = hl*vl
            else: 
                q[0,i,j] = hr
                q[1,i,j] = hr*ur
                q[2,i,j] = hr*vr
    
    grid.q=q


    
def shallow2D(use_PETSc=False,kernel_language='Fortran',iplot=True,userController=True,petscPlot=False,htmlplot=False,outdir='./_output'):
    #===========================================================================
    # Import libraries
    #===========================================================================
    from petclaw.grid import Dimension
    from petclaw.grid import Grid
    from pyclaw.solution import Solution
    from petclaw.evolve.clawpack import PetClawSolver2D
    from pyclaw.controller import Controller
    from petclaw import plot

    #===========================================================================
    # Initialize grids, then initialize the solution associated to the grid and
    # finally initialize aux array
    #===========================================================================

    # Grid:
    xlower = -2.5
    xupper = 2.5
    mx = 50
    
    ylower = -2.5
    yupper = 2.5
    my = 50
    x = Dimension('x',xlower,xupper,mx,mthbc_lower=1,mthbc_upper=1)
    y = Dimension('y',ylower,yupper,my,mthbc_lower=1,mthbc_upper=1)
    grid = Grid([x,y])

    grid.meqn = 3  # Number of equations
    grid.mbc = 2   # Number of ghost cells
    grid.t = 0.    # Initial time


    # Parameters
    grav = 1.0
    grid.aux_global['grav'] = grav
    from classic2 import cparam
    for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)

    # Initial solution
    # ================
    # Initialize petsc Structures for q. 
    # This must be called before grid.x.center and such can be accessed.
    grid.init_q_petsc_structures()

    # Characteristcs of the dam break problem
    damRadius = 0.5
    hl = 2.
    ul = 0.
    vl = 0.
    hr = 1.
    ur = 0.
    vr = 0.
    
    qinit(grid,hl,ul,vl,hr,ur,vl,damRadius) # This function is defined above

    init_solution = Solution(grid)



    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    solver = PetClawSolver2D()
    solver.kernel_language = 'Fortran'
    solver.mwaves = 3
    if solver.kernel_language =='Python': solver.set_riemann_solver('shallow_roe')
    solver.mthlim = [4]*solver.mwaves



    #===========================================================================
    # Setup controller and controller paramters
    #===========================================================================
    claw = Controller()
    claw.keep_copy = True
    claw.tfinal = 5.0
    claw.solutions['n'] = init_solution
    claw.solver = solver
    claw.outdir = outdir



    #===========================================================================
    # Solve the problem
    #===========================================================================
    status = claw.run()

  

    #===========================================================================
    # Plot results
    #===========================================================================
    if iplot:     plot.plotInteractive(outdir=outdir,format=claw.output_format)
    if htmlplot:  plot.plotHTML(outdir=outdir,format=claw.output_format)




if __name__=="__main__":
    import sys
    if len(sys.argv)>1:
        from petclaw.util import _info_from_argv
        args, kwargs = _info_from_argv(sys.argv)
        error=shallow2D(*args,**kwargs)
        print 'Error: ',error
    else: shallow2D()






