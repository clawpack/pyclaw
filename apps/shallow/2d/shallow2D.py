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

    grid.zeros_q()
    for i,x in enumerate(xCenter):
        for j,y in enumerate(yCenter):
            #r = np.sqrt(xCenter[i]**2 + yCenter[j]**2)
            if np.sqrt(x**2+y**2)<=radDam: 
                grid.q[0,i,j] = hl
                grid.q[1,i,j] = hl*ul
                grid.q[2,i,j] = hl*vl
            else: 
                grid.q[0,i,j] = hr
                grid.q[1,i,j] = hr*ur
                grid.q[2,i,j] = hr*vr

    
def shallow2D(use_petsc=False,kernel_language='Fortran',iplot=True,htmlplot=False,outdir='./_output',solver_type='classic'):
    #===========================================================================
    # Import libraries
    #===========================================================================
    import numpy as np

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D()
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D()

    solver.mwaves = 3
    if solver.kernel_language =='Python': solver.set_riemann_solver('shallow_roe')
    solver.mthlim = pyclaw.limiters.MC

    solver.mthbc_lower[0] = pyclaw.BC.outflow
    solver.mthbc_upper[0] = pyclaw.BC.outflow
    solver.mthbc_lower[1] = pyclaw.BC.outflow
    solver.mthbc_upper[1] = pyclaw.BC.outflow

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
    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    grid = pyclaw.Grid([x,y])

    grid.meqn = 3  # Number of equations
    grid.mbc = solver.mbc   # Number of ghost cells

    # Parameters
    grav = 1.0
    grid.aux_global['grav'] = grav

    # Initial solution
    # ================
    # Riemann states of the dam break problem
    damRadius = 0.5
    hl = 2.
    ul = 0.
    vl = 0.
    hr = 1.
    ur = 0.
    vr = 0.
    
    qinit(grid,hl,ul,vl,hr,ur,vl,damRadius) # This function is defined above

    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = 5.0
    claw.solution = pyclaw.Solution(grid)
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






