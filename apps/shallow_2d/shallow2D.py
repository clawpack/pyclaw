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
#import pdb  # Debugger


def qinit(state,hl,ul,vl,hr,ur,vr,radDam):
    x0=0.
    y0=0.
    xCenter = state.grid.x.center
    yCenter = state.grid.y.center
    Y,X = np.meshgrid(yCenter,xCenter)
    r = np.sqrt((X-x0)**2 + (Y-y0)**2)
    state.q[0,:,:] = hl*(r<=radDam) + hr*(r>radDam)
    state.q[1,:,:] = hl*ul*(r<=radDam) + hr*ur*(r>radDam)
    state.q[2,:,:] = hl*vl*(r<=radDam) + hr*vr*(r>radDam)

    
def shallow2D(use_petsc=False,iplot=0,htmlplot=False,outdir='./_output',solver_type='classic'):
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

    solver.num_waves = 3
    solver.limiters = pyclaw.limiters.tvd.MC

    solver.bc_lower[0] = pyclaw.BC.outflow
    solver.bc_upper[0] = pyclaw.BC.reflecting
    solver.bc_lower[1] = pyclaw.BC.outflow
    solver.bc_upper[1] = pyclaw.BC.reflecting
    solver.dim_split=1

    #===========================================================================
    # Initialize grid and state, then initialize the solution associated to the 
    # state and finally initialize aux array
    #===========================================================================

    # Grid:
    xlower = -2.5
    xupper = 2.5
    mx = 150
    ylower = -2.5
    yupper = 2.5
    my = 150
    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    grid = pyclaw.Grid([x,y])

    num_eqn = 3  # Number of equations
    state = pyclaw.State(grid,num_eqn)

    grav = 1.0 # Parameter (global auxiliary variable)
    state.aux_global['grav'] = grav

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
    
    qinit(state,hl,ul,vl,hr,ur,vl,damRadius) # This function is defined above

    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = 2.5
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir = outdir
    claw.nout = 10

    #===========================================================================
    # Solve the problem
    #===========================================================================
    status = claw.run()

    #===========================================================================
    # Plot results
    #===========================================================================
    if iplot:     plot.interactive_plot(outdir=outdir,format=claw.output_format)
    if htmlplot:  plot.html_plot(outdir=outdir,format=claw.output_format)


if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(shallow2D)
    print 'Error: ', output





