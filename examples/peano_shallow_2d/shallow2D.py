#!/usr/bin/env python
# encoding: utf-8

"""
2D shallow water equations.
"""
#===========================================================================
# Import libraries
#===========================================================================

import numpy as np
#from petclaw import plot
#import pdb  # Debugger

def init(state):
    # Initial solution
    # ================
    # Riemann states of the dam break problem
    radDam = 0.2
    hl = 2.
    ul = 0.
    vl = 0.
    hr = 1.
    ur = 0.
    vr = 0.
    
    x0=0.5
    y0=0.5
    xCenter = state.grid.x.centers
    yCenter = state.grid.y.centers
    
    Y,X = np.meshgrid(yCenter,xCenter)
    r = np.sqrt((X-x0)**2 + (Y-y0)**2)
    state.q[0,:,:] = hl*(r<=radDam) + hr*(r>radDam)
    state.q[1,:,:] = hl*ul*(r<=radDam) + hr*ur*(r>radDam)
    state.q[2,:,:] = hl*vl*(r<=radDam) + hr*vr*(r>radDam)


    
def shallow2D(use_petsc=False,outdir='./_output',solver_type='classic', disable_output=False):
    #===========================================================================
    # Import libraries
    #===========================================================================
    import numpy as np
    import clawpack.peanoclaw as peanoclaw
    import clawpack.riemann as riemann

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        import clawpack.pyclaw as pyclaw

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    subdivisionFactor = 6
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D(riemann.shallow_roe_with_efix_2D)
        solver.limiters = pyclaw.limiters.tvd.MC
        solver.dimensional_split=1
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(riemann.shallow_roe_with_efix_2D)
    peanoSolver = peanoclaw.Solver(solver, (1./3.)/subdivisionFactor, init)
    
    solver.dt_initial = 1.0

    solver.bc_lower[0] = pyclaw.BC.wall
    solver.bc_upper[0] = pyclaw.BC.wall
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.wall
    
    #===========================================================================
    # Initialize domain and state, then initialize the solution associated to the 
    # state and finally initialize aux array
    #===========================================================================

    # Domain:
    from clawpack.pyclaw import geometry
    print((geometry.__file__))
    xlower = 0.0
    xupper = 1.0
    mx = subdivisionFactor
    ylower = 0.0
    yupper = 1.0
    my = subdivisionFactor
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = geometry.Domain([x,y])

    num_eqn = 3  # Number of equations
    state = pyclaw.State(domain,num_eqn)

    grav = 1.0 # Parameter (global auxiliary variable)
    state.problem_data['grav'] = grav
    
    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = 0.1
    claw.solution = peanoclaw.solution.Solution(state,domain)
    claw.solver = peanoSolver
    claw.outdir = outdir
    if disable_output:
        claw.output_format = None
    claw.num_output_times = 5

    return claw


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(shallow2D)
    print('Error: ', output)





