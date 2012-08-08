#!/usr/bin/env python
# encoding: utf-8

from __future__ import division

"""
2D shallow water equations.
"""
#===========================================================================
# Import libraries
#===========================================================================

import numpy as np

def qinit(state,hl,ul,vl,hr,ur,vr,radDam):
    x0=0.5
    y0=0.5
    xCenter = state.grid.x.centers
    yCenter = state.grid.y.centers
    
    Y,X = np.meshgrid(yCenter,xCenter)
    r = np.sqrt((X-x0)**2 + (Y-y0)**2)
    state.q[0,:,:] = hl*(r<=radDam) + hr*(r>radDam)
    state.q[1,:,:] = hl*ul*(r<=radDam) + hr*ur*(r>radDam)
    state.q[2,:,:] = hl*vl*(r<=radDam) + hr*vr*(r>radDam)
    
def refinement_criterion(state):
    if((state.q[0,:,:].max() - state.q[0,:,:].min()) > 0.2):
    	return 1.0/36.0
    else:
    	return 1.0/18.0

    
def shallow2D(use_petsc=False,iplot=0,htmlplot=False,outdir='./_output',solver_type='classic',amr_type=None):
    #===========================================================================
    # Import libraries
    #===========================================================================
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D()
        solver.limiters = pyclaw.limiters.tvd.MC
        solver.dimensional_split=1
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D()

    from clawpack import riemann
    solver.rp = riemann.rp2_shallow_roe_with_efix
    solver.num_waves = 3

    solver.bc_lower[0] = pyclaw.BC.wall
    solver.bc_upper[0] = pyclaw.BC.wall
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.wall

    #===========================================================================
    # Initialize domain and state, then initialize the solution associated to the 
    # state and finally initialize aux array
    #===========================================================================

    # resolution of each grid
    mgrid = 40

    # number of initial AMR grids in each dimension
    msubgrid = 3

    if amr_type is not None:
        m = mgrid
    else:
        # number of Domain grid cells expressed as the product of
        # grid resolution and the number of AMR sub-grids for
        # easy comparison between the two methods
        m = mgrid*msubgrid
    
    mx = m
    my = m

    # Domain:
    xlower = 0
    xupper = 1
    ylower = 0
    yupper = 1

    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    domain = pyclaw.Domain([x,y])

    num_eqn = 3  # Number of equations
    state = pyclaw.State(domain,num_eqn)

    grav = 1.0 # Parameter (global auxiliary variable)
    state.problem_data['grav'] = grav

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

    # this closure is used by AMR-style codes
    def qinit_callback(state):
        qinit(state,hl,ul,vl,hr,ur,vl,damRadius)

    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = 2.5

    if amr_type is not None:        
        if amr_type == 'peano':
            import clawpack.peanoclaw as amrclaw
            claw.solver = amrclaw.Solver(solver
                                        ,1/(mgrid*msubgrid)
                                        ,qinit_callback
                                        #,refinement_criterion=refinement_criterion
                                        )
            claw.solution = amrclaw.Solution(state, domain)
        else:
            raise Exception('unsupported amr_type %s' % amr_type)
    else:
        claw.solver = solver
        claw.solution = pyclaw.Solution(state,domain)
        
    claw.keep_copy = True
    claw.outdir = outdir
    claw.num_output_times = 10

    #===========================================================================
    # Solve the problem
    #===========================================================================
    status = claw.run()

    #===========================================================================
    # Plot results
    #===========================================================================
    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(shallow2D)





