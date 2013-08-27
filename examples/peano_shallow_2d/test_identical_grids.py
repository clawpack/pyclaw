from nose.plugins.attrib import attr

if __name__=="__main__":
    import nose
    nose.main()
    
    
def qinit(state,hl,ul,vl,hr,ur,vr,radDam):
    import numpy as np
    x0=0.5
    y0=0.5
    xCenter = state.grid.x.centers
    yCenter = state.grid.y.centers
    Y,X = np.meshgrid(yCenter,xCenter)
    r = np.sqrt((X-x0)**2 + (Y-y0)**2)
    state.q[0,:,:] = hl*(r<=radDam) + hr*(r>radDam)
    state.q[1,:,:] = hl*ul*(r<=radDam) + hr*ur*(r>radDam)
    state.q[2,:,:] = hl*vl*(r<=radDam) + hr*vr*(r>radDam)
    
def setup_solver():
    import pyclaw
    import peanoclaw
    solver = pyclaw.ClawSolver2D()
    solver.limiters = pyclaw.limiters.tvd.MC
    solver.dimensional_split=1

    import riemann
    solver.rp = riemann.rp2_shallow_roe_with_efix
    solver.num_waves = 3

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_lower[1] = pyclaw.BC.wall
    solver.bc_upper[1] = pyclaw.BC.wall
    
    return solver

@attr(petsc=False)
@attr(peanoclaw=True)
def test_3x3_grid():
    r"""This test simply solves a 3x3 grid, once with PyClaw and once as one patch 
    with PeanoClaw. In the end it checks if the resulting qbcs match.
    """
    #===========================================================================
    # Import libraries
    #===========================================================================
    import numpy as np
    import pyclaw
    import peanoclaw

    #===========================================================================
    # Setup solver and solver parameters for PyClaw run
    #===========================================================================
    pyclaw_solver = setup_solver()
    
    #===========================================================================
    # Setup solver and solver parameters for PeanoClaw run
    #===========================================================================
    peanoclaw_solver = setup_solver()
    peano_solver = peanoclaw.Solver(peanoclaw_solver, 1.0)

    #===========================================================================
    # Initialize domain and state, then initialize the solution associated to the 
    # state and finally initialize aux array
    #===========================================================================

    # Domain:
    xlower = 0.0
    xupper = 1.0
    mx = 3
    ylower = 0.0
    yupper = 1.0
    my = 3
    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    domain = pyclaw.Domain([x,y])

    num_eqn = 3  # Number of equations
    pyclaw_state = pyclaw.State(domain,num_eqn)
    peanoclaw_state = pyclaw.State(domain, num_eqn)

    grav = 1.0 # Parameter (global auxiliary variable)
    pyclaw_state.problem_data['grav'] = grav
    peanoclaw_state.problem_data['grav'] = grav

    # Initial solution
    # ================
    # Riemann states of the dam break problem
    damRadius = 0.2
    hl = 2.
    ul = 0.
    vl = 0.
    hr = 1.
    ur = 0.
    vr = 0.
    
    qinit(pyclaw_state,hl,ul,vl,hr,ur,vl,damRadius) # This function is defined above
    qinit(peanoclaw_state,hl,ul,vl,hr,ur,vl,damRadius) # This function is defined above

    tfinal = 1.0
    #===========================================================================
    # Set up controller and controller parameters for PyClaw run
    #===========================================================================
    pyclaw_controller = pyclaw.Controller()
    pyclaw_controller.tfinal = tfinal
    pyclaw_controller.solution = pyclaw.Solution(pyclaw_state,domain)
    pyclaw_controller.solver = pyclaw_solver
    pyclaw_controller.num_output_times = 1
    
    pyclaw_controller.run()
    
    #===========================================================================
    # Set up controller and controller parameters for PyClaw run
    #===========================================================================
    peanoclaw_controller = pyclaw.Controller()
    peanoclaw_controller.tfinal = tfinal
    peanoclaw_controller.solution = pyclaw.Solution(peanoclaw_state,domain)
    peanoclaw_controller.solver = peano_solver
    peanoclaw_controller.num_output_times = 1
    
    peanoclaw_controller.run()
    
    assert(np.max(np.abs(pyclaw_solver.qbc - peanoclaw_solver.qbc)) < 1e-9)
    
    
