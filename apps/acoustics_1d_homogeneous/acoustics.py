#!/usr/bin/env python
# encoding: utf-8
    
def acoustics(use_petsc=False,kernel_language='Fortran',solver_type='classic',iplot=False,htmlplot=False,outdir='./_output',weno_order=5):
    """
    This example solves the 1-dimensional acoustics equations in a homogeneous
    medium.
    """
    from numpy import sqrt, exp, cos
    from riemann import rp1_acoustics

    #=================================================================
    # Import the appropriate classes, depending on the options passed
    #=================================================================
    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver1D()
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
        solver.weno_order=weno_order
    else: raise Exception('Unrecognized value of solver_type.')

    #========================================================================
    # Instantiate the solver and define the system of equations to be solved
    #========================================================================
    solver.kernel_language=kernel_language
    from riemann import rp_acoustics
    solver.num_waves=rp_acoustics.num_waves

    if kernel_language=='Python': 
        solver.rp = rp_acoustics.rp_acoustics_1d
    else:
        from riemann import rp1_acoustics
        solver.rp = rp1_acoustics

    solver.limiters = pyclaw.limiters.tvd.MC

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic

    #========================================================================
    # Instantiate the domain and set the boundary conditions
    #========================================================================
    x = pyclaw.Dimension('x',0.0,1.0,100)
    domain = pyclaw.Domain(x)
    num_eqn = 2
    state = pyclaw.State(domain,num_eqn)

    #========================================================================
    # Set problem-specific variables
    #========================================================================
    rho = 1.0
    bulk = 1.0

    state.problem_data['rho']=rho
    state.problem_data['bulk']=bulk
    state.problem_data['zz']=sqrt(rho*bulk)
    state.problem_data['cc']=sqrt(bulk/rho)
 

    #========================================================================
    # Set the initial condition
    #========================================================================
    xc=domain.grid.x.centers
    beta=100; gamma=0; x0=0.75
    state.q[0,:] = exp(-beta * (xc-x0)**2) * cos(gamma * (xc - x0))
    state.q[1,:] = 0.

    solver.dt_initial=domain.grid.delta[0]/state.problem_data['cc']*0.1

    #========================================================================
    # Set up the controller object
    #========================================================================
    claw = pyclaw.Controller()
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.keep_copy = True
    claw.num_output_times = 5

    claw.tfinal = 1.0

    # Solve
    status = claw.run()

    # Plot results
    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

    return claw

if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(acoustics)

def test_1d_acoustics():

    from pyclaw.util import test_app_variants
    from pyclaw.util import check_diff
    import numpy as np

    def verify_classic(claw):
        q0=claw.frames[0].state.q.reshape([-1])
        qfinal=claw.frames[claw.num_output_times].state.q.reshape([-1])
        dx=claw.solution.domain.grid.delta[0]
        test = dx*np.sum(np.abs(qfinal-q0))
        expected = 0.00104856594174
        return check_diff(test, expected, abstol=1e-5)

    def verify_sharpclaw(claw):
        q0=claw.frames[0].state.q.reshape([-1])
        qfinal=claw.frames[claw.num_output_times].state.q.reshape([-1])
        dx=claw.solution.domain.grid.delta[0]
        test = dx*np.sum(np.abs(qfinal-q0))
        expected = 0.000298879563857
        return check_diff(test, expected, abstol=1e-5)

    test_app_variants(acoustics, verify_classic, python_kernel=True, solver_type='classic')
    test_app_variants(acoustics, verify_sharpclaw, python_kernel=True, solver_type='sharpclaw')
