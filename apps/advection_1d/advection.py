#!/usr/bin/env python
# encoding: utf-8
def advection(kernel_language='Python',iplot=False,htmlplot=False,
              use_petsc=False,solver_type='classic', weno_order=5,
              outdir='./_output'):
    """
    Example python script for solving the 1d advection equation.
    """
    import numpy as np

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
        solver.weno_order=weno_order
    else:
        solver = pyclaw.ClawSolver1D()

    solver.kernel_language = kernel_language
    from clawpack.riemann import rp_advection
    solver.num_waves = rp_advection.num_waves
    if solver.kernel_language=='Python':
        solver.rp = rp_advection.rp_advection_1d
    else:
        from clawpack import riemann
        solver.rp = riemann.rp1_advection

    solver.bc_lower[0] = 2
    solver.bc_upper[0] = 2

    x = pyclaw.Dimension('x',0.0,1.0,100)
    domain = pyclaw.Domain(x)
    num_eqn = 1
    state = pyclaw.State(domain,num_eqn)
    state.problem_data['u']=1.

    grid = state.grid
    xc=grid.x.centers
    beta=100; gamma=0; x0=0.75
    state.q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver

    if outdir is not None:
        claw.outdir = outdir
    else:
        claw.output_format = None

    claw.tfinal =1.0
    status = claw.run()

    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

    return claw

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(advection)
