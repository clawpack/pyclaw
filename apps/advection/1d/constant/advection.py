#!/usr/bin/env python
# encoding: utf-8

"""
1D shallow advection equation with constant convective velocity.
"""


def advection(kernel_language='Fortran',iplot=False,htmlplot=False,use_petsc=False,solver_type='sharpclaw',outdir='./_output'):
    #===========================================================================
    # Import libraries
    #===========================================================================
    import numpy as np

    if use_petsc:
        import petclaw as pyclaw
    else:
        import pyclaw

    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D()
    else:
        solver = pyclaw.ClawSolver1D()

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    solver.kernel_language = kernel_language

    from riemann import rp_advection
    solver.mwaves = rp_advection.mwaves
    if solver.kernel_language=='Python': 
        solver.rp = rp_advection.rp_advection_1d

    solver.bc_lower[0] = 2
    solver.bc_upper[0] = 2

    # Time integrator
    solver.time_integrator='SSP104'
    solver.dt_initial = 4.0/1000.
    solver.cfl_max = 4.1
    solver.cfl_desired = 4.0
    solver.dt_variable=False
    #print solver.dt_initial

    mm=1000
    #===========================================================================
    # Initialize grids and then initialize the solution associated to the grid
    #===========================================================================
    #x = pyclaw.Dimension('x',0.0,1.0,100)
    x=pyclaw.Dimension('x',0.0,1.0,mm)
    grid = pyclaw.Grid(x)
    meqn = 1
    state = pyclaw.State(grid,meqn)
    state.aux_global['u']=1.

    xc=grid.x.center
    #beta=100; gamma=0; x0=0.5
    #state.q[0,:] = np.exp(-beta * (xc-x0)**2) * np.cos(gamma * (xc - x0))
    xx=np.linspace(0.0,1.0,mm+1)
    #for i in range(mm):
    #    a=xx[i]
    #    b=xx[i+1]
    #    state.q[0,i]=mm*(np.sin(2*(xx[i+1]-0.5)*np.pi)-np.sin(2*(xx[i]-0.5)*np.pi))/(2*np.pi)
    #np.savetxt("test_sol.txt",state.q[0,:])
    state.q[0,:]=mm*(np.sin(2*(xx[1:]-0.5)*np.pi)-np.sin(2*(xx[0:-1]-0.5)*np.pi))/(2*np.pi)
   # np.savetxt("ext_1000.txt",state.q[0,:])
    #state.q[0,:]= np.cos(2*(xc-0.5)*np.pi);


    #===========================================================================
    # Setup controller and controller paramters
    #===========================================================================
    #claw = pyclaw.Controller()
    #claw.solution = pyclaw.Solution(state)
    #claw.solver = solver
    #claw.outdir = outdir
    #claw.tfinal = 1.0

    claw = pyclaw.Controller()
    claw.keep_copy = True
    claw.outstyle = 1
    claw.nout = 10
    claw.tfinal = 3
    claw.solution = pyclaw.Solution(state)
    claw.solver = solver
    claw.outdir = outdir
    
   
#===========================================================================
    # Solve the problem
    #===========================================================================
    status = claw.run()
    sol = claw.frames[claw.nout].state.q[0,:]
    np.savetxt("M104_4.txt", sol)
    
    #===========================================================================
    # Plot results
    #===========================================================================
    if htmlplot:  pyclaw.plot.html_plot(outdir=outdir)
    if iplot:     pyclaw.plot.interactive_plot(outdir=outdir)

if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(advection)
