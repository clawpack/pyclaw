#!/usr/bin/env python
# encoding: utf-8

def setup(state_backend='pyclaw',outdir='./_output',solver_type='classic'):
    """
    Example python script for solving 1d traffic model:

    $$ q_t + umax( q(1-q) )_x = 0.$$
    """

    import numpy as np
    from clawpack import riemann
    from clawpack.pyclaw.util import get_state_backend
    pyclaw = get_state_backend(state_backend)

    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    if solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver1D(riemann.traffic_1D)
    else:
        solver = pyclaw.ClawSolver1D(riemann.traffic_1D)

    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap

    #===========================================================================
    # Initialize domain and then initialize the solution associated to the domain
    #===========================================================================
    x = pyclaw.Dimension('x',-1.0,1.0,500)
    domain = pyclaw.Domain(x)
    num_eqn = 1
    state = pyclaw.State(domain,num_eqn)

    grid = state.grid
    xc=grid.x.centers

    state.q[0,:] = 0.75*(xc<0) + 0.1*(xc>0.) 

    state.problem_data['efix']=True
    state.problem_data['umax']=1.

    #===========================================================================
    # Setup controller and controller parameters. Then solve the problem
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal =2.0
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot
    claw.keep_copy = True

    return claw

#--------------------------
def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.ylimits = [-0.1, 1.1]
    plotaxes.title = 'q[0]'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    
    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
