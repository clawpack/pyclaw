#!/usr/bin/env python
# encoding: utf-8

import numpy as np

# material parameters
E1=1.;   p1=1.
E2=4.;   p2=4.
# interface parameters
alphax=0.5; deltax=1.
alphay=0.5; deltay=1.
# Linearity parameters
linearity_mat1=2; linearity_mat2=2
# heterogeneity type
het_type='checkerboard'
#het_type='sinusoidal'
#het_type='smooth_checkerboard'
sharpness=10

def qinit(state,A,x0,y0,varx,vary):
    r""" Set initial conditions for q."""
    x =state.grid.x.centers; y =state.grid.y.centers
    # Create meshgrid
    [yy,xx]=np.meshgrid(y,x)
    s=A*np.exp(-(xx-x0)**2/(2*varx)-(yy-y0)**2/(2*vary)) #sigma(@t=0)
    #parameters from aux
    linearity_mat=state.aux[2,:]
    E=state.aux[1,:]
    #initial condition
    state.q[0,:,:]=np.where(linearity_mat==1,1,0)*s/E+np.where(linearity_mat==2,1,0)*np.log(s+1)/E
    state.q[1,:,:]=0; state.q[2,:,:]=0

def setaux(x,y):
    r"""Creates a matrix representing every grid cell in the domain, 
    whose size is len(x),len(y)
    Each entry of the matrix contains a vector of size 3 with:
    The material density p
    The young modulus E
    A flag indicating which material the grid is made of.
    The domain pattern is a checkerboard."""

    aux = np.empty((4,len(x),len(y)), order='F')
    if het_type == 'checkerboard':
        # xfrac and yfrac are x and y relative to deltax and deltay resp.
        xfrac=x-np.floor(x/deltax)*deltax
        yfrac=y-np.floor(y/deltay)*deltay
        # create a meshgrid out of xfrac and yfrac
        [yyfrac,xxfrac]=np.meshgrid(yfrac,xfrac)
        # density 
        aux[0,:,:]=p1*(xxfrac<=alphax*deltax)*(yyfrac<=alphay*deltay)\
            +p1*(xxfrac >alphax*deltax)*(yyfrac >alphay*deltay)\
            +p2*(xxfrac >alphax*deltax)*(yyfrac<=alphay*deltay)\
            +p2*(xxfrac<=alphax*deltax)*(yyfrac >alphay*deltay)
        #Young modulus
        aux[1,:,:]=E1*(xxfrac<=alphax*deltax)*(yyfrac<=alphay*deltay)\
            +E1*(xxfrac >alphax*deltax)*(yyfrac >alphay*deltay)\
            +E2*(xxfrac >alphax*deltax)*(yyfrac<=alphay*deltay)\
            +E2*(xxfrac<=alphax*deltax)*(yyfrac >alphay*deltay)
        # linearity of material
        aux[2,:,:]=linearity_mat1*(xxfrac<=alphax*deltax)*(yyfrac<=alphay*deltay)\
            +linearity_mat1*(xxfrac >alphax*deltax)*(yyfrac >alphay*deltay)\
            +linearity_mat2*(xxfrac >alphax*deltax)*(yyfrac<=alphay*deltay)\
            +linearity_mat2*(xxfrac<=alphax*deltax)*(yyfrac >alphay*deltay)
    elif het_type == 'sinusoidal' or het_type == 'smooth_checkerboard':
        [yy,xx]=np.meshgrid(y,x)
        Amp_p=np.abs(p1-p2)/2; offset_p=(p1+p2)/2
        Amp_E=np.abs(E1-E2)/2; offset_E=(E1+E2)/2
        if het_type == 'sinusoidal':
            frec_x=2*np.pi/deltax; frec_y=2*np.pi/deltay
            fun=np.sin(frec_x*xx)*np.sin(frec_y*yy)
        else:
            fun_x=xx*0; fun_y=yy*0
            for i in xrange(0,1+int(np.ceil((x[-1]-x[0])/(deltax*0.5)))):
                fun_x=fun_x+(-1)**i*np.tanh(sharpness*(xx-deltax*i*0.5))
            for i in xrange(0,1+int(np.ceil((y[-1]-y[0])/(deltay*0.5)))):
                fun_y=fun_y+(-1)**i*np.tanh(sharpness*(yy-deltay*i*0.5))
            fun=fun_x*fun_y
        aux[0,:,:]=Amp_p*fun+offset_p
        aux[1,:,:]=Amp_E*fun+offset_E
        aux[2,:,:]=linearity_mat1
    return aux

def b4step(solver,state):
    r"""put in aux[3,:,:] the value of q[0,:,:] (eps). This is required in rptpv.f"""
    state.aux[3,:,:] = state.q[0,:,:]

    # To set to 0 1st 1/2 of the domain. Used in rect domains with PBC in x
    if state.problem_data['turnZero_half_2D']==1:
        if state.t>=state.problem_data['t_turnZero'] and state.t<=state.problem_data['t_turnZero']+1:
            Y,X = np.meshgrid(state.grid.y,state.grid.x)
            state.q = state.q * (X<state.grid.upper[0]/2)

def compute_p(state):
    state.p[0,:,:]=np.exp(state.q[0,:,:]*state.aux[1,:,:])-1

def compute_F(state):
    rho = state.aux[0,:,:]; E = state.aux[1,:,:]
    
    #Compute the entropy
    u = state.q[1,:,:]/rho
    v = state.q[2,:,:]/rho

    nrg=rho * (u**2 + v**2)/2.

    eps = state.q[0,:,:]
    sigma = np.exp(E*eps) - 1.
    sigint = (sigma-np.log(sigma+1.))/E

    dx=state.grid.delta[0]; dy=state.grid.delta[1]
    
    state.F[0,:,:] = (sigint+nrg)*dx*dy 
    state.F[1,:,:] = 10*state.F[0,:,:]
    state.F[2,:,:] = 100*state.F[0,:,:]

def gauge_pfunction(q,aux):
    p = np.exp(q[0]*aux[1])-1
    return [p,10*p]

def psystem2D(iplot=False,kernel_language='Fortran',htmlplot=False,
              use_petsc=False,outdir='./_output',solver_type='classic',
              disable_output=False):

    """
    Solve the p-system in 2D with variable coefficients
    """
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    ####################################
    ######### MAIN PARAMETERS ##########
    ####################################
    # Domain
    x_lower=0.25; x_upper=20.25
    y_lower=0.25; y_upper=20.25
    # cells per layer
    Ng=10
    mx=(x_upper-x_lower)*Ng; my=(y_upper-y_lower)*Ng
    # Initial condition parameters
    A=10.
    x0=0.25 # Center of initial perturbation
    y0=0.25 # Center of initial perturbation
    varx=0.5; vary=0.5 # Width of initial perturbation

    # Boundary conditions
    bc_x_lower=pyclaw.BC.wall; bc_x_upper=pyclaw.BC.extrap
    bc_y_lower=pyclaw.BC.wall; bc_y_upper=pyclaw.BC.extrap

    # Turning off 1st half of the domain. Useful in rect domains
    turnZero_half_2D=0 #flag
    t_turnZero=50
    num_output_times=10
    # restart options
    restart_from_frame = None

    if solver_type=='classic':
        solver = pyclaw.ClawSolver2D()
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver2D()

    if kernel_language != 'Fortran':
        raise Exception('Unrecognized value of kernel_language for 2D psystem')

    from clawpack import riemann
    solver.rp = riemann.rp2_psystem

    solver.num_waves = 2
    solver.limiters = pyclaw.limiters.tvd.superbee

    solver.bc_lower[0]=bc_x_lower
    solver.bc_upper[0]=bc_x_upper
    solver.bc_lower[1]=bc_y_lower
    solver.bc_upper[1]=bc_y_upper
    solver.aux_bc_lower[0]=bc_x_lower
    solver.aux_bc_upper[0]=bc_x_upper
    solver.aux_bc_lower[1]=bc_y_lower
    solver.aux_bc_upper[1]=bc_y_upper

    solver.fwave = True
    solver.cfl_max = 0.9
    solver.cfl_desired = 0.8
    solver.before_step = b4step
    solver.dimensional_split=False

    #controller
    claw = pyclaw.Controller()
    claw.tfinal = 40.0
    claw.solver = solver
    claw.outdir = outdir
    
    if restart_from_frame is not None:
        claw.solution = pyclaw.Solution(restart_from_frame, format='petsc',read_aux=False)
        claw.solution.state.mp = 1
        grid = claw.solution.domain.grid
        claw.solution.state.aux = setaux(grid.x.centers,grid.y.centers)
        claw.num_output_times = num_output_times - restart_from_frame
        claw.start_frame = restart_from_frame
    else:
        ####################################
        ####################################
        ####################################
        #Creation of Domain
        x = pyclaw.Dimension('x',x_lower,x_upper,mx)
        y = pyclaw.Dimension('y',y_lower,y_upper,my)
        domain = pyclaw.Domain([x,y])
        num_eqn = 3
        num_aux = 4
        state = pyclaw.State(domain,num_eqn,num_aux)
        state.mF = 3
        #Set global parameters
        state.problem_data = {}
        state.problem_data['turnZero_half_2D'] = turnZero_half_2D
        state.problem_data['t_turnZero'] = t_turnZero
        state.mp = 1

        grid = state.grid
        state.aux = setaux(grid.x.centers,grid.y.centers)
        #Initial condition
        qinit(state,A,x0,y0,varx,vary)

        claw.solution = pyclaw.Solution(state,domain)
        claw.num_output_times = num_output_times

    #claw.p_function = p_function
    if disable_output:
        claw.output_format = None
    claw.compute_F = compute_F
    state.keep_gauges = True
    grid.add_gauges([[0.25,0.25],[17.85,1.25],[3.25,18.75],[11.75,11.75]])
    solver.compute_gauge_values = gauge_pfunction
    claw.write_aux_init = False

    #Solve
    status = claw.run()

    if iplot:    pyclaw.plot.interactive_plot()
    if htmlplot: pyclaw.plot.html_plot()

    return claw.solution.state


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(psystem2D)
