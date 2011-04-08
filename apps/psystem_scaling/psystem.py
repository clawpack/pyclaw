#!/usr/bin/env python
# encoding: utf-8

import numpy as np

def qinit(grid,A,x0,y0,varx,vary):
    # Set initial conditions for q.
    # Create an array with fortran native ordering
    x =grid.x.center
    y =grid.y.center
    q=np.empty([grid.meqn,len(x),len(y)], order = 'F')
    
    # Create meshgrid
    # start form left top to right bottom
    [yy,xx]=np.meshgrid(y,x)
    s=A*np.exp(-(xx-x0)**2/(2*varx)-(yy-y0)**2/(2*vary)) #sigma(@t=0)

    #Material parameter
    E=1

    # initial conditions
    q[0,:,:]=s/E
    q[1,:,:]=0
    q[2,:,:]=0

    grid.q=q

def psystem2D(iplot=False,petscPlot=False,useController=True):
    """
    psystem in 2D with variable coefficients
    """
    import time
    time_import_start=time.time()
    from petclaw.grid import Dimension
    from petclaw.grid import Grid
    from pyclaw.solution import Solution
    from petclaw.evolve.solver import PetClawSolver2D
    from pyclaw.controller import Controller
    from petclaw import plot
    time_import_end=time.time()
    print 'Import time '+str(time_import_end-time_import_start)+' sec'

    time_run_start=time.time()
####################################
######### MAIN PARAMETERS ##########
####################################
# material parameters
    #hardcoded in the Riemman Solver
    #p1=p2=E1=E2=1
# Linearity parameters
    #hardcoded in the Riemman Solver
    #linearity_mat1=linearity_mat2=1
# Domain
    x_lower=0.0; x_upper=10.0
    y_lower=0.0; y_upper=10.0
# Grid cells
    mx=100; my=100
# Initial condition parameters
    A=5.
    x0=x_upper/2.; y0=y_upper/2.
    varx=1.0; vary=1.0
# Boundary conditions
    mthbc_x_lower=1; mthbc_x_upper=1
    mthbc_y_lower=1; mthbc_y_upper=1
# Regarding time
    nout=1
    tfinal=0.45*(x_upper-x_lower)/mx*1000
# Ghost cells
    mbc=2
# number of equations
    meqn=3
####################################
####################################
####################################
# creation of grid
    x = Dimension('x',x_lower,x_upper,mx,mthbc_lower=mthbc_x_lower,mthbc_upper=mthbc_x_upper,mbc=mbc)
    y = Dimension('y',y_lower,y_upper,my,mthbc_lower=mthbc_y_lower,mthbc_upper=mthbc_y_upper,mbc=mbc)
    grid = Grid([x,y])
    grid.meqn=meqn
    grid.mbc=mbc
    grid.t=0.0
    grid.init_q_petsc_structures() #Initialize Petsc structures

# initial conditions
    qinit(grid,A,x0,y0,varx,vary)

    inital_solution = Solution(grid)

    solver = PetClawSolver2D()
    solver.max_steps=1000000
    solver.dt_variable=True
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.order = 2
    solver.mwaves = 2
    solver.fwave = True 
    solver.mthlim = [4]*solver.mwaves

    claw = Controller()
    claw.keep_copy = True
    claw.nout = nout
    claw.outdir = './_output/'
    # The output format MUST be set to petsc!
    claw.output_format = 'petsc'
    claw.tfinal = tfinal
    claw.solutions['n'] = inital_solution
    claw.solver = solver

    # Solve
    status = claw.run()
         
    time_run_end=time.time()
    print 'Running time '+str(time_run_end-time_run_start)
    
    #eps=claw.frames[claw.nout].grid.gqVec.getArray().reshape([mx,my,grid.meqn])[:,:,0]
    #return eps

if __name__=="__main__":
    psystem2D()
