#!/usr/bin/python
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

    #parameters from aux
    linearity_mat=grid.aux[2,grid.mbc:-grid.mbc,grid.mbc:-grid.mbc]
    E=grid.aux[1,grid.mbc:-grid.mbc,grid.mbc:-grid.mbc]

    # initial conditions
    q[0,:,:]=np.where(linearity_mat==1,1,0)*s/E+np.where(linearity_mat==2,1,0)*np.log(s+1)/E
    q[1,:,:]=0
    q[2,:,:]=0

    grid.q=q

def setaux(grid, x,y,mthbcx,mthbcy,mbc=2,E1=1.,p1=1.,E2=1.,p2=1.,linearity_mat1=1,linearity_mat2=1,alphax=0.5,deltax=1.,alphay=0.5,deltay=1.):
# Creates a matrix representing every grid cell in the domain, 
#  whose size is len(x),len(y)
# Each entry of the matrix contains a vector of size 3 with:
#     The material density p
#     The young modulus E
#     A flag indicating which material the grid is made of
# The domain pattern is a checkerboard
    maux=3
    aux=np.empty([maux,len(x),len(y)],order='F')
    
    # xfrac and yfrac are x and y relative to deltax and deltay resp.
    xfrac=x-np.floor(x/deltax)*deltax
    yfrac=y-np.floor(y/deltay)*deltay
    
    # create a meshgrid out of xfrac and yfrac
    # start form left top to right bottom
    [yyfrac,xxfrac]=np.meshgrid(yfrac,xfrac)

    # density 
    aux[0,:,:]=p1*(xxfrac<=alphax*deltax)*(yyfrac<=alphay*deltay)+p1*(xxfrac>alphax*deltax)*(yyfrac>alphay*deltay)+p2*(xxfrac>alphax*deltax)*(yyfrac<=alphay*deltay)+p2*(xxfrac<=alphax*deltax)*(yyfrac>alphay*deltay)
    #Young modulus
    aux[1,:,:]=E1*(xxfrac<=alphax*deltax)*(yyfrac<=alphay*deltay)+E1*(xxfrac>alphax*deltax)*(yyfrac>alphay*deltay)+E2*(xxfrac>alphax*deltax)*(yyfrac<=alphay*deltay)+E2*(xxfrac<=alphax*deltax)*(yyfrac>alphay*deltay)
    # linearity of material
    aux[2,:,:]=linearity_mat1*(xxfrac<=alphax*deltax)*(yyfrac<=alphay*deltay)+linearity_mat1*(xxfrac>alphax*deltax)*(yyfrac>alphay*deltay)+linearity_mat2*(xxfrac>alphax*deltax)*(yyfrac<=alphay*deltay)+linearity_mat2*(xxfrac<=alphax*deltax)*(yyfrac>alphay*deltay)

#ghost cells    
    mx=len(x)-2*mbc
    my=len(y)-2*mbc
    
#first take care of the ghost cells below and above the domain
    if  grid.y.nstart == 0:
        for i in range(mbc,mbc+mx):
            for ibc in range(1,mbc+1):
                for k in range(0,maux):
                    if mthbcy==2:
                        aux[k,i,mbc-ibc]=aux[k,i,mbc+my-ibc]
                    else:
                        aux[k,i,mbc-ibc]=aux[k,i,mbc]
    
    if grid.y.nend == grid.y.n:
        for i in range(mbc,mbc+mx):
            for ibc in range(1,mbc+1):
                for k in range(0,maux):
                    if mthbcy==2:
                        aux[k,i,mbc-1+my+ibc]=aux[k,i,mbc-1+ibc]
                    else:
                        aux[k,i,mbc-1+my+ibc]=aux[k,i,mbc-1+my]


#now take care of ghost cells at the left and right of the domain
    if grid.x.nstart == 0: 
        for j in range(0,my+2*mbc):
            for ibc in range(1,mbc+1):
                for k in range(0,maux):
                    if mthbcx==2:
                        aux[k,mbc-ibc,j]=aux[k,mbc+mx-ibc,j]
                    else:
                        aux[k,mbc-ibc,j]=aux[k,mbc,j]
    
    if grid.x.nend == grid.x.n:
        for j in range(0,my+2*mbc):
            for ibc in range(1,mbc+1):
                for k in range(0,maux):
                    if mthbcx==2:
                        aux[k,mbc-1+mx+ibc,j]=aux[k,mbc-1+ibc,j]
                    else:
                        aux[k,mbc-1+mx+ibc,j]=aux[k,mbc-1+mx,j]

    return aux

def b4step(solver,solutions):
    grid = solutions['n'].grids[0]
    
    #put in aux[3,:,:] the value of q[0,:,:] (eps). This is required in rptpv.f
    #Uncomment when flux differencing is working. Increase the size of aux accordingly. 
    #grid.aux[3,:,:]=grid.q[0,:,:]
    
    # To set to 0 1st 1/2 of the domain. Used in rect domains with PBC in x
    if grid.aux_global['turnZero_half_2D']==1:
        if grid.t>=grid.aux_global['t_turnZero'] and grid.t<=grid.aux_global['t_turnZero']+1:
            grid.q[:,0:np.floor(grid.aux_global['mx']/2),:]=0
            solutions['n'].grid.q=grid.q

def psystem2D(iplot=False,petscPlot=True,useController=True):
    """
    psystem in 2D with variable coefficients
    """

    from petclaw.grid import Dimension
    from petclaw.grid import Grid
    from pyclaw.solution import Solution
    from petclaw.evolve.solver import PetClawSolver2D
    from pyclaw.controller import Controller
    from petclaw import plot

####################################
######### MAIN PARAMETERS ##########
####################################
# material parameters
    E1=1.;   p1=1.
    E2=10.;   p2=10.
# interface parameters
    alphax=0.5; deltax=1.
    alphay=0.5; deltay=1.
# Linearity parameters
    linearity_mat1=2; linearity_mat2=2
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
# Turning off 1st half of the domain. Useful in rect domains
    turnZero_half_2D=0 #flag
    t_turnZero=0.5
# Regarding time
    nout=4
    tfinal=2.0
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

 #Set global parameters
    grid.aux_global = {}
    grid.aux_global['turnZero_half_2D'] = turnZero_half_2D
    grid.aux_global['t_turnZero'] = t_turnZero
    grid.aux_global['mx'] = mx

# setaux
    xghost=grid.x.centerghost
    yghost=grid.y.centerghost
    grid.aux=setaux(grid,xghost,yghost,mthbc_x_lower,mthbc_y_lower,mbc,E1,p1,E2,p2,linearity_mat1,linearity_mat2,alphax,deltax,alphay,deltay)

# initial conditions
    qinit(grid,A,x0,y0,varx,vary)

# parameters that go to the Riemann Solver
#    grid.aux_global['linearity_mat1']= linearity_mat1
#    grid.aux_global['linearity_mat2']= linearity_mat2
#    from dimsp2 import cparam
#    for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)

    inital_solution = Solution(grid)

    solver = PetClawSolver2D()
    solver.max_steps=1000000
    solver.dt_variable=True
    solver.cfl_max = 0.5
    solver.cfl_desired = 0.45
    solver.order = -1
    solver.mwaves = 2
    solver.fwave = True 
    solver.start_step = b4step
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

    if petscPlot:
        plot.plotPetsc(claw)

    if iplot:
        plot.plotInteractive()

    #eps=claw.frames[claw.nout].grid.gqVec.getArray().reshape([mx,my,grid.meqn])[:,:,0]
    #return eps


if __name__=="__main__":
    psystem2D()
