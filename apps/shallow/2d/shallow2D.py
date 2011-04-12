#!/usr/bin/env python
# encoding: utf-8

"""
2D shallow water equations.
"""
#===========================================================================
# Import libraries
#===========================================================================

import numpy as np
import pdb


def qinit(grid,hl,ul,vl,hr,ur,vr,radDam):
    # Create an array with fortran native ordering
    xCenter = grid.x.center
    yCenter = grid.y.center
    Y,X = np.meshgrid(yCenter,xCenter)

    q=np.empty([grid.meqn,len(xCenter),len(yCenter)], order = 'F')
    for i,x in enumerate(xCenter):
        for j,y in enumerate(yCenter):
            #r = np.sqrt(xCenter[i]**2 + yCenter[j]**2)
            if x<0.: 
                q[0,i,j] = hl
                q[1,i,j] = hl*ul
                q[2,i,j] = hl*vl
            else: 
                q[0,i,j] = hr
                q[1,i,j] = hr*ur
                q[2,i,j] = hr*vr
    
    grid.q=q



#def setaux(maxlength):
#    """
#    This function "allocate" memory for some arrays that are needed in both rpn2sw.f
#    and rpt2sw.f Riemann solver. These arrays are not really aux variable in the strict sense.
#    Therefore, this solution might not be the optimal one.
#    """
#    
#    aux=np.zeros([4,maxlength],order='F')
#
#    return aux



    
def shallow2D(iplot=False,petscPlot=False,useController=True,htmlplot=True):
    #===========================================================================
    # Import libraries
    #===========================================================================
    from petsc4py import PETSc   
    from petclaw.grid import Grid
    from petclaw.grid import Dimension
    from pyclaw.solution import Solution
    from petclaw.evolve.clawpack import PetClawSolver2D
    from pyclaw.controller import Controller
    from petclaw import plot

    #===========================================================================
    # Initialize grids, then initialize the solution associated to the grid and
    # finally initialize aux array
    #===========================================================================

    # Grid:
    xlower = -2.5
    xupper = 2.5
    mx = 100
    
    ylower = -2.5
    yupper = 2.5
    my = 100
    x = Dimension('x',xlower,xupper,mx,mthbc_lower=1,mthbc_upper=1)
    y = Dimension('y',ylower,yupper,my,mthbc_lower=1,mthbc_upper=1)
    grid = Grid([x,y])

    grid.meqn = 3  # Number of equations
    grid.mbc = 3   # Number of ghost cells
    grid.t = 0.    # Initial time

    # Parameters
    #grid.aux_global['grav'] = 1.0
    #from step1 import cparam
    #for key,value in grid.aux_global.iteritems(): setattr(cparam,key,value)

    # Initial solution
    # ================
    # Initialize petsc Structures for q. 
    # This must be called before grid.x.center and such can be accessed.
    grid.init_q_petsc_structures()

    radDam = 0.5
    hl = 2.
    ul = 0.
    vl = 0.
    hr = 1.
    ur = 0.
    vr = 0.
    
    qinit(grid,hl,ul,vl,hr,ur,vl,radDam) # This function is defined above

    # aux array
#    xCenter = grid.x.center
#    yCenter = grid.y.center
#    maxNbrCells = max(len(xCenter),len(yCenter))
#    grid.aux = setaux(maxNbrCells)

    init_solution = Solution(grid)



    #===========================================================================
    # Setup solver and solver parameters
    #===========================================================================
    kernelsType = 'F'
    solver = PetClawSolver2D(kernelsType = kernelsType)
    solver.mwaves = 3
    if kernelsType =='P': solver.set_riemann_solver('shallow_roe')
    solver.mthlim = [4]*solver.mwaves



    #===========================================================================
    # Setup controller and controller paramters
    #===========================================================================
    claw = Controller()
    claw.keep_copy = True
    claw.output_format = 'petsc' # The output format MUST be set to petsc!!
    claw.tfinal = 2.0
    claw.solutions['n'] = init_solution
    claw.solver = solver


    #===========================================================================
    # Solve the problem
    #===========================================================================
    status = claw.run()


    pdb.set_trace()

    if htmlplot: plot.plotHTML()
    if petscPlot: plot.plotPetsc(output_object)
    if iplot: plot.plotInteractive()



if __name__=="__main__":
    import sys
    from petclaw.util import _info_from_argv
    args, kwargs = _info_from_argv(sys.argv)
    error=shallow2D(*args,**kwargs)
    print 'Error: ',error






