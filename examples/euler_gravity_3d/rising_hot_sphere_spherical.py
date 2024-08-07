#!/usr/bin/env python
# encoding: utf-8

""" 
Test problem demonstrating 3D hot-sphere rising in an stabilized atmosphere
   in spherical coordinates

This problem evolves the 3D Euler equations using an F-wave
   method, with gravitational source term modifications.
The primary variables are: 
   density (rho), x,y, and z momentum (rho*u,rho*v,rho*w), and energy.
"""
import numpy as np
from clawpack.riemann.mappedGrid import euler3d_mappedgrid as mg

# Test for MPI, and set sizes accordingly
try:
    from mpi4py import MPI
    mpiAvailable = True
except ImportError:
    import warnings
    warnings.warn('mpi4py is not available')
    mpiAvailable = False

if mpiAvailable:
    mpiRank = MPI.COMM_WORLD.Get_rank()
    mpiSize = MPI.COMM_WORLD.Get_size()
else:
    mpiRank = 0
    mpiSize = 1

# Constants
gamma = 1.4 # Ratio of specific heats
gamma1 = gamma - 1.
gR = 980.665 # Acceleration due to gravity [cm/s**2] in Radial Direction
kBoltzmann = 1.3807e-16 # Boltzmann constant [erg/K]
nAvogadro = 6.0221e23 # Avogadro's number [1/mol]
rEarth = 637120000.0 # Radius of Earth [cm]

# Hot Sphere Parameters
rSphere = 40.e5 # Radius of Sphere
TSphere = 2.7e4 # Temperature of Sphere Perturbation
altSphere = rEarth + 150.e5 # Altitude of Sphere
theta0 = np.pi/4.0 # Center of Sphere (theta angle)
phi0 = 0.15 # Center of Sphere (phi angle)
xSphere = altSphere*np.sin(theta0)*np.cos(phi0) # Sphere Center x-coord
ySphere = altSphere*np.sin(theta0)*np.sin(phi0) # Sphere Center y-coord
zSphere = altSphere*np.cos(theta0) # Sphere Center z-coord

# Grid Parameters
mxyz = [160,160,80] # Number of Grid Cells
xyzMin = [theta0-0.15, phi0-0.15, rEarth + 80.0e5  ] # Domain limits (min)
xyzMax = [theta0+0.15, phi0+0.15, rEarth + 950.0e5 ] # Domain limits (max)
mapType = "Spherical"
z0 = xyzMin[2] - rEarth # Bottom of Grid - Distance Above Earth Surface
zN = xyzMax[2] - rEarth # Top of Grid - Distance Above Earth Surface

# Gravity Terms
gravityTerm = True # Turn Gravity Term On or Off in Riemann Solver
gravityEflux = False # Turn Gravity Term in Energy Flux On/Off
gFlux = 0
if gravityEflux: gFlux = 1

#-----------------------------------------------------------------------
# Description:
#   Equilibrium atmosphere
#
# Inputs:
#   p0[mz+2*mbc]      : pressure (1D array)
#   rho0[mz+2*mbc]    : density (1D array)
#   Mavg[mz+2*mbc]    : average molecular mass (1D array)
# 
# Input/Outputs:
#   p0,rho0,Mavg      : 1D z-column initialization of p0 and rho0
#-----------------------------------------------------------------------
def setEquilibriumAtmosphere(p0,rho0,Mavg):

    p0 = [1.28255457e+02,2.45768842e+01,4.14947876e+00,6.29750420e-01,1.01220380e-01,2.64133921e-02,1.22941741e-02,7.08667395e-03,4.52931611e-03,3.07286214e-03,2.16905463e-03,1.57652477e-03,1.17092484e-03,8.84611067e-04,6.77691403e-04,5.25138237e-04,4.10841768e-04,3.24102394e-04,2.57470120e-04,2.05925021e-04,1.65598592e-04,1.33701518e-04,1.08364754e-04,8.82441931e-05,7.21143717e-05,5.91376054e-05,4.86178229e-05,4.00787900e-05,3.30908693e-05,2.73888126e-05,2.27031016e-05,1.88518481e-05,1.56898948e-05,1.30700401e-05,1.08991559e-05,9.09869161e-06,7.60521743e-06,6.36376491e-06,5.32972657e-06,4.46856235e-06,3.74878325e-06,3.14890785e-06,2.64613146e-06,2.22646032e-06,1.87396531e-06,1.57844875e-06,1.33028392e-06,1.12211091e-06,9.47071388e-07,7.99762122e-07,6.75921511e-07,5.71493939e-07,4.83610358e-07,4.09325094e-07,3.46744110e-07,2.93793938e-07,2.49152408e-07,2.11367113e-07,1.79432411e-07,1.52415843e-07,1.29549499e-07,1.10136422e-07,9.37086690e-08,7.97324669e-08,6.79127210e-08,5.78532722e-08,4.93172661e-08,4.20604343e-08,3.58836884e-08,3.06389102e-08,2.61608771e-08,2.23557534e-08,1.91042726e-08,1.63479490e-08,1.39976779e-08,1.19853352e-08,1.02623231e-08,8.78713846e-09,7.53940212e-09,6.46885245e-09,5.55032464e-09,4.76222864e-09,4.09020086e-09,3.51658796e-09]

    rho0 = [1.93347036e-07,4.03984315e-08,7.33795328e-09,1.16964004e-09,1.64049100e-10,2.53990286e-11,7.54287116e-12,3.40478277e-12,1.84556481e-12,1.10964372e-12,7.13581470e-13,4.81506393e-13,3.36472592e-13,2.41540079e-13,1.77156053e-13,1.32213794e-13,1.00089557e-13,7.67024111e-14,5.93930647e-14,4.64294817e-14,3.65782332e-14,2.90138753e-14,2.31378048e-14,1.85800114e-14,1.49929512e-14,1.21526733e-14,9.89015561e-15,8.07840567e-15,6.61976992e-15,5.43890503e-15,4.48202167e-15,3.70250573e-15,3.06590093e-15,2.54266886e-15,2.11283102e-15,1.75827860e-15,1.46560471e-15,1.22337830e-15,1.02239821e-15,8.55585508e-16,7.16578299e-16,6.01033981e-16,5.04419184e-16,4.23940996e-16,3.56468062e-16,2.99992883e-16,2.52633808e-16,2.12955966e-16,1.79630105e-16,1.51610996e-16,1.28075790e-16,1.08244792e-16,9.15665290e-17,7.74771188e-17,6.56137471e-17,5.55805979e-17,4.71251502e-17,3.99708405e-17,3.39261636e-17,2.88137888e-17,2.44878021e-17,2.08159094e-17,1.77092661e-17,1.50666724e-17,1.28321441e-17,1.09306468e-17,9.31730480e-18,7.94587120e-18,6.77866202e-18,5.78764327e-18,4.94156316e-18,4.22266806e-18,3.60840539e-18,3.08771188e-18,2.64374425e-18,2.26362608e-18,1.93817162e-18,1.65953699e-18,1.42386938e-18,1.22167290e-18,1.04819271e-18,8.99349679e-19,7.72429901e-19,6.64098458e-19]

    Mavg = [28.85614554,28.85337155,28.83817654,28.56226512,27.60224909,26.26692289,25.23573593,24.45469565,23.79308533,23.18781005,22.61490394,22.07318988,21.55703223,21.06778441,20.60540309,20.17202267,19.76585711,19.38847601,19.0408475, 18.71970337,18.42758099,18.16274099,17.92359740,17.70606183,17.51035814,17.33530373,17.17893585,17.03979933,16.91620578,16.80712079,16.71028376,16.62471452,16.54940299,16.48292773,16.42454596,16.37307369,16.32776306,16.28801338,16.2531155, 16.22247335,16.19551611,16.17188138,16.15108306,16.13288090,16.11686426,16.10282002,16.09046507,16.07960946,16.07007411,16.06169374,16.05433222,16.04784993,16.04215209,16.03712679,16.0327204,16.02883120,16.02540929,16.02239140,16.01973516,16.01738918,16.01531699,16.01348647,16.01187781,16.01045286,16.00919766,16.00808580,16.00710454,16.00623687,16.00546792,16.00478755,16.00418349,16.00365220,16.00317996,16.00276269,16.00239247,16.00206303,16.00176987,16.00150902,16.00127962,16.00107519,16.00089299,16.00073063,16.00058692,16.00045964]

    return p0,rho0,Mavg

#-----------------------------------------------------------------------
# Description:
#   Modify pressure to create numeric atmosphere equilibrium
#
# Inputs:
#   ze0[mz+2*mbc+1]     : cell edge grid values
#   p0[mz+2*mbc]        : pressure
#   rho0[mz+2*mbc]      : density
#
# Input/Outputs:
#   p0,rho0           : 1D z-column modification of p0 and rho0
#-----------------------------------------------------------------------
def modifyEquilibriumAtmosphere(zep0,p0,rho0):

    # Compute the delta-z (dz)
    nz = np.size(zep0)-1
    dz = np.zeros([nz],dtype='float',order='F')
    for iz in range(nz-1):
        dz[iz] = zep0[iz+1]-zep0[iz]

    # Compute modified pressure at cell centers
    iz = nz-1
    dz2 = (dz[iz]+dz[iz-1])*0.5
    p0[iz] = p0[iz] + rho0[iz]*gR*dz2
    for iz in range(nz-1,0,-1):
        dz2 = (dz[iz]+dz[iz-1])*0.5
        finterp = dz[iz-1]/(dz[iz]+dz[iz-1])
        rho_b = rho0[iz]*finterp + rho0[iz-1]*(1.-finterp)
        p0[iz-1] = p0[iz] + rho_b*gR*dz2

    return p0

#-----------------------------------------------------------------------
# Description:
#   Custom BCs for the z-direction
#-----------------------------------------------------------------------
def customBCLowerZ(state,dim,t,qbc,auxbc,mbc):
    for k in range(mbc):
        rZ = np.sqrt(xcpZ[k]**2 + ycpZ[k]**2 + zcpZ[k]**2)-rEarth
        qbc[0,:,:,k] =  rho0[k]
        qbc[1,:,:,k] =  0.
        qbc[2,:,:,k] =  0.
        qbc[3,:,:,k] =  0.
        qbc[4,:,:,k] =  p0[k]/gamma1 + qbc[0,:,:,k]*gR*rZ*gFlux

def customBCUpperZ(state,dim,t,qbc,auxbc,mbc):
    for k in range(mbc):
        rZ = np.sqrt(xcpZ[-k-1]**2 + ycpZ[-k-1]**2 + zcpZ[-k-1]**2)-rEarth
        qbc[0,:,:,-k-1] = rho0[-k-1]
        qbc[1,:,:,-k-1] = qbc[1,:,:,-mbc-1]
        qbc[2,:,:,-k-1] = qbc[2,:,:,-mbc-1]
        qbc[3,:,:,-k-1] = qbc[3,:,:,-mbc-1]
        rhov2 = (qbc[1,:,:,-k-1]**2 + qbc[2,:,:,-k-1]**2 + qbc[3,:,:,-k-1]**2)/qbc[0,:,:,-k-1]
        qbc[4,:,:,-k-1] = p0[-k-1]/gamma1 + 0.5*rhov2 + qbc[0,:,:,-k-1]*gR*rZ*gFlux

def customAuxBCLowerZ(state,dim,t,qbc,auxbc,mbc):
    auxbc[:,:,:,:mbc] = auxtmp[:,:,:,:mbc]

def customAuxBCUpperZ(state,dim,t,qbc,auxbc,mbc):
    auxbc[:,:,:,-mbc:] = auxtmp[:,:,:,-mbc:]

#-----------------------------------------------------------------------
#  Main script for solving 3D Euler equations using Clawpack/PyClaw.
#-----------------------------------------------------------------------
def euler3d(kernel_language='Fortran',solver_type='classic',\
            use_petsc=False,outdir='./_output',\
            output_format='hdf5',file_prefix='equil',disable_output=False,\
            mx=mxyz[0],my=mxyz[1],mz=mxyz[2],\
            tfinal=64.0,num_output_times=1):

    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='classic':
        solver = pyclaw.ClawSolver3D()
        solver.dimensional_split = True
        solver.transverse_waves = 20
        solver.limiters = pyclaw.limiters.tvd.minmod
        solver.num_ghost = 2
        solver.order = 2
        solver.fwave = True
    elif solver_type=='sharpclaw':
        solver = pyclaw.SharpClawSolver3D()
    else:
        raise Exception('Unrecognized solver_type.')

    import logging
    solver.logger.setLevel(logging.DEBUG)
    from clawpack import riemann
    solver.rp = riemann.euler_mapgrid_3D
    solver.num_eqn = 5
    solver.num_waves = 3
    solver.cfl_max = 0.6
    solver.cfl_desired = 0.5
    solver.dt_initial = 1.e-0
    solver.max_steps = 10000

    # Initialize Domain
    x = pyclaw.Dimension(0.0,1.0,mx,name='x')
    y = pyclaw.Dimension(0.0,1.0,my,name='y')
    z = pyclaw.Dimension(0.0,1.0,mz,name='z')
    domain = pyclaw.Domain([x,y,z])

    num_aux = 15
    state = pyclaw.State(domain,solver.num_eqn,num_aux)
    state.problem_data['gamma'] = gamma
    state.problem_data['g_r'] = gR
    state.problem_data['gravity'] = gravityTerm
    state.problem_data['gravityflux'] = gravityEflux

    # Grids
    mbc = solver.num_ghost
    grid = state.grid

    # Computational Grid Sizes
    dxc = domain.grid.delta[0]
    dyc = domain.grid.delta[1]
    dzc = domain.grid.delta[2]
    pmx, pmy, pmz = grid.num_cells[0], grid.num_cells[1], grid.num_cells[2]

    # Computational Grid Centers and nodes
    centers = grid.c_centers # centers (Comp.)
    centersBC = grid.c_centers_with_ghost(mbc) # centers w Ghost (Comp.)
    nodesBC = grid.c_nodes_with_ghost(mbc) # nodes w Ghost (Comp.)

    # Grid Centers Without Boundary Cells (1D Slice) - Comp. and Phys.
    xcc = grid.x.centers # x centers (Comp.)
    ycc = grid.y.centers # y centers (Comp.)
    zcc = grid.z.centers # z centers (Comp.)
    xcp,ycp,zcp = mg.mapc2pwrapper(xcc,ycc,zcc,pmz,xyzMin,xyzMax,mapType)
    
    # Grid Centers Without Boundary Cells (3D Arrays)
    Xcc,Ycc,Zcc = centers[0][:][:][:],centers[1][:][:][:],centers[2][:][:][:]
    Xcp,Ycp,Zcp = mg.mapc2pwrapper(Xcc,Ycc,Zcc,pmz,xyzMin,xyzMax,mapType)
    Xcp = np.reshape(Xcp,[pmx,pmy,pmz],order='F') # x centers (Phys.)
    Ycp = np.reshape(Ycp,[pmx,pmy,pmz],order='F') # y centers (Phys.)
    Zcp = np.reshape(Zcp,[pmx,pmy,pmz],order='F') # z centers (Phys.)

    # Grid nodes With Boundary Cells (1D Slice along z)- Comp. and Phys.
    xecZ = nodesBC[0][0][0][:] # x nodes along z (Comp.)
    yecZ = nodesBC[1][0][0][:] # y nodes along z (Comp.)
    zecZ = nodesBC[2][0][0][:] # z nodes along z (Comp.)
    xepZ,yepZ,zepZ = mg.mapc2pwrapper(xecZ,yecZ,zecZ,pmz,xyzMin,xyzMax,mapType)

    # Grid Centers With Boundary Cells (1D Slice along z) - Comp. and Phys.
    global xcpZ, ycpZ, zcpZ
    xccZ = centersBC[0][0][0][:] # x centers along z (Comp.)
    yccZ = centersBC[1][0][0][:] # y centers along z (Comp.)
    zccZ = centersBC[2][0][0][:] # z centers along z (Comp.)
    xcpZ,ycpZ,zcpZ = mg.mapc2pwrapper(xccZ,yccZ,zccZ,pmz,xyzMin,xyzMax,mapType)
 
    if np.sqrt(xepZ[0]**2+yepZ[0]**2+zepZ[0]**2)-rEarth <= 0:
        print("WARNING: z may go below Earth's surface"," zepZ: ",zepZ[0:10])

    #-----------------------------------------------------------------------
    # Create vectors for 1D pressure and density column with boundary cells
    #-----------------------------------------------------------------------
    mz0 = pmz+2*mbc
    global p0, rho0, Mavg
    p0 = np.zeros([mz0],dtype='float',order='F')
    rho0 = np.zeros([mz0],dtype='float',order='F')
    Mavg = np.zeros([mz0],dtype='float',order='F')

    # Set the equilibrium pressure such that dp/dz = -rho*gR
    p0,rho0,Mavg = setEquilibriumAtmosphere(p0,rho0,Mavg)

    # Modify the equilibrium such that dp/dz = -rho*gR is held numerically
    altEdgesAboveEarth = np.sqrt(xepZ**2 + yepZ**2 + zepZ**2) - rEarth
    p0 = modifyEquilibriumAtmosphere(altEdgesAboveEarth,p0,rho0)

    # Set the auxiliary variables
    xlower,ylower,zlower = nodesBC[0][0][0][0],nodesBC[1][0][0][0],nodesBC[2][0][0][0]
    dxc,dyc,dzc = domain.grid.delta[0],domain.grid.delta[1],domain.grid.delta[2]

    global auxtmp
    auxtmp = np.zeros([num_aux,pmx+2*mbc,pmy+2*mbc,pmz+2*mbc],dtype='float',order='F')
    auxtmp = mg.setauxiliaryvariables(num_aux,mbc,pmx,pmy,pmz,xlower,ylower,zlower,dxc,dyc,dzc,xyzMin,xyzMax,mapType)
    state.aux[:,:,:,:] = auxtmp[:,mbc:-mbc,mbc:-mbc,mbc:-mbc]

    # Set Index for Capcaity Function in state.aux (Python 0-based)
    state.index_capa = 12 

    # Set the state variables (Initial Conditions)

    # Initialize p,T,velSqrd
    p = np.zeros([pmx,pmy,pmz],dtype='float',order='F')
    T = np.zeros([pmx,pmy,pmz],dtype='float',order='F')
    velSqrd = np.zeros([pmx,pmy,pmz],dtype='float',order='F')

    # Density
    for i in range(np.size(p,0)):
        for j in range(np.size(p,1)):
            # NEEDS TO BE FIXED WHEN MPI SLICES NORMAL TO Z
            state.q[0,i,j,:] = rho0[mbc:pmz+mbc] 
    
    # Momentum
    state.q[1,:,:,:] = 0. # x-momentum (rho*u)
    state.q[2,:,:,:] = 0. # y-momentum (rho*v)
    state.q[3,:,:,:] = 0. # z-momentum (rho*w)

    # Velocity Squared (u**2+v**2+w**2)
    velSqrd[:,:,:] = (state.q[1,:,:,:]**2+state.q[2,:,:,:]**2+state.q[3,:,:,:]**2)/state.q[0,:,:,:]**2

    # Energy
    for i in range(np.size(p,0)):
        for j in range(np.size(p,1)):
            # NEEDS TO BE FIXED WHEN MPI SLICES NORMAL TO Z
            p[i,j,:] = p0[mbc:pmz+mbc]
    state.q[4,:,:,:] = p/gamma1 + 0.5*state.q[0,:,:,:]*velSqrd + gFlux*state.q[0,:,:,:]*gR*(np.sqrt(Xcp[:,:,:]**2+Ycp[:,:,:]**2+Zcp[:,:,:]**2)-rEarth)

    #-----------------------------------------------------------------------
    # Add Perturbation
    #-----------------------------------------------------------------------
    # Define Temperature for Perturbation
    T = p/state.q[0,:,:,:]
    L = np.sqrt((Xcp-xSphere)**2+(Ycp-ySphere)**2+(Zcp-zSphere)**2)
    for i in range(pmx):
        for j in range(pmy):
            for k in range(pmz):
                if L[i,j,k] <= rSphere:
                    # Temperature Perturbation
                    mu = Mavg[k+mbc]/nAvogadro
                    T[i,j,k] += TSphere*(kBoltzmann/mu)*(1.0-L[i,j,k]/rSphere)
                    p[i,j,k] = T[i,j,k]*state.q[0,i,j,k]
    state.q[4,:,:,:] = p/gamma1 + 0.5*state.q[0,:,:,:]*velSqrd + gFlux*state.q[0,:,:,:]*gR*(np.sqrt(Xcp[:,:,:]**2+Ycp[:,:,:]**2+Zcp[:,:,:]**2)-rEarth)

    # Setup Boundary Conditions

    # X - Boundary Conditions
    solver.bc_lower[0] = pyclaw.BC.extrap
    solver.bc_upper[0] = pyclaw.BC.extrap

    # Y - Boundary Conditions
    solver.bc_lower[1] = pyclaw.BC.extrap
    solver.bc_upper[1] = pyclaw.BC.extrap

    # Z - Boundary Conditions
    solver.bc_lower[2] = pyclaw.BC.custom
    solver.bc_upper[2] = pyclaw.BC.custom
    solver.user_bc_lower = customBCLowerZ
    solver.user_bc_upper = customBCUpperZ

    # Aux - Boundary Conditions
    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[1] = pyclaw.BC.extrap
    solver.aux_bc_upper[1] = pyclaw.BC.extrap
    solver.aux_bc_lower[2] = pyclaw.BC.custom
    solver.aux_bc_upper[2] = pyclaw.BC.custom
    solver.user_aux_bc_lower = customAuxBCLowerZ
    solver.user_aux_bc_upper = customAuxBCUpperZ

    # Solver Parameters
    claw = pyclaw.Controller()
    claw.verbosity = 4
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.output_format = output_format
    claw.output_file_prefix = file_prefix
    claw.keep_copy = False
    if disable_output:
        claw.output_format = None
    claw.tfinal = tfinal
    claw.num_output_times = num_output_times
    claw.outdir = outdir

    return claw

# __main__()
if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(euler3d)
