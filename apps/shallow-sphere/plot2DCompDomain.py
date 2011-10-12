# Import some libraries
import pyclaw
import shallow-4-Rossby-Haurwitz-wave as rh
import numpy as np

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt


# Nondimensionalized radius of the earth
Rsphere = 1.0

def contourLineSphere(fileName='fort.q0000',path='./_output'):
    """
    This function plot the contour lines on a spherical surface for the shallow
    water equations solved on a sphere.
    """
      
    # Concatenate path and file name
    pathFileName = path + "/" + fileName

    # Open file
    f = file(pathFileName,"r")

    # Read first six lines. They are not used later!
    unsed = f.readline()  # grid_number
    unused = f.readline() # AMR_level

    # Read mx, my, xlow, ylow, dx and dy
    line = f.readline()
    sline = line.split()
    mx = int(sline[0])

    line = f.readline()
    sline = line.split()
    my = int(sline[0])

    line = f.readline()
    sline = line.split()
    xlower = float(sline[0])

    line = f.readline()
    sline = line.split()
    ylower = float(sline[0])


    line = f.readline()
    sline = line.split()
    dx = float(sline[0])

    line = f.readline()
    sline = line.split()
    dy = float(sline[0])


    # Grid:
    # ====
    xupper = xlower + mx * dx
    yupper = ylower + my * dy

    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    grid = pyclaw.Grid([x,y])


    # Override default mapc2p function
    # ================================
    grid.mapc2p = rh.mapc2p_sphere_vectorized

    # Compute the physical coordinates of each cell's center
    # ======================================================
    grid.compute_p_center(recompute=True)
    xp = grid._p_center[0]
    yp = grid._p_center[1]
    zp = grid._p_center[2]

    grid.compute_c_center(recompute=True)
    xc = grid._c_center[0]
    yc = grid._c_center[1]
    
    # Define arrays of conserved variables
    h = np.zeros((mx,my))
    hu = np.zeros((mx,my))
    hv = np.zeros((mx,my))
    hw = np.zeros((mx,my))
    
    # Read solution
    for j in range(my):
        tmp = np.fromfile(f,dtype='float',sep=" ",count=4*mx)
        tmp = tmp.reshape((mx,4))
        h[:,j] = tmp[:,0]
        hu[:,j] = tmp[:,1]
        hv[:,j] = tmp[:,2]
        hw[:,j] = tmp[:,3]


    # Plot solution in the computational domain
    # #########################################

    # Fluid height
    plt.figure()
    CS = plt.contour(xc,yc,h)
    plt.title('Fluid height (computational domain)')
    plt.xlabel('xc')
    plt.ylabel('yc')
    plt.clabel(CS, inline=1, fontsize=10)
    plt.show()


if __name__=="__main__":
    from pyclaw.util import run_app_from_main
    output = run_app_from_main(contourLineSphere)
    print 'Error: ',output
