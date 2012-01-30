#!/usr/bin/env python
# encoding: utf-8

"""
This script plots the solution of the shallow water on a sphere in the 
rectangular computational domain. The user can specify the name of the solution
file and its path. If these two information are not given, the script checks 
whether the solution fort.q0000 in ./_output exist and plots it. If it it does
not exist a error message is printed at screen.
The file must be ascii and clawpack format.

This function shows how to read and plot the solution stored in an ascii file 
written by pyclaw.
"""

# Import some libraries
import pyclaw
import shallow_4_Rossby_Haurwitz_wave as rh
import numpy as np

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt

import sys

# Nondimensionalized radius of the earth
Rsphere = 1.0

def contourLineSphere(fileName='fort.q0000',path='./_output'):
    """
    This function plot the contour lines on a spherical surface for the shallow
    water equations solved on a sphere.
    """  
 
    # Open file
    # =========
    
    # Concatenate path and file name
    pathFileName = path + "/" + fileName

    try:
        f = file(pathFileName,"r")
    except IOError as e:
        print("({})".format(e))
        sys.exit()


    # Read file header
    # ================
    # The information contained in the first two lines are not used.
    unsed = f.readline()  # patch_number
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


    # Patch:
    # ====
    xupper = xlower + mx * dx
    yupper = ylower + my * dy

    x = pyclaw.Dimension('x',xlower,xupper,mx)
    y = pyclaw.Dimension('y',ylower,yupper,my)
    patch = pyclaw.Patch([x,y])


    # Override default mapc2p function
    # ================================
    patch.mapc2p = rh.mapc2p_sphere_vectorized  


    # Compute the physical coordinates of each cell's centers
    # ======================================================
    patch.compute_p_centers(recompute=True)
    xp = patch._p_centers[0]
    yp = patch._p_centers[1]
    zp = patch._p_centers[2]

    patch.compute_c_centers(recompute=True)
    xc = patch._c_centers[0]
    yc = patch._c_centers[1]
    
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
    # =========================================

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
