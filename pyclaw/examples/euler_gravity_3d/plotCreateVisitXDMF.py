#!/usr/bin/env python
# encoding: utf-8
""" 
!-----------------------------------------------------------------------
! Description: Creates XDMF (.xmf) files to go along with clawpack .hdf
!   files for use visualizing within VisIt (https://visit.llnl.gov).
!-----------------------------------------------------------------------
"""
from __future__ import absolute_import
from __future__ import print_function
import os,sys
import h5py
from numpy import *
from visitHDF5XDMF import *
sys.path.append(os.getcwd())
set_printoptions(precision=14)

from rising_hot_sphere import *
#from rising_hot_sphere_spherical import *

# Parameters
neqn = 5
var_names = ['rho','rhou','rhov','rhow','energy']
mesh_name = 'mesh.h5'
out_dir = './_output'
file_prefix = 'equil'
mx = mxyz[0]
my = mxyz[1]
mz = mxyz[2]

# Create Grid
xc,yc,zc = linspace(0.0,1.0,mx),linspace(0.0,1.0,my),linspace(0.0,1.0,mz)
Xc,Yc,Zc = meshgrid(xc,yc,zc,indexing='ij')
Xp,Yp,Zp = mg.mapc2pwrapper(Xc,Yc,Zc,mz,xyzMin,xyzMax,mapType)
Xp,Yp,Zp = reshape(Xp,[mx,my,mz],order='F'),reshape(Yp,[mx,my,mz],order='F'),reshape(Zp,[mx,my,mz],order='F')

# Make list of hdf files (assumes these are the output)
import glob
cwd = os.getcwd()
os.chdir(out_dir)
files = sorted(glob.glob('*.hdf'))
os.chdir(cwd)

# Read First File Attributes
filename = out_dir + '/' + files[0]
print(filename)
file0 = h5py.File(filename,'r')
dset0 = file0['/patch1']
coord_names = []
for coord in dset0.attrs['dimensions']:
    coord_names.append(coord)
ndims = size(coord_names)

# Write Grid to HDF5 File
writeHDF53DS(Xp,Yp,Zp,coord_names,mesh_name,out_dir)

# Write XMF Files for Each Output File
for fname in files:
    f = h5py.File(out_dir+'/'+fname,'r')
    dset = f['/patch1']
    time = float(dset.attrs['t'])
    print("Writing: ",fname[:-4] + '.xmf')
    writeXDMF3DSMeshHyper(mx,my,mz,var_names,fname,mesh_name,out_dir,time)
