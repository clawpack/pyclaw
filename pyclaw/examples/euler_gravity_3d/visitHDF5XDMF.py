#!/usr/bin/env python
# encoding: utf-8

""" 
!-----------------------------------------------------------------------
! HDF5 and XDMF routines for creating output files for VisIt.
!-----------------------------------------------------------------------
"""
from __future__ import absolute_import
import os,sys
import h5py
from numpy import *

# Write a 2D HDF5 File for Rectilinear Mesh
def writeHDF52D(x,y,filename,out_dir):
    # Write H5 Grid Data
    h5_fname = out_dir + '/' + filename
    f = h5py.File(h5_fname,'w')
    dset = f.create_dataset("x",(x.size,),'f8')
    dset[...] = x
    dset = f.create_dataset("y",(y.size,),'f8')
    dset[...] = y
    f.close()
        
# Write a 3D HDF5 File for Rectilinear Mesh
def writeHDF53D(x,y,z,filename,out_dir):
    # Write H5 Grid Data
    h5_fname = out_dir + '/' + filename
    f = h5py.File(h5_fname,'w')
    dset = f.create_dataset("x",(x.size,),'f8')
    dset[...] = x
    dset = f.create_dataset("y",(y.size,),'f8')
    dset[...] = y
    dset = f.create_dataset("z",(z.size,),'f8')
    dset[...] = z
    f.close()

# Write a 3D HDF5 File for Curvilinear Mesh
def writeHDF53DS(x,y,z,names,filename,out_dir):
    # Write H5 Grid Data
    h5_fname = out_dir + '/' + filename
    f = h5py.File(h5_fname,'w')
    dset = f.create_dataset(names[0],(size(x,0),size(y,1),size(z,2),),'f8')
    dset[:] = x
    dset = f.create_dataset(names[1],(size(x,0),size(y,1),size(z,2),),'f8')
    dset[:] = y
    dset = f.create_dataset(names[2],(size(x,0),size(y,1),size(z,2),),'f8')
    dset[:] = z
    f.close()

# Write a XDMF (.xmf) File for 2D Rectilinear Mesh
def writeXDMF2D(t0,nx,ny,filename,mesh_name,out_dir,simTime):
    # Write XDMF (.xmf) file
    prefix = filename[:-7]
    xmfname = out_dir + '/' + filename[:-3] + '.xmf'
    f = open(xmfname,'w')
    f.write('<?xml version="1.0" ?> \n')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []> \n')
    f.write('<Xdmf Version="2.0"> \n')
    f.write(' <Domain> \n')
    f.write(' \n')
    f.write('   <!-- Logical Mesh  -->\n')
    f.write(' \n')
    f.write('   <Grid Name="mesh2D" GridType="Uniform">\n')
    f.write('     <Time TimeType="single" Value="'+str(simTime)+'"/>\n')
    f.write('     <Topology TopologyType="2DRectMesh" NumberOfElements="'+str(nx)+' '+str(ny)+'"/>\n')
    f.write('     <Geometry GeometryType="VxVy">\n')
    f.write('       <DataItem Dimensions="'+str(nx)+'" NumberType="Float" Precision="8" Format="HDF">\n')
    f.write('         '+mesh_name+':/x \n')
    f.write('       </DataItem>\n')
    f.write('       <DataItem Dimensions="'+str(ny)+'" NumberType="Float" Precision="8" Format="HDF">\n')
    f.write('         '+mesh_name+':/y \n')
    f.write('       </DataItem>\n')
    f.write('     </Geometry>\n')
    f.write('     <Attribute Name="'+prefix[:-1]+'" AttributeType="Scalar" Center="Node">\n')
    f.write('       <DataItem Dimensions="'+str(nx)+' '+str(ny)+'" NumberType="Float" Precision="8" Format="HDF">\n')
    f.write('         '+filename+':/'+prefix[:-1]+' \n')
    f.write('       </DataItem>\n')
    f.write('     </Attribute>\n')
    f.write('   </Grid>\n')
    f.write(' </Domain> \n')
    f.write('</Xdmf>')
    f.close()

# Write a XDMF (.xmf) File for 3D Rectilinear Mesh
def writeXDMF3D(nx,ny,nz,var_names,filename,mesh_name,out_dir,simTime):
    # Write XDMF (.xmf) file
    xmfname = out_dir + '/' + filename[:-3] + '.xmf'
    f = open(xmfname,'w')
    f.write('<?xml version="1.0" ?> \n')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []> \n')
    f.write('<Xdmf Version="2.0"> \n')
    f.write(' <Domain> \n')
    f.write(' \n')
    f.write('   <!-- Logical Mesh  -->\n')
    f.write(' \n')
    f.write('   <Grid Name="mesh" GridType="Uniform">\n')
    f.write('     <Time TimeType="single" Value="'+str(simTime)+'"/>\n')
    f.write('     <Topology TopologyType="3DRectMesh" NumberOfElements="'+str(nz)+' '+str(ny)+' '+str(nx)+'"/>\n')
    f.write('     <Geometry GeometryType="VxVyVz">\n')
    f.write('       <DataItem Dimensions="'+str(nx)+'" NumberType="Float" Precision="8" Format="HDF">\n')
    f.write('         '+mesh_name+':/x \n')
    f.write('       </DataItem>\n')
    f.write('       <DataItem Dimensions="'+str(ny)+'" NumberType="Float" Precision="8" Format="HDF">\n')
    f.write('         '+mesh_name+':/y \n')
    f.write('       </DataItem>\n')
    f.write('       <DataItem Dimensions="'+str(nz)+'" NumberType="Float" Precision="8" Format="HDF">\n')
    f.write('         '+mesh_name+':/z \n')
    f.write('       </DataItem>\n')
    f.write('     </Geometry>\n')
    for var_name in var_names:
        f.write('     <Attribute Name="'+var_name+'" AttributeType="Scalar" Center="Node">\n')
        f.write('       <DataItem Dimensions="'+str(nz)+' '+str(ny)+' '+str(nx)+'" NumberType="Float" Precision="8" Format="HDF">\n')
        f.write('         '+filename+':/'+var_name+' \n')
        f.write('       </DataItem>\n')
        f.write('     </Attribute>\n')
    f.write('   </Grid>\n')
    f.write(' </Domain> \n')
    f.write('</Xdmf>')
    f.close()

# Write a XDMF (.xmf) File for 3D Curvilinear Mesh
def writeXDMF3DSMesh(nx,ny,nz,var_names,filename,mesh_name,out_dir,simTime):
    # Write XDMF (.xmf) file
    xmfname = out_dir + '/' + filename[:-3] + '.xmf'
    f = open(xmfname,'w')
    f.write('<?xml version="1.0" ?> \n')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []> \n')
    f.write('<Xdmf Version="2.0"> \n')
    f.write(' <Domain> \n')
    f.write(' \n')
    f.write('   <!-- Logical Mesh  -->\n')
    f.write(' \n')
    f.write('   <Grid Name="mesh" GridType="Uniform">\n')
    f.write('     <Time TimeType="single" Value="'+str(simTime)+'"/>\n')
    f.write('     <Topology TopologyType="3DSMesh" NumberOfElements="'+str(nx)+' '+str(ny)+' '+str(nz)+'"/>\n')
    f.write('     <Geometry GeometryType="X_Y_Z">\n')
    f.write('       <DataItem Dimensions="'+str(nx)+' '+str(ny)+' '+str(nz)+'" NumberType="Float" Precision="8" Format="HDF">\n')
    f.write('         '+mesh_name+':/x \n')
    f.write('       </DataItem>\n')
    f.write('       <DataItem Dimensions="'+str(nx)+' '+str(ny)+' '+str(nz)+'" NumberType="Float" Precision="8" Format="HDF">\n')
    f.write('         '+mesh_name+':/y \n')
    f.write('       </DataItem>\n')
    f.write('       <DataItem Dimensions="'+str(nx)+' '+str(ny)+' '+str(nz)+'" NumberType="Float" Precision="8" Format="HDF">\n')
    f.write('         '+mesh_name+':/z \n')
    f.write('       </DataItem>\n')
    f.write('     </Geometry>\n')
    for var_name in var_names:
        f.write('     <Attribute Name="'+var_name+'" AttributeType="Scalar" Center="Node">\n')
        #f.write('       <DataItem Dimensions="'+str(nz*ny*nx)+'" NumberType="Float" Precision="8" Format="HDF">\n')
        f.write('       <DataItem Dimensions="'+str(nx)+' '+str(ny)+' '+str(nz)+'" NumberType="Float" Precision="8" Format="HDF">\n')
        f.write('         '+filename+':/'+var_name+' \n')
        f.write('       </DataItem>\n')
        f.write('     </Attribute>\n')
    f.write('   </Grid>\n')
    f.write(' </Domain> \n')
    f.write('</Xdmf>')
    f.close()

# Write a XDMF (.xmf) File for 3D Rectilinear Mesh with HyperSlab Data 
def writeXDMF3DSMeshHyper(nx,ny,nz,var_names,filename,mesh_name,out_dir,simTime):
    # Write XDMF (.xmf) file
    xmfname = out_dir + '/' + filename[:-4] + '.xmf'
    f = open(xmfname,'w')
    f.write('<?xml version="1.0" ?> \n')
    f.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []> \n')
    f.write('<Xdmf Version="2.0"> \n')
    f.write(' <Domain> \n')
    f.write(' \n')
    f.write('   <!-- Logical Mesh  -->\n')
    f.write(' \n')
    f.write('   <Grid Name="mesh" GridType="Uniform">\n')
    f.write('     <Time TimeType="single" Value="'+str(simTime)+'"/>\n')
    f.write('     <Topology TopologyType="3DSMesh" NumberOfElements="'+str(nx)+' '+str(ny)+' '+str(nz)+'"/>\n')
    f.write('     <Geometry GeometryType="X_Y_Z">\n')
    f.write('       <DataItem Dimensions="'+str(nx)+' '+str(ny)+' '+str(nz)+'" NumberType="Float" Precision="8" Format="HDF">\n')
    f.write('         '+mesh_name+':/x \n')
    f.write('       </DataItem>\n')
    f.write('       <DataItem Dimensions="'+str(nx)+' '+str(ny)+' '+str(nz)+'" NumberType="Float" Precision="8" Format="HDF">\n')
    f.write('         '+mesh_name+':/y \n')
    f.write('       </DataItem>\n')
    f.write('       <DataItem Dimensions="'+str(nx)+' '+str(ny)+' '+str(nz)+'" NumberType="Float" Precision="8" Format="HDF">\n')
    f.write('         '+mesh_name+':/z \n')
    f.write('       </DataItem>\n')
    f.write('     </Geometry>\n')
    ieqn = 0
    neqn = size(var_names)
    for var_name in var_names:
        f.write('     <Attribute Name="'+var_name+'" AttributeType="Scalar" Center="Node">\n')
        f.write('       <DataItem ItemType="HyperSlab" Dimensions="1 '+str(nx)+' '+str(ny)+' '+str(nz)+'" Type="HyperSlab">\n')
        f.write('         <DataItem Dimensions="3 4" Format="XML">\n')
        f.write('           '+str(ieqn)+' '+str(0)+' '+str(0)+' '+str(0)+' \n')
        f.write('           '+str(1)+' '+str(1)+' '+str(1)+' '+str(1)+' \n')
        f.write('           '+str(1)+' '+str(nx)+' '+str(ny)+' '+str(nz)+' \n')
        f.write('         </DataItem>\n')
        f.write('         <DataItem Dimensions="'+str(neqn)+' '+str(nx)+' '+str(ny)+' '+str(nz)+'" NumberType="Float" Precision="8" Format="HDF">\n')
        f.write('           '+filename+':/patch1/q\n')
        f.write('         </DataItem>\n')
        f.write('       </DataItem>\n')
        f.write('     </Attribute>\n')
        ieqn+=1
    f.write('   </Grid>\n')
    f.write(' </Domain> \n')
    f.write('</Xdmf>')
    f.close()
