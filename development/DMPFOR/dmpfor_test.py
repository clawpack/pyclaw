from petsc4py import PETSc
import numpy as np
import DMPFOR
from six.moves import range


global_nx =3
global_ny =2
dof=4

da = PETSc.DA().create(dim=2,
dof=dof,
sizes=[global_nx, global_ny], 
#periodic_type = PETSc.DA.PeriodicType.GHOSTED_XYZ,
#stencil_type=self.STENCIL,
#stencil_width=2,
comm=PETSc.COMM_WORLD)


gVec = da.createGlobalVector()
lVec = da.createLocalVector()




ranges = da.getRanges()

nx_start = ranges[0][0]
nx_end = ranges[0][1]
ny_start = ranges[1][0]
ny_end = ranges[1][1]

nx = nx_end - nx_start
ny = ny_end - ny_start


q = np.empty((dof, nx, ny), order='F')

for i in range(0,nx):
    for j in range(0,ny):
        for k in range(0,dof):
            q[k,i,j] = k+10*i+100*j

gVec.array = q

q = gVec.array.reshape((dof, nx, ny), order='F')

print("da array from python")
print(q)


print("da array from fortran")
DMPFOR.dmpfor(q,dof,nx,ny)


print("da array from python after rolling axes using rollaxis")
rolled_q_1 = np.rollaxis(q,0,3)
rolled_q_1 = np.reshape(rolled_q_1,(nx,ny,dof),order='F')
print(rolled_q_1)
print("da array from fortran after rolling axes using rollaxis")
DMPFOR.dmpfor(rolled_q_1,nx,ny,dof)


print("da array from python after rolling axes using element by element copy")
rolled_q_2 = np.empty((nx,ny,dof),order='F')
for i in range(0,nx):
    for j in range(0,ny):
        for k in range(0,dof):
            rolled_q_2[i,j,k] = q[k,i,j]
print(rolled_q_2)
print("da array from fortran after rolling axes using element by element copy")
DMPFOR.dmpfor(rolled_q_2,nx,ny,dof)







