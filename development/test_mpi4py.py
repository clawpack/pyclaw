from __future__ import absolute_import
from __future__ import print_function
from mpi4py import MPI  
comm = MPI.COMM_WORLD 
size = comm.Get_size()
rank = comm.Get_rank()
x = rank
print('x before', x)
max_x =comm.reduce( sendobj=x, op=MPI.MAX,  root=0)
x = comm.bcast(max_x, root=0)
print('x after', x)

