from __future__ import absolute_import
from __future__ import print_function
from numpy import *
from FRS import vec_rp
q = random.random((2,800))
aux = random.random((2,800))
waves, s = vec_rp(2,q[:,0:799],q[:,1:800],aux[:,0:799],aux[:,1:800])
print(waves)
print(s)



