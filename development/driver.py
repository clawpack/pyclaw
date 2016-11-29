#!/usr/bin/env python
from __future__ import absolute_import
from __future__ import print_function
from petsc4py import PETSc
from six.moves import range

class PetCLAW:
  def advection1D(self, N):
    '''David: If you put in the linear algebra that you need (comments), I will move it to code'''
    da = PETSc.DA().create([N])
    da.view()
    f = da.createGlobalVector()
    f.view()
    a = f.getArray()
    for i in range(f.getSize()):
      a[i] = i*i
    f.view()
    return

  def run(self):
    self.advection1D(10)
    print('Done')
    return

if __name__ == '__main__':
  PetCLAW().run()
