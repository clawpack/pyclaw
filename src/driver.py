#!/usr/bin/env python
from petsc4py import PETSc

class PetCLAW:
  def advection1D(self, N):
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
    print 'Done'
    return

if __name__ == '__main__':
  PetCLAW().run()
