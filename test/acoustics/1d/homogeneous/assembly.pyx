from petsc4py.PETSc cimport Vec,  PetscVec
from petsc4py.PETSc cimport Mat,  PetscMat
from petsc4py.PETSc cimport TS,   PetscTS

from petsc4py.PETSc import Error

cdef extern from "petsccompat.h":
    pass

cdef extern from "assemblyimpl.h":
    ctypedef struct Params:
        double cc, zz, xmin, xmax
    int FormIJacobian (PetscTS ts, double t, PetscVec X, PetscVec Xdot, double a, PetscMat J, PetscMat P, Params *p)

cdef extern from "petscdadef.h":
    pass

def formIJacobian(TS ts, double t, Vec X, Vec Xdot, double a, Mat J, Mat P, double cc, double zz, double xmin, double xmax):
    cdef int ierr
    cc,zz,xmin,xmax = map(float,(cc,zz,xmin,xmax))
    cdef Params p = dict(cc=cc,zz=zz,xmin=xmin,xmax=xmax)
    ierr = FormIJacobian(ts.ts, t, X.vec, Xdot.vec, a, J.mat, P.mat, &p)
    if ierr != 0: raise Error(ierr)
    if J != P: J.assemble() # for matrix-free operator
    return Mat.Structure.SAME_NONZERO_PATTERN
