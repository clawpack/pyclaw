#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <petscts.h>

typedef struct Params {
  double cc, zz;                /* speed and entry in eigenvector */
  double xmin, xmax;            /* grid extents */
} Params;

PetscErrorCode FormIJacobian(TS ts, PetscReal t, Vec X, Vec Xdot, PetscReal, Mat J, Mat P, Params *p);

#endif /* !ASSEMBLY_H */
