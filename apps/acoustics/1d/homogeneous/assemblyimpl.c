#include "assemblyimpl.h"
#include <petscdmda.h>

typedef struct {
  PetscScalar p,u;
} Node;

static inline PetscScalar UpwindFlux1D(PetscReal u,PetscScalar hL,PetscScalar hR)
{return (u > 0) ? hL*u : hR*u;}

#undef  __FUNCT__
#define __FUNCT__ "FormIJacobian"
PetscErrorCode FormIJacobian(TS ts, PetscReal t, Vec X, Vec Xdot, PetscReal a, Mat J, Mat P, Params *p)
{
  PetscErrorCode ierr;
  DM da;
  DMDALocalInfo info;
  PetscReal hx,c = p->cc,z = p->zz;
  Node *x;
  Vec Xloc;
  PetscInt i;

  PetscFunctionBegin;
  ierr = TSGetDM(ts,&da);CHKERRQ(ierr);
  ierr = DMDAGetLocalInfo(da,&info);CHKERRQ(ierr);
  hx = (p->xmax - p->xmin) / info.mx;

  /* This is a linear problem so it does not use the state to assemble the Jacobian. */
  ierr = DMGetLocalVector(da,&Xloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalBegin(da,X,INSERT_VALUES,Xloc);CHKERRQ(ierr);
  ierr = DMGlobalToLocalEnd(da,X,INSERT_VALUES,Xloc);CHKERRQ(ierr);
  ierr = DMDAVecGetArray(da,Xloc,&x);CHKERRQ(ierr);

  for (i=info.xs; i<info.xs+info.xm; i++) {
    MatStencil rows[1] = {{0,0,i,0}}, cols[3] = {{0,0,i-1,0},{0,0,i,0},{0,0,i+1,0}};
    /* Flux: q_t + (A q)_x
     * z = sqrt(rho*bulk), c = sqrt(rho/bulk)
     * Spectral decomposition: A = R * D * Rinv
     * [    cz] = [-z   z] [-c    ] [ 1/2z  1/2]
     * [c/z   ] = [ 1   1] [     c] [-1/2z  1/2]
     *
     * We decompose this into the left-traveling waves Al = R * D^- Rinv
     * and the right-traveling waves Ar = R * D^+ * Rinv
     * Multiplying out these expressions produces the following two matrices
     */
    PetscReal
      Al[2][2] = {{-c/2     , c*z/2  },
                  {c/(2*z)  , -c/2   }}, /* Left traveling waves */
      Ar[2][2] = {{c/2      , c*z/2  },
                  {c/(2*z)  , c/2    }}; /* Right traveling waves */
    PetscScalar v[2][3][2];
    PetscInt j,k;
    ierr = PetscMemzero(v,sizeof v);CHKERRQ(ierr);
    for (j=0; j<2; j++) {
      for (k=0; k<2; k++) {
        v[j][0][k] += -Ar[j][k] / hx;
        v[j][1][k] += (Ar[j][k] - Al[j][k]) / hx;
        v[j][2][k] += Al[j][k] / hx;
      }
    }
    /* Add shift for time integration */
    v[0][1][0] += a;
    v[1][1][1] += a;
    ierr = MatSetValuesBlockedStencil(P,1,rows,3,cols,&v[0][0][0],INSERT_VALUES);CHKERRQ(ierr);
  }
  ierr = DMRestoreLocalVector(da,&Xloc);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(P,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(P,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
