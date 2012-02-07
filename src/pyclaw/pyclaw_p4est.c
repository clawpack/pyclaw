
#include <p4est.h>
#include <p8est.h>

/* TODO: This should be called elsewhere in a more sensible place */
void
pyclaw_MPI_Init (void)
{
  int argc;
  char **argv;

  argc = 0;
  argv = NULL;

  MPI_Init (&argc, &argv);
}

/* TODO: This should be called elsewhere in a more sensible place */
void
pyclaw_MPI_Finalize (void)
{
  MPI_Finalize ();
}

/* 2D p4est routines */

p4est_t *
pyclaw_p4est_new (void)
{
  p4est_connectivity_t * conn;
  p4est_t              * p4est;

  conn = p4est_connectivity_new_unitsquare ();
  p4est = p4est_new (MPI_COMM_WORLD, conn, 0, NULL, NULL);

  /* Contrary to common p4est usage we rely on conn kept alive in p4est */
  return p4est;
}

void
pyclaw_p4est_destroy (p4est_t * p4est)
{
  p4est_connectivity_t * conn;

  conn = p4est->connectivity;
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);
}

/* 3D p4est routines */

p8est_t *
pyclaw_p8est_new (void)
{
  p8est_connectivity_t * conn;
  p8est_t              * p8est;

  conn = p8est_connectivity_new_unitcube ();
  p8est = p8est_new (MPI_COMM_WORLD, conn, 0, NULL, NULL);

  /* Contrary to common p4est usage we rely on conn kept alive in p4est */
  return p8est;
}

void
pyclaw_p8est_destroy (p8est_t * p8est)
{
  p8est_connectivity_t * conn;

  conn = p8est->connectivity;
  p8est_destroy (p8est);
  p8est_connectivity_destroy (conn);
}
