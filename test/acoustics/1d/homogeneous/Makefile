PYCLAWMAKE = $(PYCLAW)/Makefile.common

# This makefile uses the Fortran compiler and flags from PETSc. This
# means that PETSc needs to be configured with Fortran support which
# would otherwise not be necessary. However, using arbitrary environment
# variables and hoping it all works out is way more fragile than I want
# any part of.

FCOMPILE_SINGLE = ${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS}

RP_SOURCE =  $(RIEMANN)/src/rp1_acoustics.f

F2PY = f2py
PYTHON = python
MODULE = assembly

all: classic1.so sharpclaw1.so ${MODULE}.so

include ${PETSC_DIR}/conf/variables

include $(PYCLAWMAKE)

