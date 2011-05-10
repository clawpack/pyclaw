PYCLAWMAKE = $(PYCLAW)/Makefile.common

include $(PYCLAWMAKE)

# This makefile uses the Fortran compiler and flags from PETSc. This
# means that PETSc needs to be configured with Fortran support which
# would otherwise not be necessary. However, using arbitrary environment
# variables and hoping it all works out is way more fragile than I want
# any part of.

FCOMPILE_SINGLE = ${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS}

LIBOBJ = \
   limiter.o \
   philim.o

STEP1_SOURCE =  ${ONE_D_CLASSIC}/step1.f $(RIEMANN)/src/rp1_acoustics.f

F2PY = f2py
PYTHON = python
MODULE = assembly

all: step1.so flux1.so ${MODULE}.so

step1.so: $(STEP1_SOURCE) $(LIBOBJ)
	${F2PY} -m step1 -c $^
	
flux1.so: $(ONE_D_SHARPCLAW)/flux1.f90 $(SOURCE)
	f2py -m flux1 -c $(ONE_D_SHARPCLAW)/flux1.f90 $(RIEMANN)/src/rp1_acoustics.f $(ONE_D_SHARPCLAW)/reconstruct.f90 $(ONE_D_SHARPCLAW)/workspace.f90

${MODULE}.so: ${MODULE}.pyx ${MODULE}impl.h ${MODULE}impl.c petscdadef.h
	CC='${PCC}' CFLAGS='${PCC_FLAGS}' F90='${FC_LINKER}' FFLAGS='${FFLAGS}' LDFLAGS='${FC_LINKER_FLAGS}' ${PYTHON} setup.py -q build_ext --inplace
	${RM} -r build ${MODULE}_wrap.c

limiter.o: $(ONE_D_CLASSIC)/limiter.f
	$(FCOMPILE_SINGLE) -c -o $@ $^

philim.o: $(ONE_D_CLASSIC)/philim.f
	$(FCOMPILE_SINGLE) -c -o $@ $^

clean:
	rm -f *.o *.so *.pyc

clobber: clean
	rm -f -r _output/
	rm -f -r _plots

include ${PETSC_DIR}/conf/variables

