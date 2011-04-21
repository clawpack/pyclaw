PETCLAWMAKE = $(PETCLAW)/Makefile.common

include $(PETCLAWMAKE)

# This makefile uses the Fortran compiler and flags from PETSc. This
# means that PETSc needs to be configured with Fortran support which
# would otherwise not be necessary. However, using arbitrary environment
# variables and hoping it all works out is way more fragile than I want
# any part of.

FCOMPILE_SINGLE = ${FC} -c ${FC_FLAGS} ${FFLAGS} ${FCPPFLAGS}

LIBOBJ = \
   limiter.o \
   philim.o

STEP1_SOURCE =  ${ONE_D_LIBPATH}/step1.f rp1.f

F2PY = f2py
PYTHON = python
MODULE = assembly

all: step1.so ${MODULE}.so

${MODULE}.so: ${MODULE}.pyx ${MODULE}impl.h ${MODULE}impl.c petscdadef.h
	CC='${PCC}' CFLAGS='${PCC_FLAGS}' F90='${FC_LINKER}' FFLAGS='${FFLAGS}' LDFLAGS='${FC_LINKER_FLAGS}' ${PYTHON} setup.py -q build_ext --inplace
	${RM} -r build ${MODULE}_wrap.c

step1.so: $(STEP1_SOURCE) $(LIBOBJ)
	${F2PY} -m step1 -c $^

limiter.o: $(ONE_D_LIBPATH)/limiter.f
	$(FCOMPILE_SINGLE) -c -o $@ $^

philim.o: $(ONE_D_LIBPATH)/philim.f
	$(FCOMPILE_SINGLE) -c -o $@ $^

clean:
	rm -f *.o *.so

clobber:
	rm -f *.o *.so
	rm -f -R *output*

include ${PETSC_DIR}/conf/variables

