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

SOURCE =  $(RIEMANN)/src/rp1_acoustics.f

F2PY = f2py

step1.so: ${ONE_D_CLASSIC}/step1.f $(LIBOBJ) $(SOURCE)
	${F2PY} -m step1 -c $^

flux1.so: $(ONE_D_SHARPCLAW)/flux1.f90 reconstruct.o $(SOURCE)
	f2py -m flux1 -c $(ONE_D_SHARPCLAW)/flux1.f90 $(RIEMANN)/src/rp1_acoustics.f $(ONE_D_SHARPCLAW)/reconstruct.f90 $(ONE_D_SHARPCLAW)/workspace.f90

limiter.o: $(ONE_D_CLASSIC)/limiter.f
	$(FCOMPILE_SINGLE) -c -o $@ $^

philim.o: $(ONE_D_CLASSIC)/philim.f
	$(FCOMPILE_SINGLE) -c -o $@ $^

reconstruct.o:
	$(FC) -o reconstruct.o -c $(ONE_D_SHARPCLAW)/reconstruct.f90

clean:
	rm -f *.o *.so *.pyc

clobber: clean
	rm -f -r _output/
	rm -f -r _plots

include ${PETSC_DIR}/conf/variables

