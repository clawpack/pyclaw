PYCLAWMAKE = $(PYCLAW)/Makefile.common

include $(PYCLAWMAKE)

LIBOBJ = reconstruct.o

SOURCE =   $(RIEMANN)/src/rp1_acoustics.f $(ONE_D_SHARPCLAW)/reconstruct.f90

flux1.so: $(ONE_D_SHARPCLAW)/flux1.f90 $(LIBOBJ) $(SOURCE)
	f2py -m flux1 -c $(ONE_D_SHARPCLAW)/flux1.f90 $(SOURCE) $(ONE_D_SHARPCLAW)/workspace.f90

reconstruct.o:
	$(FC) -o reconstruct.o -c $(ONE_D_SHARPCLAW)/reconstruct.f90

ClawParams.o:
	$(FC) -o ClawParams.o -c $(ONE_D_SHARPCLAW)/ClawParams.f90

clean:
	rm -f *.o *.so *.log *.mod *.pyc

clobber: clean
	rm -f -R *output*
