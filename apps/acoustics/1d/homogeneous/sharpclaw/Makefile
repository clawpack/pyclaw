PETCLAWMAKE = $(PETCLAW)/Makefile.common

include $(PETCLAWMAKE)

LIBOBJ = reconstruct.o

SOURCE =   ../rp1.f $(ONE_D_LIBPATH)/sharpclaw/reconstruct.f90

flux1.so: $(ONE_D_LIBPATH)/sharpclaw/flux1.f90 $(LIBOBJ) $(SOURCE)
	f2py -m flux1 -c $(ONE_D_LIBPATH)/sharpclaw/flux1.f90 $(SOURCE) 

reconstruct.o:
	$(FC) -o reconstruct.o -c $(ONE_D_LIBPATH)/sharpclaw/reconstruct.f90

ClawParams.o:
	$(FC) -o ClawParams.o -c $(ONE_D_LIBPATH)/sharpclaw/ClawParams.f90

clean:
	rm -f *.o *.so *.log *.mod *.pyc

clobber: clean
	rm -f -R *output*
