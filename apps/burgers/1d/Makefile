PYCLAWMAKE = $(PYCLAW)/Makefile.common

include $(PYCLAWMAKE)

LIBOBJ = \
    limiter.o\
    philim.o\

SOURCE = $(RIEMANN)/src/rp1_burgers.f90

step1.so: $(ONE_D_CLASSIC)/step1.f $(SOURCE) $(LIBOBJ) 
	f2py  -m step1 -c $(ONE_D_CLASSIC)/step1.f $(SOURCE) $(LIBOBJ)
 
limiter.o:
	$(FC) -o limiter.o -c $(ONE_D_CLASSIC)/limiter.f

philim.o:
	$(FC) -o philim.o -c $(ONE_D_CLASSIC)/philim.f


clean:
	rm -f *.o *.so
	
clobber:
	rm -f *.o *.so
	rm -f -R *output*
