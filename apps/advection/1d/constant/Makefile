PETCLAWMAKE = $(PETCLAW)/Makefile.common

include $(PETCLAWMAKE)

ONE_D_LIBPATH = $(PETCLAW)/src/fortran/1d

LIBOBJ = \
    limiter.o\
    philim.o\

SOURCE = rp1.f

step1.so: $(ONE_D_LIBPATH)/step1.f $(SOURCE) $(LIBOBJ) 
	f2py  -m step1 -c $(ONE_D_LIBPATH)/step1.f $(SOURCE) $(LIBOBJ)
 
limiter.o:
	$(FC) -o limiter.o -c $(ONE_D_LIBPATH)/limiter.f

philim.o:
	$(FC) -o philim.o -c $(ONE_D_LIBPATH)/philim.f


clean:
	rm -f *.o *.so
	
clobber:
	rm -f *.o *.so
	rm -f -R *output*


