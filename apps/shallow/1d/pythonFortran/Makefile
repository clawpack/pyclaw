PETCLAWMAKE = $(PETCLAW)/Makefile.common

include $(PETCLAWMAKE)

LIBOBJ = \
   limiter.o \
   philim.o

SOURCE =   $(RIEMANN)/src/rp1_shallow_roe-with-efix.f

step1.so: $(ONE_D_CLASSIC)/step1.f $(LIBOBJ) $(SOURCE)
	f2py -m step1 -c $(ONE_D_CLASSIC)/step1.f $(SOURCE) $(LIBOBJ) 

limiter.o:
	$(FC) -o limiter.o -c $(ONE_D_CLASSIC)/limiter.f

philim.o:
	$(FC) -o philim.o -c $(ONE_D_CLASSIC)/philim.f

clean:
	rm -f *.o *.so

clobber:
	rm -f *.o *.so
	rm -f -R *output*

