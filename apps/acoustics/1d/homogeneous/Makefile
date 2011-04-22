PETCLAWMAKE = $(PETCLAW)/Makefile.common

include $(PETCLAWMAKE)

LIBOBJ = \
   limiter.o \
   philim.o

SOURCE =   rp1.f

step1.so: $(ONE_D_LIBPATH)/step1.f $(LIBOBJ) $(SOURCE)
	echo $(ONE_D_LIBPATH)
	f2py -m step1 -c $(ONE_D_LIBPATH)/step1.f $(SOURCE) $(LIBOBJ) 

limiter.o:
	$(FC) -o limiter.o -c $(ONE_D_LIBPATH)/limiter.f

philim.o:
	$(FC) -o philim.o -c $(ONE_D_LIBPATH)/philim.f

clean:
	rm -f *.o *.so *.pyc

clobber: clean
	rm -f -r _output/
	rm -f -r _plots

