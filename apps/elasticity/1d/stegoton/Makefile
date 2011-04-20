PETCLAWMAKE = $(PETCLAW)/Makefile.common

include $(PETCLAWMAKE)

#FC=f77
#FLAGS=-qextname
#F2PYFLAGS=--fcompiler=bgpgnu  --compiler=mpicc

LIBOBJ = \
    limiter.o\
    philim.o\

SOURCE = rp1nel.f

step1.so : $(ONE_D_LIBPATH)/step1fw.f $(SOURCE) $(LIBOBJ)
	f2py $(F2PYFLAGS) -m step1 -c $(ONE_D_LIBPATH)/step1fw.f $(SOURCE) $(LIBOBJ) 
 
limiter.o:
	$(FC) $(FLAGS) -o limiter.o -c $(ONE_D_LIBPATH)/limiter.f

philim.o:
	$(FC) $(FLAGS) -o philim.o -c $(ONE_D_LIBPATH)/philim.f

clean:
	rm *.o *.so
	
