PYCLAWMAKE = $(PYCLAW)/Makefile.common

include $(PYCLAWMAKE)

#FC=f77
#FLAGS=-qextname
#F2PYFLAGS=--fcompiler=bgpgnu  --compiler=mpicc

LIBOBJ = \
    limiter.o\
    philim.o\

SOURCE = $(RIEMANN)/src/rp1_nonlinear_elasticity.f

all: step1.so

step1.so : $(ONE_D_CLASSIC)/step1fw.f $(SOURCE) $(LIBOBJ)
	f2py $(F2PYFLAGS) -m step1 -c $(ONE_D_CLASSIC)/step1fw.f $(SOURCE) $(LIBOBJ) 
 
limiter.o:
	$(FC) $(FLAGS) -o limiter.o -c $(ONE_D_CLASSIC)/limiter.f

philim.o:
	$(FC) $(FLAGS) -o philim.o -c $(ONE_D_CLASSIC)/philim.f

clean:
	rm *.o *.so
	
