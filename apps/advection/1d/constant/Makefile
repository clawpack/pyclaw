FLAGS=-qextname
F2PYFLAGS=--fcompiler=bgpgnu  --compiler=mpicc 

step1.so : step1.f limiter.o philim.o rp1.o
	f2py $(F2PYFLAGS) -m step1 -c --noopt --debug step1.f limiter.o rp1.o philim.o 
 
limiter.o:
	$(F77) $(FLAGS) -o limiter.o -c limiter.f

philim.o:
	$(F77) $(FLAGS) -o philim.o -c philim.f

rp1.o:
	$(F77) $(FLAGS) -o rp1.o -c rp1.f

clean:
	rm -f *.o *.so
	
clobber:
	rm -f *.o *.so
	rm -f -R *output*
