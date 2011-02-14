
step1.so : step1.f limiter.o philim.o rp1.o
	f2py  -m step1 -c --noopt --debug step1.f limiter.o rp1.o philim.o 
 
limiter.o:
	gfortran -o limiter.o -c limiter.f

philim.o:
	gfortran -o philim.o -c philim.f

rp1.o:
	gfortran -o rp1.o -c rp1.f

clean:
	rm -f *.o *.so
	
clobber:
	rm -f *.o *.so
	rm -f -R *output*
