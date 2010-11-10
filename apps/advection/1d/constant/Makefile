

step1.so : step1.f limiter.o philim.o rp1.o
	f2py --fcompiler=bgpgnu  --compiler=mpicc -m step1 -c --noopt --debug step1.f limiter.o rp1.o philim.o 
 
limiter.o:
	$(F77) -O0 -g -qextname -o limiter.o -c limiter.f

philim.o:
	$(F77) -O0 -g -qextname -o philim.o -c philim.f

rp1.o:
	$(F77) -O0 -g -qextname -o rp1.o -c rp1.f

clean:
	rm *.o *.so
	
