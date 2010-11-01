

step1.so : step1.f limiter.o philim.o rp1.o
	f2py -m step1 -c  step1.f limiter.o rp1.o philim.o 
	
limiter.o:
	$(FC) -o limiter.o -c limiter.f

philim.o:
	$(FC) -o philim.o -c philim.f

rp1.o:
	$(FC) -o rp1.o -c rp1.f

clean:
	rm *.o *.so
	
