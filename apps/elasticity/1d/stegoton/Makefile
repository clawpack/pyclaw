#FC=f77
FLAGS=-qextname
F2PYFLAGS=--fcompiler=bgpgnu  --compiler=mpicc

step1.so : step1fw.f limiter.o philim.o rp1nel.o
	f2py $(F2PYFLAGS) -m step1 -c step1fw.f limiter.f rp1nel.f philim.f 
 
limiter.o:
	$(FC) $(FLAGS) -o limiter.o -c limiter.f

philim.o:
	$(FC) $(FLAGS) -o philim.o -c philim.f

rp1nel.o:
	$(FC) $(FLAGS) -o rp1nel.o -c rp1nel.f

clean:
	rm *.o *.so
	
