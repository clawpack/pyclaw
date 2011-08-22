all : limiter.o philim.o 

limiter.o:
	$(FC) $(FFLAGS) -o limiter.o -c limiter.f95

philim.o:
	$(FC) $(FFLAGS) -o philim.o -c philim.f95

clean:
	rm *.o *.so

