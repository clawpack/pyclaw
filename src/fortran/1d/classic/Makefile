all : limiter.o philim.o 
	
limiter.o:
	$(FC) $(FFLAGS) -o limiter.o -c limiter.f

philim.o:
	$(FC) $(FFLAGS) -o philim.o -c philim.f

clean:
	rm *.o *.so
	
