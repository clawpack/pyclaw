

LIBOBJ = \
   $(PETCLAW)/src/fortran/1d/limiter.o \
   $(PETCLAW)/src/fortran/1d/philim.o

OBJ =   rp1.f

step1.so: $(PETCLAW)/src/fortran/1d/step1.f $(LIBOBJ) $(OBJ)
	f2py -m step1 -c $(PETCLAW)/src/fortran/1d/step1.f $(LIBOBJ) $(OBJ) 
        
rp1.o:
	gfortran -o rp1.o -c rp1.f

clean:
	rm -f *.o *.so
	
clobber:
	rm -f *.o *.so
	rm -f -R *output*
