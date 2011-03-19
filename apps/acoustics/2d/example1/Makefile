
FC ?= gfortran


TWO_D_LIBPATH = $(PETCLAW)/src/fortran/2d/nonInterleaved/
LIBOBJ = \
   step2ds.o\
   flux2.o\
   limiter.o\
   philim.o\

SOURCE =   rpn2ac.f rpt2ac.f

dimsp2.so: $(TWO_D_LIBPATH)dimsp2.f $(SOURCE) $(LIBOBJ)
	f2py -m dimsp2 -c $(TWO_D_LIBPATH)dimsp2.f $(LIBOBJ) $(SOURCE)  
        
limiter.o: 
	$(FC) -o limiter.o -c $(TWO_D_LIBPATH)limiter.f

philim.o:
	$(FC) -o philim.o -c $(TWO_D_LIBPATH)philim.f

flux2.o:
	$(FC) -o flux2.o -c $(TWO_D_LIBPATH)flux2.f

step2ds.o:
	$(FC) -o step2ds.o -c $(TWO_D_LIBPATH)step2ds.f

clean:
	rm -f *.o *.so
	
clobber:
	rm -f *.o *.so
	rm -f -R *output*
