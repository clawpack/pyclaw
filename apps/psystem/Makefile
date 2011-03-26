
FC ?= gfortran


TWO_D_LIBPATH = $(PETCLAW)/src/fortran/2d/interleaved/
LIBOBJ = \
   step2ds.o\
   flux2.o\
   limiter.o\
   philim.o\

SOURCE =   rpn2pv.f rpt2pv.f

dimsp2.so: $(TWO_D_LIBPATH)dimsp2.pyf $(SOURCE) $(LIBOBJ)
	f2py -c $(TWO_D_LIBPATH)dimsp2.pyf $(TWO_D_LIBPATH)dimsp2.f $(LIBOBJ) $(SOURCE)  

#dimsp2.pyf: dimsp2.f $(SOURCE) 
#	f2py -m dimsp2 -h dimsp2.pyf  dimsp2.f $(SOURCE)
        
#dimsp2.so: dimsp2.pyf $(SOURCE) $(LIBOBJ)
#	f2py -c dimsp2.pyf dimsp2.f $(LIBOBJ) $(SOURCE)  

 
limiter.o: 
	$(FC) -o limiter.o -c $(TWO_D_LIBPATH)limiter.f

philim.o:
	$(FC) -o philim.o -c $(TWO_D_LIBPATH)philim.f

flux2.o:
	$(FC) -o flux2.o -c $(TWO_D_LIBPATH)flux2fw.f

step2ds.o:
	$(FC) -o step2ds.o -c $(TWO_D_LIBPATH)step2ds.f

clean:
	rm -f *.o *.so
	
clobber:
	rm -f *.o *.so
	rm -f -R *output*
