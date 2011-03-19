
FC ?= gfortran



objects = step2ds.o  flux2.o  limiter.o  philim.o rpn2ac.f rpt2ac.f
step1.so: dimsp2.f $(objects)
	f2py -m dimsp2 -c dimsp2.f $(objects)  
        
limiter.o: limiter.f
	$(FC) -o limiter.o -c limiter.f

philim.o:
	$(FC) -o philim.o -c philim.f

flux2.o:
	$(FC) -o flux2.o -c flux2.f

step2ds.o:
	$(FC) -o step2ds.o -c step2ds.f

rpt2.o:
	$(FC) -o rpt2ac.o -c rpt2ac.f

rpn2.o:
	$(FC) -o rpn2ac.o -c rpn2ac.f

clean:
	rm -f *.o *.so
	
clobber:
	rm -f *.o *.so
	rm -f -R *output*
