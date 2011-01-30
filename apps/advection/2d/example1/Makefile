


objects = step2ds.o  flux2.o  limiter.o  philim.o rpt2.f rpn2.f
step1.so: dimsp2.f $(objects)
	f2py -m dimsp2 -c dimsp2.f $(objects)  
        
limiter.o: limiter.f
	gfortran -o limiter.o -c limiter.f

philim.o:
	gfortran -o philim.o -c philim.f

flux2.o:
	gfortran -o flux2.o -c flux2.f

step2ds.o:
	gfortran -o step2ds.o -c step2ds.f

rpt2.o:
	gfortran -o rpt2.o -c rpt2.f

rpn2.o:
	gfortran -o rpn2.o -c rpn2.f

clean:
	rm -f *.o *.so
	
clobber:
	rm -f *.o *.so
	rm -f -R *output*
