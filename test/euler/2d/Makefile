PYCLAWMAKE = $(PYCLAW)/Makefile.common

include $(PYCLAWMAKE)

LIBOBJ = \
   step2ds.o\
   flux2.o\
   limiter.o\
   philim.o\

SOURCE =   $(RIEMANN)/src/rpn2_euler_5wave.f $(RIEMANN)/src/rpt2_euler_5wave.f

dimsp2.so: $(TWO_D_CLASSIC)/dimsp2.f $(SOURCE) $(LIBOBJ)
	f2py -m dimsp2 -c $(TWO_D_CLASSIC)/dimsp2.f $(LIBOBJ) $(SOURCE)  
 
limiter.o: 
	$(FC) -o limiter.o -c $(ONE_D_CLASSIC)/limiter.f

philim.o:
	$(FC) -o philim.o -c $(ONE_D_CLASSIC)/philim.f

flux2.o:
	$(FC) -o flux2.o -c $(TWO_D_CLASSIC)/flux2.f

step2ds.o:
	$(FC) -o step2ds.o -c $(TWO_D_CLASSIC)/step2ds.f

clean:
	rm -f *.o *.so *.log *.pyc
	
clobber: clean
	rm -f -R *_output*
	rm -f -R *_plots*
