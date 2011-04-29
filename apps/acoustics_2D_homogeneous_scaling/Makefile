PETCLAWMAKE = $(PETCLAW)/Makefile.common

include $(PETCLAWMAKE)

FC = mpif77 -g -fPIC -Wall -funderscoring -fno-second-underscore -fPIC -O3 -funroll-loops
F2PYFLAGS = --fcompiler=bgpgnu

LIBOBJ = \
   step2ds.o\
   flux2.o\
   limiter.o\
   philim.o\

SOURCE =   $(RIEMANN)/src/rpn2_acoustics.f $(RIEMANN)/src/rpt2_acoustics.f

dimsp2.so: $(TWO_D_CLASSIC)/dimsp2.f $(SOURCE) $(LIBOBJ)
	f2py $(F2PYFLAGS) -m dimsp2 -c $(TWO_D_CLASSIC)/dimsp2.f $(LIBOBJ) $(SOURCE)  
 
limiter.o: 
	$(FC) -o limiter.o -c $(ONE_D_CLASSIC)/limiter.f

philim.o:
	$(FC) -o philim.o -c $(ONE_D_CLASSIC)/philim.f

flux2.o:
	$(FC) -o flux2.o -c $(TWO_D_CLASSIC)/flux2.f

step2ds.o:
	$(FC) -o step2ds.o -c $(TWO_D_CLASSIC)/step2ds.f

clean:
	rm -f *.o *.so *.log
	
clobber: clean
	rm -f -R *_output*
	rm -f -R *_plots*

