PETCLAWMAKE = $(PETCLAW)/Makefile.common

include $(PETCLAWMAKE)

LIBOBJ = \
   step2ds.o\
   flux2.o\
   limiter.o\
   philim.o\

SOURCE =   rpn2sw.f rpt2sw.f

dimsp2.so: $(TWO_D_CLASSIC)/dimsp2.f $(SOURCE) $(LIBOBJ)
	f2py -m dimsp2 -c $(TWO_D_CLASSIC)/dimsp2.f $(LIBOBJ) $(SOURCE)  
 
#dimsp2.so: $(TWO_D_CLASSIC)/dimsp2.pyf $(SOURCE) $(LIBOBJ)
#	f2py -c $(TWO_D_CLASSIC)/dimsp2.pyf $(TWO_D_CLASSIC)/dimsp2.f $(LIBOBJ) $(SOURCE)  

#dimsp2.pyf: dimsp2.f $(SOURCE) 
#	f2py -m dimsp2 -h dimsp2.pyf  dimsp2.f $(SOURCE)
        
#dimsp2.so: dimsp2.pyf $(SOURCE) $(LIBOBJ)
#	f2py -c dimsp2.pyf dimsp2.f $(LIBOBJ) $(SOURCE)  

 
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

