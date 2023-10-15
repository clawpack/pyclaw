F2PY = f2py

CLASSIC = ../../src/pyclaw/classic

TWO_D_CLASSIC_SOURCES = step2qcor.f90 qcor.f90 $(CLASSIC)/step2ds.f90 $(CLASSIC)/flux2.f90 $(CLASSIC)/limiter.f90 $(CLASSIC)/philim.f90

all: 
	make classic2.so
	make problem.so

classic2.so: $(TWO_D_CLASSIC_SOURCES) $(RP_SOURCE) 
	${F2PY} -m classic2 -c $(TWO_D_CLASSIC_SOURCES) $(RP_SOURCE) 	

problem.so: mapc2p.f90 setaux.f90 qinit.f90 src2.f90
	$(F2PY) -m problem -c $^

clean:
	rm -f *.o *.so *.pyc *.log

clobber: clean
	rm -rf _output/
	rm -rf _plots/
