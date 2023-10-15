F2PY = f2py

SHARPCLAW = ../../src/pyclaw/sharpclaw

ONE_D_SHARPCLAW_SOURCES = evec.f90 $(SHARPCLAW)/flux1.f90 $(SHARPCLAW)/ClawParams.f90 $(SHARPCLAW)/workspace.f90 $(SHARPCLAW)/weno.f90 $(SHARPCLAW)/reconstruct.f90

all: 
	make sharpclaw1.so
	make euler_tfluct.so

sharpclaw1.so: $(ONE_D_SHARPCLAW_SOURCES) 
	${F2PY} -m sharpclaw1 -c $(ONE_D_SHARPCLAW_SOURCES) 

euler_tfluct.so: euler_tfluct.f90
	f2py -m $(basename $(notdir $@)) -c $^

clean:
	rm -f *.o *.so *.pyc *.log

clobber: clean
	rm -rf _output/
	rm -rf _plots/
