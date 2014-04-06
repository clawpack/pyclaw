F2PY = f2py

SHARPCLAW = ../../src/pyclaw/sharpclaw
RIEMANN = ../../../riemann/src

ONE_D_SHARPCLAW_SOURCES = evec.f90 tfluct.f90 $(SHARPCLAW)/flux1.f90 $(SHARPCLAW)/ClawParams.f90 $(SHARPCLAW)/workspace.f90 $(SHARPCLAW)/weno.f90 $(SHARPCLAW)/reconstruct.f90

ONE_D_RIEMANN_SOURCES = $(RIEMANN)/rp1_euler_with_efix.f90

all: 
	make sharpclaw1.so
	make euler_with_efix_1D.so

sharpclaw1.so: $(ONE_D_SHARPCLAW_SOURCES) 
	${F2PY} -m sharpclaw1 -c $(ONE_D_SHARPCLAW_SOURCES) 

euler_with_efix_1D.so: $(ONE_D_RIEMANN_SOURCES)
	${F2PY} -m euler_with_efix_1D -c $(ONE_D_RIEMANN_SOURCES)

clean:
	rm -f *.o *.so *.pyc *.log

clobber: clean
	rm -rf _output/
	rm -rf _plots/
