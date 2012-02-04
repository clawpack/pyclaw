PYCLAWMAKE = $(PYCLAW)/Makefile.common

# No transverse solver implemented for KPP, so acoustics is included
# here as a hack.  But this will fail for multi-D algorithm!
RP_SOURCE =   $(RIEMANN)/src/rpn2_kpp.f90 $(RIEMANN)/src/rpt2_dummy.f90

all: 
	make classic2.so
	make sharpclaw2.so

include $(PYCLAWMAKE)

