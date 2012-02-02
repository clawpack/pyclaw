PYCLAWMAKE = $(PYCLAW)/Makefile.common

RP_SOURCE =   $(RIEMANN)/src/rpn2_vc_acoustics.f90 $(RIEMANN)/src/rpt2_vc_acoustics.f90

all: 
	make classic2.so
	make sharpclaw2.so

include $(PYCLAWMAKE)

