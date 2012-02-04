PYCLAWMAKE = $(PYCLAW)/Makefile.common

RP_SOURCE =   $(RIEMANN)/src/rpn3_vc_acoustics.f90 $(RIEMANN)/src/rpt3_vc_acoustics.f90 $(RIEMANN)/src/rptt3_vc_acoustics.f90

all: 
	make classic3.so
	make sharpclaw2.so

include $(PYCLAWMAKE)

