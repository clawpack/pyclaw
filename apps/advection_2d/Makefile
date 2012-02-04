PYCLAWMAKE = $(PYCLAW)/Makefile.common

RP_SOURCE = $(RIEMANN)/src/rpt2_advection.f90 $(RIEMANN)/src/rpn2_advection.f90

all: 
	make classic2.so
	make sharpclaw2.so

include $(PYCLAWMAKE)

