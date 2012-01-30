PYCLAWMAKE = $(PYCLAW)/Makefile.common

RP_SOURCE =   $(RIEMANN)/src/rpn2_euler_5wave.f90 $(RIEMANN)/src/rpt2_euler_5wave.f90

all: 
	make classic2.so
	make sharpclaw2.so

include $(PYCLAWMAKE)

