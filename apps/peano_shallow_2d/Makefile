PYCLAWMAKE = $(PYCLAW)/Makefile.common

RP_SOURCE =   $(RIEMANN)/src/rpn2_shallow_roe-with-efix.f90 $(RIEMANN)/src/rpt2_shallow_roe-with-efix.f90 

all: 
	make classic2.so
	make sharpclaw2.so

include $(PYCLAWMAKE)

