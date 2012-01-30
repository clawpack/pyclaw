PYCLAWMAKE = $(PYCLAW)/Makefile.common

FC ?= gfortran

RP_SOURCE =   $(RIEMANN)/src/rpn2_vc_advection.f90 $(RIEMANN)/src/rpt2_vc_advection.f90

all: 
	make classic2.so
	make sharpclaw2.so

include $(PYCLAWMAKE)
