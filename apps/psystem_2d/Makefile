PYCLAWMAKE = $(PYCLAW)/Makefile.common

FC ?= gfortran

RP_SOURCE =   $(RIEMANN)/src/rpn2_psystem.f90 $(RIEMANN)/src/rpt2_psystem.f90

all: 
	make classic2.so
	make sharpclaw2.so

include $(PYCLAWMAKE)
