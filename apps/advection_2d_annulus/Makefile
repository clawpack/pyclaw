PYCLAWMAKE = $(PYCLAW)/Makefile.common

FC ?= gfortran

RP_SOURCE =   $(RIEMANN)/src/rpn2_vc_advection.f90 $(RIEMANN)/src/rpt2_vc_advection.f90

all: classic2.so sharpclaw2.so

include $(PYCLAWMAKE)
