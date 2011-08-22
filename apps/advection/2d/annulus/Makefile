PYCLAWMAKE = $(PYCLAW)/Makefile.common

FC ?= gfortran

RP_SOURCE =   $(RIEMANN)/src/rpn2_vc_advection.f95 $(RIEMANN)/src/rpt2_vc_advection.f95

all: classic2.so sharpclaw2.so

include $(PYCLAWMAKE)
