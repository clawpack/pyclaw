PYCLAWMAKE = $(PYCLAW)/Makefile.common
ADDMAKE = ./Makefile.add

RP_SOURCE =   $(RIEMANN)/src/rpn2swsphere.f $(RIEMANN)/src/rpt2swq.f

all: classic2.so sharpclaw2.so qinit.so

include $(PYCLAWMAKE)
include $(ADDMAKE)
