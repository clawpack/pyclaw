# Common Makefile used by all problems
PYCLAWMAKE = $(PYCLAW)/Makefile.common

# The following Makefile called, Makefile.add, contains the
# additional fortran routine that are use for this specific test 
# case.
ADDMAKE = ./Makefile.add

RP_SOURCE =   $(RIEMANN)/src/rpn2swsphere.f $(RIEMANN)/src/rpt2swq.f

all: classic2.so sharpclaw2.so init.so

include $(PYCLAWMAKE)
include $(ADDMAKE)
