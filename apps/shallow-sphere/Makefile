# Common Makefile used by all problems
PYCLAWMAKE = $(PYCLAW)/Makefile.common

# Source file used by the Makefile.common
RP_SOURCE =   $(RIEMANN)/src/rpn2swsphere.f $(RIEMANN)/src/rpt2swq.f

all: classic2.so sharpclaw2.so 

include $(PYCLAWMAKE)


# The following Makefile called, Makefile.add, contains the
# ADDitional fortran routine that are use for this specific test 
# case.
ADDMAKE = ./Makefile.add

# Source file used to initialize q and aux
INITIALIZE_SOURCE = mapc2m.f setaux.f qinit.f

all: init.so

include $(ADDMAKE)


