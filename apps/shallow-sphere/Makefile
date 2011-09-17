# Common Makefile used by all problems
######################################
PYCLAWMAKE = $(PYCLAW)/Makefile.common

# Source file used by the Makefile.common
#########################################
RP_SOURCE =   $(RIEMANN)/src/rpn2swsphere.f $(RIEMANN)/src/rpt2swq.f

all: classic2.so sharpclaw2.so init.so


# Source files and commands used to initialize q and aux
########################################################
INITIALIZE_SOURCE = mapc2m.f qinit.f

init.so: $(INITIALIZE_SOURCE)
	$(F2PY) -m init -c $^

include $(PYCLAWMAKE)
