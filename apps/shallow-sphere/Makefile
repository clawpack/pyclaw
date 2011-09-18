# Common Makefile used by all problems
######################################
PYCLAWMAKE = $(PYCLAW)/Makefile.common

# Source file used by the Makefile.common
#########################################
RP_SOURCE =   $(RIEMANN)/src/rpn2swsphere.f $(RIEMANN)/src/rpt2swq.f

TWO_D_CLASSIC_SOURCES = step2qcor.f flux2.o limiter.o philim.o

all: classic2.so sharpclaw2.so problem.so


# Source files and commands used to initialize q and aux and define
# the source routine function.
########################################################
INITIALIZE_SOURCE = mapc2m.f setaux.f qinit.f src2.f

problem.so: $(INITIALIZE_SOURCE)
	$(F2PY) -m problem -c $^

include $(PYCLAWMAKE)
