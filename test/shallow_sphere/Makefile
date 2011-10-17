# Common Makefile used by all problems
# ====================================
PYCLAWMAKE = $(PYCLAW)/Makefile.common

# Source file used by the Makefile.common
# =======================================
RP_SOURCE =   $(RIEMANN)/src/rpn2swsphere.f $(RIEMANN)/src/rpt2swq.f

all: classic2.so problem.so


# Override step2.f with a new function that contains a call to an additional
# function, i.e. qcor.f
# ==========================================================================
override TWO_D_CLASSIC_SOURCES = step2qcor.f qcor.o flux2.o limiter.o philim.o

qcor.o: qcor.f
	$(FC) $(FFLAGS) -o qcor.o -c qcor.f


# Include common Makefile
# =======================
include $(PYCLAWMAKE)

# Problem's source Fortran files
INITIALIZE_SOURCE = mapc2p.f setaux.f qinit.f src2.f

problem.so: $(INITIALIZE_SOURCE)
	$(F2PY) -m problem -c $^

