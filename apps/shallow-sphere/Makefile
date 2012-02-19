# Common Makefile used by all problems
# ====================================
PYCLAWMAKE = $(PYCLAW)/Makefile.common

all: 
	make classic2.so
	make problem.so

# Override step2.f with a new function that contains a call to an additional
# function
# ==========================================================================
override TWO_D_CLASSIC_SOURCES = step2qcor.f90 qcor.f90 $(CLASSIC)/step2ds.f90 $(CLASSIC)/flux2.f90 $(CLASSIC)/limiter.f90 $(CLASSIC)/philim.f90

# Include common Makefile
# =======================
include $(PYCLAWMAKE)

# Problem's source Fortran files
INITIALIZE_SOURCE = mapc2p.f90 setaux.f90 qinit.f90 src2.f90

problem.so: $(INITIALIZE_SOURCE)
	$(F2PY) -m problem -c $^

