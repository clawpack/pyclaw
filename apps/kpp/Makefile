PYCLAWMAKE = $(PYCLAW)/Makefile.common

# No transverse solver implemented for KPP, so acoustics is included
# here as a hack.  But this will fail for multi-D algorithm!
RP_SOURCE =   $(RIEMANN)/src/rpn2_kpp.f90 $(RIEMANN)/src/rpt2_dummy.f90

all: classic2.so sharpclaw2.so

include $(PYCLAWMAKE)

