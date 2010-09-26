
# Library version of Makefile
# ---------------------------

# Only "make .objs" and "make .htmls" are useful in this directory.

# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/util/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.

# Note: This version uses dynamic memory allocation from f90 and so requires
# an f90 compiler such as gfortran.

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# Add any desired compiler flags such as -g here:
FFLAGS ?= -O2

CLAW_LIB = $(shell pwd)

CLAW_SOURCES = \
  $(CLAW_LIB)/amr2ez.f \
  $(CLAW_LIB)/bc2amr.f \
  $(CLAW_LIB)/setaux.f \
  $(CLAW_LIB)/setprob.f \
  $(CLAW_LIB)/b4step2.f \
  $(CLAW_LIB)/advanc.f \
  $(CLAW_LIB)/bound.f \
  $(CLAW_LIB)/stepgrid.f \
  $(CLAW_LIB)/estdt.f \
  $(CLAW_LIB)/errest.f \
  $(CLAW_LIB)/flag2refine.f  \
  $(CLAW_LIB)/allowflag.f  \
  $(CLAW_LIB)/spest.f  \
  $(CLAW_LIB)/igetsp.f \
  $(CLAW_LIB)/reclam.f \
  $(CLAW_LIB)/birect.f \
  $(CLAW_LIB)/cleanup.f \
  $(CLAW_LIB)/colate.f \
  $(CLAW_LIB)/bufnst.f \
  $(CLAW_LIB)/errf1.f \
  $(CLAW_LIB)/gfixup.f \
  $(CLAW_LIB)/filval.f \
  $(CLAW_LIB)/fixcapaq.f \
  $(CLAW_LIB)/filpatch.f \
  $(CLAW_LIB)/prefilp.f \
  $(CLAW_LIB)/flglvl.f \
  $(CLAW_LIB)/fluxad.f \
  $(CLAW_LIB)/fluxsv.f \
  $(CLAW_LIB)/ginit.f \
  $(CLAW_LIB)/grdfit.f \
  $(CLAW_LIB)/intfil.f \
  $(CLAW_LIB)/moment.f \
  $(CLAW_LIB)/nestck.f \
  $(CLAW_LIB)/prepc.f \
  $(CLAW_LIB)/prepf.f \
  $(CLAW_LIB)/projec.f \
  $(CLAW_LIB)/signs.f \
  $(CLAW_LIB)/findcut.f \
  $(CLAW_LIB)/smartbis.f \
  $(CLAW_LIB)/putnod.f \
  $(CLAW_LIB)/putsp.f \
  $(CLAW_LIB)/regrid.f \
  $(CLAW_LIB)/setgrd.f \
  $(CLAW_LIB)/setuse.f \
  $(CLAW_LIB)/stst1.f \
  $(CLAW_LIB)/tick.f \
  $(CLAW_LIB)/trimbd.f \
  $(CLAW_LIB)/update.f \
  $(CLAW_LIB)/nodget.f \
  $(CLAW_LIB)/upbnd.f  \
  $(CLAW_LIB)/basic.f \
  $(CLAW_LIB)/outval.f \
  $(CLAW_LIB)/copysol.f \
  $(CLAW_LIB)/outvar.f \
  $(CLAW_LIB)/outmsh.f \
  $(CLAW_LIB)/outtre.f \
  $(CLAW_LIB)/domain.f  \
  $(CLAW_LIB)/cellave.f \
  $(CLAW_LIB)/fdisc.f \
  $(CLAW_LIB)/fss.f \
  $(CLAW_LIB)/zeroin.f \
  $(CLAW_LIB)/setflags.f \
  $(CLAW_LIB)/shiftset.f \
  $(CLAW_LIB)/conck.f \
  $(CLAW_LIB)/domshrink.f \
  $(CLAW_LIB)/domprep.f \
  $(CLAW_LIB)/domup.f \
  $(CLAW_LIB)/domcopy.f \
  $(CLAW_LIB)/coarsen.f \
  $(CLAW_LIB)/auxcoarsen.f \
  $(CLAW_LIB)/intcopy.f \
  $(CLAW_LIB)/preintcopy.f \
  $(CLAW_LIB)/icall.f \
  $(CLAW_LIB)/preicall.f \
  $(CLAW_LIB)/step2.f \
  $(CLAW_LIB)/flux2.f \
  $(CLAW_LIB)/flux2fw.f \
  $(CLAW_LIB)/inlinelimiter.f \
  $(CLAW_LIB)/limiter.f \
  $(CLAW_LIB)/philim.f \
  $(CLAW_LIB)/cstore.f \
  $(CLAW_LIB)/qad.f \
  $(CLAW_LIB)/src2.f \
  $(CLAW_LIB)/src1d.f \
  $(CLAW_LIB)/saveqc.f \
  $(CLAW_LIB)/check.f \
  $(CLAW_LIB)/restrt.f \
  $(CLAW_LIB)/valout.f \
  $(CLAW_LIB)/check4nans.f \
  $(CLAW_LIB)/check_hdf.f \
  $(CLAW_LIB)/restrt_hdf.f \
  $(CLAW_LIB)/valout_hdf.f \
  $(CLAW_LIB)/opendatafile.f \
  $(CLAW_LIB)/dumpgauge.f \
  $(CLAW_LIB)/setgauges.f \
  $(CLAW_LIB)/quick_sort1.f \
  $(CLAW_LIB)/rpn2from1d.f \
  $(CLAW_LIB)/rpt2.f \
  $(CLAW_LIB)/init_alloc.f90 \
  $(CLAW_LIB)/restrt_alloc.f90 \
  $(CLAW_LIB)/resize_alloc_static.f90 \
  $(CLAW_LIB)/resize_alloc.f90 \

CLAW_COMMONS = \
  $(CLAW_LIB)/call.i \
  $(CLAW_LIB)/gauges.i \
  $(CLAW_LIB)/regions.i

default: .objs

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


### DO NOT remove this line - make depends on it ###

