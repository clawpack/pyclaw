
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/util/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.


# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = amrclaw                  # Clawpack package to use
CLAW_EXE = xamr                     # Executable to create
CLAW_setrun_file = setrun.py        # File containing function to make data
CLAW_OUTDIR = _output               # Directory for output
CLAW_setplot_file = setplot.py      # File containing function to set plots
CLAW_PLOTDIR = _plots               # Directory for plots

# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# Add any desired compiler flags such as -g here:
FFLAGS = 


# ---------------------------------
# List of sources for this program:
# ---------------------------------

    
CLAW_SOURCES = \
  ../qinit.f \
  setprob.f \
  ../rpn2.f \
  ../rpt2.f


# Clawpack library to be used:
CLAW_LIB = $(CLAW)/amrclaw/2d/lib

# Include files (containing common blocks):
CLAW_INCLUDE = $(CLAW_LIB)/call.i

    
CLAW_LIBSOURCES = \
  $(CLAW_LIB)/amr2ez.f \
  $(CLAW_LIB)/bc2amr.f \
  $(CLAW_LIB)/setaux.f \
  $(CLAW_LIB)/b4step2.f \
  $(CLAW_LIB)/qad.f \
  $(CLAW_LIB)/src2.f \
  $(CLAW_LIB)/src1d.f \
  $(CLAW_LIB)/advanc.f \
  $(CLAW_LIB)/bound.f \
  $(CLAW_LIB)/stepgrid.f \
  $(CLAW_LIB)/auxcoarsen.f \
  $(CLAW_LIB)/fixcapaq.f \
  $(CLAW_LIB)/estdt.f \
  $(CLAW_LIB)/igetsp.f \
  $(CLAW_LIB)/reclam.f \
  $(CLAW_LIB)/birect.f \
  $(CLAW_LIB)/cleanup.f \
  $(CLAW_LIB)/colate.f \
  $(CLAW_LIB)/bufnst.f \
  $(CLAW_LIB)/spest.f \
  $(CLAW_LIB)/flag2refine.f \
  $(CLAW_LIB)/allowflag.f \
  $(CLAW_LIB)/errest.f \
  $(CLAW_LIB)/errf1.f \
  $(CLAW_LIB)/gfixup.f \
  $(CLAW_LIB)/filval.f \
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
  $(CLAW_LIB)/prepf.f \
  $(CLAW_LIB)/prepc.f \
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
  $(CLAW_LIB)/setflags.f \
  $(CLAW_LIB)/shiftset.f \
  $(CLAW_LIB)/conck.f \
  $(CLAW_LIB)/domshrink.f \
  $(CLAW_LIB)/domprep.f \
  $(CLAW_LIB)/domup.f \
  $(CLAW_LIB)/domcopy.f \
  $(CLAW_LIB)/coarsen.f \
  $(CLAW_LIB)/intcopy.f \
  $(CLAW_LIB)/preintcopy.f \
  $(CLAW_LIB)/icall.f \
  $(CLAW_LIB)/preicall.f \
  $(CLAW_LIB)/step2.f \
  $(CLAW_LIB)/flux2.f \
  $(CLAW_LIB)/limiter.f \
  $(CLAW_LIB)/philim.f \
  $(CLAW_LIB)/cstore.f \
  $(CLAW_LIB)/saveqc.f \
  $(CLAW_LIB)/valout.f \
  $(CLAW_LIB)/check.f \
  $(CLAW_LIB)/restrt.f \
  $(CLAW_LIB)/setgauges.f \
  $(CLAW_LIB)/dumpgauge.f \
  $(CLAW_LIB)/quick_sort1.f \
  $(CLAW_LIB)/opendatafile.f \
  $(CLAW_LIB)/init_alloc.f90 \
  $(CLAW_LIB)/restrt_alloc.f90 \
  $(CLAW_LIB)/resize_alloc_static.f90 \
# $(CLAW_LIB)/resize_alloc.f90 \

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


### DO NOT remove this line - make depends on it ###
    
