
# Makefile for Clawpack code in this directory.
# This version only sets the local files and frequently changed
# options, and then includes the standard makefile pointed to by CLAWMAKE.
CLAWMAKE = $(CLAW)/util/Makefile.common

# See the above file for details and a list of make options, or type
#   make .help
# at the unix prompt.

# Adjust these variables if desired:
# ----------------------------------

CLAW_PKG = Classic                  # Clawpack package to use
CLAW_EXE = xclaw                    # Executable to create
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
  driver.f \
  qinit.f \
  rp1.f \
  setprob.f 


# Clawpack libraries to be used:
CLAW_LIB = $(CLAW)/clawpack/1d/lib
    
CLAW_LIBSOURCES = \
  $(CLAW_LIB)/claw1ez.f \
  $(CLAW_LIB)/bc1.f \
  $(CLAW_LIB)/setaux.f \
  $(CLAW_LIB)/b4step1.f \
  $(CLAW_LIB)/claw1.f \
  $(CLAW_LIB)/step1.f \
  $(CLAW_LIB)/copyq1.f \
  $(CLAW_LIB)/limiter.f \
  $(CLAW_LIB)/philim.f \
  $(CLAW_LIB)/src1.f \
  $(CLAW_LIB)/out1.f \
  $(CLAW_LIB)/opendatafile.f



#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


### DO NOT remove this line - make depends on it ###
    
