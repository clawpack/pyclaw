
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


# Environment variable FC should be set to fortran compiler, e.g. gfortran
FC ?= gfortran   # default if not set as environment variable
# Add any desired compiler flags such as -g here:
FFLAGS =

CLAW_LIB = $(shell pwd)

CLAW_SOURCES = \
  $(CLAW_LIB)/claw2ez.f \
  $(CLAW_LIB)/bc2.f \
  $(CLAW_LIB)/setprob.f \
  $(CLAW_LIB)/setaux.f \
  $(CLAW_LIB)/b4step2.f \
  $(CLAW_LIB)/claw2.f \
  $(CLAW_LIB)/step2.f \
  $(CLAW_LIB)/step2ds.f \
  $(CLAW_LIB)/dimsp2.f \
  $(CLAW_LIB)/flux2.f \
  $(CLAW_LIB)/flux2fw.f \
  $(CLAW_LIB)/copyq2.f \
  $(CLAW_LIB)/inlinelimiter.f \
  $(CLAW_LIB)/limiter.f \
  $(CLAW_LIB)/philim.f \
  $(CLAW_LIB)/src2.f \
  $(CLAW_LIB)/rpn2.f \
  $(CLAW_LIB)/rpt2.f \
  $(CLAW_LIB)/cellave.f \
  $(CLAW_LIB)/fdisc.f \
  $(CLAW_LIB)/out2.f \
  $(CLAW_LIB)/out2_hdf.f \
  $(CLAW_LIB)/restart2.f \
  $(CLAW_LIB)/restart2_hdf.f \
  $(CLAW_LIB)/opendatafile.f

default: .objs

#-------------------------------------------------------------------
# Include Makefile containing standard definitions and make options:
include $(CLAWMAKE)


### DO NOT remove this line - make depends on it ###

