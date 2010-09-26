#
# Local application rules for clawpack
#   This rule file should be applied to an actual executable which will link
#   to one of the clawpack packages.
#
# ============================================================================
#  Available rules:
#   $(EXECUTABLE) = Build executable
#   lib = Build required libraries
#   new = Clean source directories include libraries and recompile
#   objects = Make library and local object and module files
#   htmls = Build html help files
#   data = Build all data files
#   output = Run executable and produce output defined in setrun.py
#   plots = Create plots for the output of the executable
#   clean = Remove object and module files along with the compiled binary
#   clobber = Perform clean on library as well as local directory
#   all = Rule for making all htmls and plotting the results of a run
#   help = Display this help message
# ============================================================================

# General makefile settings
SHELL = /bin/sh
INSTALL_PROGRAM ?= $(INSTALL)
INSTALL_DATA ?= $(INSTALL) -m 644

# Clawpack settings
CLAW_PKG ?= claw
LIB_DIM ?= 1d
# Can be defined if you want a particular version
# VERSION ?= 4.4.0
EXECUTABLE ?= x$(CLAW_PKG)
SETRUN_FILE ?= setrun.py
SETPLOT_FILE ?= setplot.py
OUT_DIR ?= _output
PLOT_DIR ?= _plots
MOD_PATH ?= $(curdir)

# Assume the following lists of source and data files
SRC ?=
MOD_SRC ?=
MODULES ?=
DATA_FILES ?= claw.data

# Convert the source lists into object and .mod files
OBJ = $(subst .f,.o, $(subst .f90,.o, $(SRC)))
MOD_OBJ = $(subst .f90,.o,$(MOD_SRC))
MOD_FILES = $(addprefix $(MOD_PATH),$(addsuffix .mod,$(MODULES)))

# Fortran compiler:  FC may be set as an environment variable or in make
# file that 'includes' this one.
FC ?= gfortran
LINK ?= $(FC)
# Path to version of python to use:  May need to use something other than
# the system default in order for plotting to work.  Can set PY as
# environment variable or in make file that 'includes' this one.
PY ?= python
PYTHON = $(PY)

# ============================================================================
# Flags for compiling and linking, all settings can be overridden via the 
# FFLAGS, LFLAGS, and INCLUDE variables
ALL_FFLAGS ?=
ALL_LFLAGS ?=
ifdef INCLUDE
	ALL_INCLUDE = -I$(INCLUDE)
else
	ALL_INCLUDE = 
endif
ifeq ($(FC),gfortran)
	ALL_FFLAGS += -J$(MOD_PATH)
else 
	ifeq ($(FC),ifort)
		ALL_FFLAGS += -module $(MOD_PATH)
	endif
endif

# Library paths, geoclaw needs two libraries
LIB_PATH ?= $(CLAW)/lib
LIB_SRC ?= $(CLAW)/$(CLAW_PKG)/$(LIB_DIM)/lib
ifdef VERSION
	LIB_NAME = $(CLAW_PKG)$(LIB_DIM).$(VERSION)
else
	LIB_NAME = $(CLAW_PKG)$(LIB_DIM)
endif
ifeq ($(CLAW_PKG),claw)
	ALL_LFLAGS += -l$(LIB_NAME) -L$(LIB_PATH)
	ALL_INCLUDE += -I. -I$(MOD_PATH) -I$(LIB_SRC)  
endif
ifeq ($(CLAW_PKG),amrclaw)
	ALL_LFLAGS += -l$(LIB_NAME) -L$(LIB_PATH)
	ALL_INCLUDE += -I. -I$(MOD_PATH) -I$(LIB_SRC) 
endif
ifeq ($(CLAW_PKG),geoclaw)
	AMRCLAW_PATH = $(CLAW)/amrclaw/2d/lib
	ALL_LFLAGS += -l$(LIB_NAME) -lamrclaw2d -L$(LIB_PATH) 
	ALL_INCLUDE += -I. -I$(MOD_PATH) -I$(LIB_SRC) -I$(AMRCLAW_PATH) 
endif


# Add libraries

# Add on user definitions
ifdef FFLAGS
	ALL_FFLAGS += $(FFLAGS)
endif
ifdef LFLAGS
	ALL_LFLAGS += $(LFLAGS)
endif

#----------------------------------------------------------------------------
# Valid suffixes
.SUFFIXES:
.SUFFIXES: .f90 .f .mod .o .data

# how to make .o files from .f files:
%.o : %.f90 ; $(FC) -c  $< -o $@ $(ALL_INCLUDE) $(ALL_FFLAGS)
%.o : %.f ; $(FC) -c $< -o $@ $(ALL_INCLUDE) $(ALL_FFLAGS) 

# Targets that do not correspond to file names:
.PHONY: lib, new, htmls, data, output, plots, clean, cleanlib, clobber, all	

#-----------------------------------------------------------------------------
# Rules for building the executable

$(EXECUTABLE): lib $(MOD_OBJ) $(MOD_FILES) $(OBJ) $(MAKEFILE_LIST) ;
	$(LINK) $(MOD_OBJ) $(OBJ) -o $(EXECUTABLE) $(ALL_INCLUDE) $(ALL_LFLAGS)

lib:
ifeq ($(CLAW_PKG),geoclaw)
	$(MAKE) -C $(AMRCLAW_PATH) lib
	$(MAKE) -C $(AMRCLAW_PATH) install
endif
	$(MAKE) -C $(LIB_SRC) lib
	$(MAKE) -C $(LIB_SRC) install
	
new: clean
ifeq ($(CLAW_PKG),geoclaw)
	$(MAKE) -C $(AMRCLAW_PATH) new
	$(MAKE) -C $(AMRCLAW_PATH) install
endif
	$(MAKE) -C $(LIB_SRC) new
	$(MAKE) -C $(LIB_SRC) install
	$(MAKE) $(EXECUTABLE)
	
objects: lib $(MOD_OBJ) $(OBJ)

#-----------------------------------------------------------------------------
# Rules for building html help files

# Command to create *.html files from *.f etc:
CC2HTML = $(PYTHON) $(CLAW)/doc/clawcode2html.py --force 

# make list of html files to be created by 'make .htmls':
HTML = \
  $(subst .f,.f.html,$(wildcard *.f)) \
  $(subst .f95,.f95.html,$(wildcard *.f95)) \
  $(subst .m,.m.html,$(wildcard *.m)) \
  $(subst .py,.py.html,$(wildcard *.py)) \
  $(subst .data,.data.html,$(wildcard *.data)) \
  $(subst .txt,.html,$(wildcard *.txt)) \
  Makefile.html

# Rules to make html files:  
# e.g. qinit.f --> qinit.f.html
%.f.html : %.f ; $(CC2HTML) $<              
%.f95.html : %.f95 ; $(CC2HTML) $<
%.m.html : %.m ; $(CC2HTML) $<
%.py.html : %.py ; $(CC2HTML) $<
%.data.html : %.data ; $(CC2HTML) $<
Makefile.html : Makefile ; $(CC2HTML) $<    
# drop .txt extension, e.g. README.txt --> README.html
%.html : %.txt ; $(CC2HTML) --dropext $<    

htmls: $(HTML) ;

#----------------------------------------------------------------------------

# Make data files needed by Fortran code:

$(DATA_FILES): $(SETRUN_FILE)
	$(PYTHON) $(SETRUN_FILE)

data: $(SETRUN_FILE) $(DATA_FILES) $(MAKEFILE_LIST)

#----------------------------------------------------------------------------
# Run the code and put fort.* files into subdirectory named output:
# runclaw will execute setrun.py to create data files and determine
# what executable to run, e.g. xclaw or xamr.
output: $(EXECUTABLE) $(DATA_FILES) $(MAKEFILE_LIST)
	$(PYTHON) $(CLAW)/python/pyclaw/runclaw.py  $(EXECUTABLE) $(OUT_DIR)

#----------------------------------------------------------------------------
# Python plotting rules

# Rule to make the plots into subdirectory specified by CLAW_PLOTDIR,
# using data in subdirectory specified by CLAW_OUTDIR and the plotting
# commands specified in CLAW_setplot_file.
ifdef COMSPEC   
  PLOT_CMD := C:\Python25\python.exe C:\cygwin$(CLAW)\python\pyclaw\plotters\plotclaw.py
else
  PLOT_CMD := $(PYTHON) $(CLAW)/python/pyclaw/plotters/plotclaw.py
endif

plots: $(SETPLOT_FILE) $(MAKEFILE_LIST) ;
	$(PLOT_CMD) $(OUT_DIR) $(PLOT_DIR) $(SETPLOT_FILE)

#----------------------------------------------------------------------------

# Clean up options:

clean:
	-rm -f $(MOD_OBJ) $(MOD_FILES) $(OBJ) $(EXECUTABLE)

cleanlib:
	$(MAKE) -C $(LIB_SRC) clean

clobber: clean cleanlib
	-rm -f -r $(OUT_DIR) $(PLOT_DIR) $(DATA_FILES) $(HTML)
	-rm -f -r eagleplots
	-rm -f fort.* *.pyc *.db pyclaw.log
	
#----------------------------------------------------------------------------

# Do multiple things:

all: htmls output plots

#----------------------------------------------------------------------------

help: 
	@echo '   "make ','$(EXECUTABLE),'" = Build executable'
	@echo '   "make lib" = Build required libraries'
	@echo '   "make new" = Clean source directories include libraries and recompile'
	@echo '   "make objects" = Make library and local object and module files'
	@echo '   "make htmls" = Build html help files'
	@echo '   "make data" = Build all data files'
	@echo '   "make output" = Run executable and produce output defined in setrun.py'
	@echo '   "make plots" = Create plots for the output of the executable'
	@echo '   "make clean" = Remove object and module files along with the compiled binary'
	@echo '   "make clobber" = Perform clean on library as well as local directory'
	@echo '   "make all" = Rule for making all htmls and plotting the results of a run'
	@echo '   "make help" = Display this help message'

# Previous make rules for compatibility
# .help: help
# 
# .objs: objects
# 	
# .exe: $(CLAW_EXE)
# 	
# .data: $(DATA_FILES)
# 	
# .output: output
# 
# .plots: plots
# 	
# .program: program
# 	
# .htmls: htmls