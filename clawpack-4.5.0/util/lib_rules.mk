#
# Library rules for clawpack
#
# ============================================================================
#  Available rules:
#  	all = Build execubtable
#  	lib = Build main library
#  	new = Clean and then rebuild library
#  	install = Install library in $(INSTALL_PATH)
#   objects = Build all library objects and associated modules
#  	clean = Deletes library objects, module files, and library
#  	clobber = Same as clean
#  	help = Prints this help message
# ============================================================================

# General makefile settings
SHELL = /bin/sh
INSTALL_PROGRAM ?= $(INSTALL)
INSTALL_DATA ?= $(INSTALL) -m 644

# Assume the following lists of source
LIB_SRC ?=
LIB_MOD_SRC ?=
LIB_MODULES ?=

# Library settings
NAME ?= claw
INSTALL_PATH ?= $(CLAW)/lib
MOD_PATH ?= $(curdir)
MAJOR_VERSION ?= 4
MINOR_VERSION ?= 5.0

# ============================================================================
# Default compiler and settings
FC ?= gfortran
LINK ?= $(FC)

# Base compiler and linker flags
ALL_FFLAGS =
ALL_LFLAGS =
ALL_INCLUDE =

# Dynamic library flags (archetecture dependent), note the platforms here are
# found by sys.platform in ptyhon, check
UNAME = $(shell uname)
ifeq ($(UNAME), Linux)
	LIB_FULL_NAME = lib$(NAME).so.$(MAJOR_VERSION).$(MINOR_VERSION)
	LIB_MAJOR_NAME = lib$(NAME).so.$(MAJOR_VERSION)
	LIB_SHORT_NAME = lib$(NAME).so
	LIB_INSTALL_NAME = -Wl,-soname,$(LIB_MAJOR_NAME)
	VERSION_FLAGS =
	
	# Compiler specific flags
	ifeq ($(FC),gfortran)
		ALL_FFLAGS += -J$(MOD_PATH) -fPIC
		ALL_LFLAGS += -shared
	endif
else 
	ifeq ($(UNAME), Darwin)
		# Darwin, has a few more options and calls its libraries something
		# different but is taken care of here.  Note that we have to use the 
		# flag -flat_namespace in order for it to act like linux, see the ld 
		# man page for more information on the two-level namespace OS X uses
		LIB_FULL_NAME = lib$(NAME).$(MAJOR_VERSION).$(MINOR_VERSION).dylib
		LIB_MAJOR_NAME = lib$(NAME).$(MAJOR_VERSION).dylib
		LIB_SHORT_NAME = lib$(NAME).dylib
		LIB_INSTALL_NAME = -install_name $(LIB_MAJOR_NAME)
		VERSION_FLAGS = -compatibility_version $(MAJOR_VERSION) -current_version $(MAJOR_VERSION).$(MINOR_VERSION)
	
		# Compiler specific flags
		ifeq ($(FC),gfortran)
			ALL_FFLAGS += -J$(MOD_PATH) -fPIC
			ALL_LFLAGS += -dynamiclib -flat_namespace -undefined suppress
		else 
			ifeq ($(FC),ifort) 
				ALL_FFLAGS += -module $(MOD_PATH)
				ALL_LFLAGS += -dynamiclib -flat_namespace -undefined suppress
			endif
		endif
	endif
endif
# Loop through INCLUDE list to produce ALL_INCLUDE
ALL_INCLUDE = $(addprefix -I, $(INCLUDE) . $(MOD_PATH))

# Add on user definitions at the end
FFLAGS ?=
LFLAGS ?=
ALL_FFLAGS += $(FFLAGS)
ALL_LFLAGS += $(LFLAGS)

# ============================================================================
# Fortran rules and file suffixes
.SUFFIXES:
.SUFFIXES: .f90 .f .mod .o
%.o : %.f90 ; $(FC) -c -o $@ $(ALL_FFLAGS) $(ALL_INCLUDE) $<
%.o : %.f ; $(FC) -c -o $@ $(ALL_FFLAGS) $(ALL_INCLUDE) $<

# ============================================================================
# Convert the source lists and lib name to object lists and correct paths
LIB_OBJ = $(subst .f,.o, $(subst .f90,.o, $(LIB_SRC)))
LIB_MOD_OBJ = $(subst .f90,.o,$(LIB_MOD_SRC))
LIB_MOD_FILES = $(addprefix $(LIB_INCLUDE),$(addsuffix .mod,$(LIB_MODULES)))

# ============================================================================
# Targets
.PHONY.: all, lib, new, install, clean, clobber, help

all: $(LIB_FULL_NAME)
	$(MAKE) install

lib: $(LIB_FULL_NAME)
	
new:
	$(MAKE) clean 
	$(MAKE) $(LIB_FULL_NAME)

install: $(INSTALL_PATH)/$(LIB_FULL_NAME)
	
objects: $(LIB_MOD_OBJ) $(LIB_MOD_FILES) $(LIB_OBJ)

$(INSTALL_PATH)/$(LIB_FULL_NAME): $(LIB_FULL_NAME)
	@if [ ! -d ${INSTALL_PATH} ]; then \
	  mkdir $(INSTALL_PATH); \
	fi
	-cp $(LIB_FULL_NAME) $(INSTALL_PATH)
	@if [ ! -e ${INSTALL_PATH}/${LIB_MAJOR_NAME} ]; then \
	  ln -s $(INSTALL_PATH)/$(LIB_FULL_NAME) $(INSTALL_PATH)/$(LIB_MAJOR_NAME); \
	fi
	@if [ ! -e ${INSTALL_PATH}/${LIB_SHORT_NAME} ]; then \
	  ln -s $(INSTALL_PATH)/$(LIB_FULL_NAME) $(INSTALL_PATH)/$(LIB_SHORT_NAME); \
	fi

$(LIB_FULL_NAME): $(LIB_MOD_OBJ) $(LIB_MOD_FILES) $(LIB_OBJ)  
	$(LINK) $(LIB_INSTALL_NAME) $(VERSION_FLAGS) $(LIB_MOD_OBJ) $(LIB_OBJ) -o $(LIB_FULL_NAME) $(ALL_INCLUDE) $(ALL_LFLAGS)

$(LIB_MOD_FILES): $(LIB_MOD_OBJ)

clean::
	-rm -f $(LIB_OBJ) $(LIB_MOD_OBJ) $(LIB_MOD_FILES) $(LIB_FULL_NAME)
	
clobber:: clean

help::
	@echo "Available rules:"
	@echo "  all = Build execubtable"
	@echo "  lib = Build main library"
	@echo "  new = Clean and then rebuild library"
	@echo "  install = Install library in $(INSTALL_PATH)"
	@echo "  objects = Build all library objects and associated modules"
	@echo "  clean = Deletes library objects, module files, and library"
	@echo "  clobber = Same as clean"
	@echo "  help = Prints this help message"

### DO NOT remove this line - make depends on it ###