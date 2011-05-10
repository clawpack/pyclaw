# ============================================================================
#  Program:     Riemann Solver Evaluation Test
#  Created:     2010-04-20
#  Author:      Kyle Mandli
# ============================================================================
#  Available rules:
#  Available rules:"
#    all - Build execubtable"
#    $(EXECUTABLE) - Build main executable
#    clean - Deletes local objects
#    help - Prints this help message
# ============================================================================

# Compiler settings and flags
FC ?= gfortran
LINK ?= $(FC)
ALL_INCLUDE = 
ifdef INCLUDE
	ALL_INCLUDE += -I$(INCLUDE)
endif
ALL_FFLAGS = $(ALL_INCLUDE) $(FFLAGS) -openmp
ALL_LFLAGS = $(ALL_INCLUDE) $(LFLAGS) -openmp
EXECUTABLE = rp_test

# ============================================================================
# Fortran rules and file suffixes
.SUFFIXES:
.SUFFIXES: .f90 .mod .o
%.o : %.f90 ; $(FC) -c -o $@ $(ALL_FFLAGS) $<

# ============================================================================
# Source lists
LOCAL_SRC = main.f90

# ============================================================================
# Object lists
LOCAL_OBJ = $(subst .f90,.o, $(LOCAL_SRC))

# ============================================================================
# Targets
.PHONY.: all, clean, help

all: $(EXECUTABLE)

$(EXECUTABLE): $(LOCAL_OBJ)
	$(LINK) $(LOCAL_OBJ) -o $(EXECUTABLE) $(ALL_LFLAGS)

clean:
	-rm -f $(LOCAL_OBJ) $(EXECUTABLE) *.mod

help:
	@echo "Available rules:"
	@echo "  all = Build execubtable"
	@echo "  "$(EXECUTABLE)" = Build main executable"
	@echo "  clean = Deletes local objects"
	@echo "  help = Prints this help message"

### DO NOT remove this line - make depends on it ###