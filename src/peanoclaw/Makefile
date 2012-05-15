UNAME := $(shell uname)

# Default settings
#------------------
ifndef DIMENSIONS
DIMENSION := 2
endif
ifndef BUILD_MODE
# Options: debug, release
BUILD_MODE := release
endif

# Get OS settings
#-----------------
ifeq ($(UNAME), Linux)
LIB_EXTENSION := so
endif
ifeq ($(UNAME), Darwin)
LIB_EXTENSION := lib
endif

# Use Scons
all:
ifdef PEANOCLAW_SRC
	echo Using Scons as build system
	cd $(PEANOCLAW_SRC); scons target=peanoclaw build=$(BUILD_MODE) dim=$(DIMENSION) -j2
	cp $(PEANOCLAW_SRC)/build/peanoclaw/$(BUILD_MODE)/dim$(DIMENSION)/multicore_no/parallel_no/gcc/libpeano-claw-2d.$(LIB_EXTENSION) .
endif

clean:
ifdef PEANOCLAW_SRC
	echo Using Scons as build system
	cd $(PEANOCLAW_SRC); scons target=peanoclaw build=$(BUILD_MODE) dim=$(DIMENSION) -c
endif

#Use Make
all-with-make:
ifdef PEANOCLAW_SRC
	echo Using Make as build system
	cp $(PEANOCLAW_SRC)/peano/applications/peanoclaw/makefile $(PEANOCLAW_SRC)/makefile-peanoclaw
	cd $(PEANOCLAW_SRC); make init; make build
	cp $(PEANOCLAW_SRC)/libpeano-claw-2d.$(LIB_EXTENSION) .
endif

clean-with-make:
ifdef PEANOCLAW_SRC
	echo Using Make as build system
	cp $(PEANOCLAW_SRC)/peano/applications/peanoclaw/makefile $(PEANOCLAW_SRC)/makefile-peanoclaw
	cd $(PEANOCLAW_SRC); make clean
endif

