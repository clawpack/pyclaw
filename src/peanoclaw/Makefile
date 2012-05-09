UNAME := $(shell uname)

# Get OS settings
#-----------------
ifeq ($(UNAME), Linux)
LIB_EXTENSION := so
endif
ifeq ($(UNAME), Darwin)
LIB_EXTENSION := lib
endif

# Check wether make or scons should be used
#-------------------------------------------
ifneq ("$(wildcard $(PEANOCLAW_SRC)/makefile1)","")

all:
ifdef PEANOCLAW_SRC
	echo Make
	echo $(wildcard $(PEANOCLAW_SRC)/makefile)
endif

else

# Use Scons
all:
ifdef PEANOCLAW_SRC
	echo Scons
	#cd $(PEANOCLAW_SRC); scons target=peanoclaw build=release -j2;
	#cp $(PEANOCLAW_SRC)/build/peanoclaw/debug/dim2/multicore_no/parallel_no/gcc/libpeano-claw-2d.$(LIB_EXTENSION) .
endif

endif
