all: _iso_c_advection.so iso_c_classic1.so

FCFLAGS:=$(shell python -c "import numpy.distutils.fcompiler as fc; fc.log.set_verbosity(-1); f = fc.new_fcompiler(); f.customize(); print ' '.join(f.get_flags_f90())")
LDFLAGS:=$(shell python -c "import numpy.distutils.fcompiler as fc; fc.log.set_verbosity(-1); f = fc.new_fcompiler(); f.customize(); print ' '.join(f.get_flags_linker_so())")

build/f90:
	mkdir -p build/f90

%.o %.mod: %.f90
	gfortran $FCFLAGS -c -o $@ $<

_iso_c_advection.so: iso_c/rp1_advection.f90
	gfortran ${FCFLAGS} ${LDFLAGS} -o $@ $<

iso_c_classic1.so: iso_c/philim.f90 iso_c/limiter.f90 iso_c/step1.f90
	gfortran ${FCFLAGS} ${LDFLAGS} -o $@ $^

clean:
	rm _iso_c_advection.so
	rm iso_c_classic1.so

show_flags:
	@echo LDFLAGS=$(LDFLAGS)
	@echo FCFLAGS=$(FCFLAGS)
