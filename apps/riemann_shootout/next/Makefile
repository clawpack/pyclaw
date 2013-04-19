all: _iso_c_advection.so iso_c_classic1.so

FCFLAGS=-fPIC -fno-second-underscore -Wall -O3 -funroll-loops
LDFLAGS=-Wall -undefined dynamic_lookup -bundle -Wl,-install_name,$@

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
