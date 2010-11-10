#FFLAGS = -O0 -g -qextname
#F2PYFLAGS = --fcompiler=bgpgnu --compiler=mpicc

SOURCES = \
  step1.f \
  rp1.f

LIBSOURCES = \
  $(PETCLAW)/src/fortran/1d/limiter.f \
  $(PETCLAW)/src/fortran/1d/philim.f

LIB = $(CURDIR)/

step1.so : $(SOURCES) $(LIBSOURCES)
	f2py -m step1 -c step1.f rp1.f $(LIBSOURCES)
	
clean:
	rm -f *.so *.log
