FC = gfortran
FFLAGS = #-qextname
F2PYFLAGS = #--fcompiler=bgpgnu --compiler=mpicc

SOURCES = \
  rp1nel.f

LIBSOURCES = \
  $(PETCLAW)/src/fortran/1d/limiter.f \
  $(PETCLAW)/src/fortran/1d/philim.f

LIB = $(CURDIR)/

# make list of .o files required from the sources above:
OBJECTS = $(subst .f,.o, $(SOURCES))
LIBOBJECTS = $(subst .f,.o, $(LIBSOURCES))

.SUFFIXES: .f .o
%.o : %.f ; $(FC) $(FFLAGS) -c -I$(LIB) $< -o $@ $(FFLAGS) 

step1.so : $(LIBOBJECTS) $(OBJECTS)
	f2py $(F2PYFLAGS) -m step1 -c step1fw.f $(OBJECTS) $(LIBOBJECTS)
	
clean:
	rm -f *.o *.so *.log
	
