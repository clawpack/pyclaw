SOURCES = \
  step1.f \
  rp1.f

LIBSOURCES = \
  $(PETCLAW)/src/fortran/1d/limiter.f \
  $(PETCLAW)/src/fortran/1d/philim.f

INCLUDE = 
LIB ?= $(CURDIR)/

# make list of .o files required from the sources above:
OBJECTS = $(subst .f,.o, $(SOURCES))
LIBOBJECTS = $(subst .f,.o, $(LIBSOURCES))

.SUFFIXES: .f .o
%.o : %.f $(INCLUDE) ; $(FC) -c -I$(LIB) $< -o $@ $(FFLAGS) 

step1.so : $(SOURCES) $(LIBSOURCES) $(OBJECTS)
	f2py -m step1 -c  step1.f rp1.o $(LIBSOURCES)
	
clean:
	rm -f *.o *.so *.log
	
