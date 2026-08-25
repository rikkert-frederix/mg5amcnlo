# Note that optimization doesn't work on all machines

include ../make_opts

LIBRARY = libcernlib.a
LIBDIR = ../../lib/
OBJECTS = cernlib.o

all: $(LIBDIR)$(LIBRARY)

%.o: %.f90 makefile
	$(FC) $(FFLAGS) -c $<

$(LIBDIR)$(LIBRARY): $(OBJECTS)
	$(RM) $(LIBRARY)
	ar cru $(LIBRARY) $(OBJECTS)
	ranlib $(LIBRARY)
	mv $(LIBRARY) $(LIBDIR)
	cp -f *.mod $(LIBDIR)

clean:
	$(RM) *.o *.mod $(LIBDIR)$(LIBRARY) $(LIBDIR)cernlib_module.mod
