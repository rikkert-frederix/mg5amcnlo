include ../Source/make_opts

LIBDIR = ../lib/
LINKLIBS = -L$(LIBDIR) -lgeneric
LIBS = $(LIBDIR)libgeneric.a

ifneq ("$(wildcard ../../Source/RUNNING)","")
    LINKLIBS += -lrunning
    LIBS += $(LIBDIR)librunning.$(libext)
endif

%.o : %.f
	$(FC) $(FFLAGS) -c -I. -I$(LIBDIR) $<
%.o : %.for
	$(FC) $(FFLAGS) -c $<
