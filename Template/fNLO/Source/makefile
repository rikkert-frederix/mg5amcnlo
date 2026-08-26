include make_opts

LIBDIR= ../lib/
BINDIR= ../bin/
PDFDIR= ./PDF/
CUTTOOLSDIR= ./CutTools/
IREGIDIR= ./IREGI/src/

GENERIC = fnlo_runtime_common.o timing_state.o $(alfas_functions).o \
	  alfas_functions_bridge.o	\
	  run_state.o kin_functions.o run_printout.o				\
	  ranmar.o setrun.o setrun_bridge.o

ifneq ("$(wildcard ./RUNNING)","")
	LINKLIBS += -lrunning
	LIBS += $(LIBDIR)librunning.$(libext) 
endif


.f.o: ; $(FC) $(FFLAGS) -c $*.f
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

all: $(LIBDIR)libdhelas.a $(LIBDIR)libgeneric.a $(LIBDIR)libpdf.a	\
    $(LIBDIR)libmodel.a \
    param_card.inc $(LIBS)
	rm -f PDF/*.o


$(LIBDIR)librunning.$(libext): RUNNING
	rm -f $(LIBDIR)librunning.a
	cd RUNNING; make

$(LIBDIR)libdhelas.a: DHELAS
	rm -f $(LIBDIR)libdhelas.a
	cd DHELAS; make

$(LIBDIR)libgeneric.a: $(GENERIC) extra_weights.o
	rm -f $(LIBDIR)libgeneric.a
	ar cru libgeneric.a $(GENERIC) extra_weights.o
	ranlib libgeneric.a
	mv libgeneric.a $(LIBDIR)
	cp -f *.mod $(LIBDIR)
	rm -f $(alfas_functions).o alfas_functions_bridge.o

timing_state.o run_state.o ranmar.o $(alfas_functions).o: \
	fnlo_runtime_common.o

setrun.o: extra_weights.o run_state.o
setrun_bridge.o: run_card.inc MODEL/coupl.inc extra_weights.o \
		 run_state.o setrun.o

alfas_functions_bridge.o: $(alfas_functions).o
run_printout.o: run_state.o

$(LIBDIR)libpdf.a: $(LIBDIR)libgeneric.a PDF
	rm -f $(LIBDIR)libpdf.a
	cd PDF; make	

$(LIBDIR)libmodel.a: MODEL param_card.inc
	rm -f $(LIBDIR)libmodel.a
	cd MODEL; make	

run_card.inc: ../Cards/run_card.dat
	../bin/aMCatNLO treatcards run

param_card.inc: ../Cards/param_card.dat
	../bin/aMCatNLO treatcards param

CutTools: $(LIBDIR)libcts.a
libcuttools: $(LIBDIR)libcts.a

IREGI: $(LIBDIR)libiregi.a
libiregi: $(LIBDIR)libiregi.a

$(LIBDIR)libcts.a: $(CUTTOOLSDIR)
	cd $(CUTTOOLSDIR); make
	ln -sf ../Source/$(CUTTOOLSDIR)includects/libcts.a $(LIBDIR)libcts.a
	ln -sf ../Source/$(CUTTOOLSDIR)includects/mpmodule.mod $(LIBDIR)mpmodule.mod

$(LIBDIR)libiregi.a: $(IREGIDIR)
	cd $(IREGIDIR); make
	ln -sf ../Source/$(IREGIDIR)libiregi.a $(LIBDIR)libiregi.a

clean4pdf:
	rm -f ../lib/libpdf.a
	rm -f ../lib/libgeneric.a

cleanCT:
	cd $(CUTTOOLSDIR); make clean; cd ..

cleanIR:
	cd $(IREGIDIR); make clean; cd ..

libmodel: $(LIBDIR)libmodel.a

libdhelas: $(LIBDIR)libdhelas.a

clean:
	rm -f *.o *.mod
	rm -f param_card.inc run_card.inc
	cd MODEL; make clean; cd ..
	cd PDF; make clean; cd ..
	cd DHELAS; make clean; cd ..
	if [ -d $(CUTTOOLSDIR) ]; then cd $(CUTTOOLSDIR); make clean; cd ..; fi
	if [ -d $(IREGIDIR) ]; then cd $(IREGIDIR); make clean; cd ..; fi
