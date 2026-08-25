# ----------------------------------------------------------------------------
#
# Makefile for PDF library
# Apr. 25 2003
#
# ----------------------------------------------------------------------------

LIBRARY	      = libpdf.$(libext)
LIBDIR        = ../../lib/

include ../make_opts

%.o: %.f90 makefile
	$(FC) $(FFLAGS) -I$(LIBDIR) -c $<

%.o: %.f makefile
	$(FC) $(FFLAGS) -I$(LIBDIR) -c $<

ifdef epdf
  $(error The fNLO template no longer supports ePDF/eMELA)
endif

ifdef lhapdf
  ifneq ($(lhapdfversion),6)
    $(error The fNLO template supports LHAPDF 6 only)
  endif
  ifeq ($(lhapdfsubversion),1) # 6.1.X
    PDF         = pdfwrap_lhapdf.o pdf_lhapdf6.o pdg2pdf_lhapdf6.o \
                  pdg2pdf_lhapdf6_bridge.o opendata.o
  else # 6.2.X and later
    CXXFLAGS   += -std=c++11
    PDF         = pdfwrap_lhapdf.o pdf_lhapdf62.o pdg2pdf_lhapdf6.o \
                  pdg2pdf_lhapdf6_bridge.o opendata.o
  endif
else
  PDF         = NNPDFDriver.o pdf.o pdfwrap.o pdfwrap_bridge.o pdg2pdf.o \
                pdg2pdf_bridge.o opendata.o
endif

all: $(LIBDIR)$(LIBRARY)

pdfwrap_bridge.o: pdfwrap.o pdf.inc ../alfas.inc
pdf.o pdfwrap.o: NNPDFDriver.o
pdg2pdf.o: pdf.o
pdg2pdf_bridge.o: pdg2pdf.o
pdg2pdf_lhapdf6_bridge.o: pdg2pdf_lhapdf6.o

$(LIBDIR)$(LIBRARY): $(PDF)
	$(call CREATELIB, $@, $(PDF))

clean:
	$(RM) *.o *.mod $(LIBDIR)$(LIBRARY)
