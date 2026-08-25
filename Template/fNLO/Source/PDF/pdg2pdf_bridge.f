c     Python-generated parton_lum_*.f files call this external entry point.
c     Everything below this boundary is implemented in Fortran modules.
      double precision function pdg2pdf(ih,ipdg,ibeam,x,xmu)
      use pdg2pdf_internal_module, only: pdg2pdf_internal_value
      implicit none
      integer ih,ipdg,ibeam,lhaid,pdfscheme
      character(len=7) pdlabel
      double precision x,xmu
      common /to_pdf/ lhaid,pdfscheme,pdlabel

      pdg2pdf=pdg2pdf_internal_value(ih,ipdg,ibeam,x,xmu,pdlabel)
      end
