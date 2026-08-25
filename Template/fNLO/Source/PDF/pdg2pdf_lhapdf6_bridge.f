c     Python-generated parton_lum_*.f files call this external entry point.
      double precision function pdg2pdf(ih,ipdg,ibeam,x,xmu)
      use pdg2pdf_lhapdf6_module, only: pdg2pdf_lhapdf6_value
      implicit none
      integer ih,ipdg,ibeam
      double precision x,xmu

      pdg2pdf=pdg2pdf_lhapdf6_value(ih,ipdg,ibeam,x,xmu)
      end
