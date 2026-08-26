c     Python-generated parton_lum_*.f files call this external entry point.
c     Everything below this boundary is implemented in Fortran modules.
      double precision function pdg2pdf(ih,ipdg,ibeam,x,xmu)
      use pdg2pdf_internal_module, only: pdg2pdf_internal_value
      use fnlo_runtime_common, only: pdlabel
      implicit none
      integer ih,ipdg,ibeam
      double precision x,xmu

      pdg2pdf=pdg2pdf_internal_value(ih,ipdg,ibeam,x,xmu,pdlabel)
      end
