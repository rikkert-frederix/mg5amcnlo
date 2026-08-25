c     Python-generated parton_lum_*.f files call these external entry points.
c     Everything below this boundary is implemented in Fortran modules.
      double precision function pdg2pdf_timed(ih,ipdg,ibeam,x,xmu)
      use timing_state, only: tr_pdf
      implicit none
      integer ih,ipdg,ibeam
      double precision x,xmu,pdg2pdf
      real time_before,time_after
      external pdg2pdf

      call cpu_time(time_before)
      pdg2pdf_timed=pdg2pdf(ih,ipdg,ibeam,x,xmu)
      call cpu_time(time_after)
      tr_pdf=tr_pdf+(time_after-time_before)
      end


      double precision function pdg2pdf(ih,ipdg,ibeam,x,xmu)
      use pdg2pdf_internal_module, only: pdg2pdf_internal_value
      implicit none
      integer ih,ipdg,ibeam,lhaid,pdfscheme
      character(len=7) pdlabel
      double precision x,xmu
      common /to_pdf/ lhaid,pdfscheme,pdlabel

      pdg2pdf=pdg2pdf_internal_value(ih,ipdg,ibeam,x,xmu,pdlabel)
      end
