c     Python-generated parton_lum_*.f files call these external entry points.
      double precision function pdg2pdf_timed(ih,ipdg,ibeam,x,xmu)
      use pdg2pdf_lhapdf6_module, only: pdg2pdf_lhapdf6_value
      use timing_state, only: tr_pdf
      implicit none
      integer ih,ipdg,ibeam
      double precision x,xmu
      real time_before,time_after

      call cpu_time(time_before)
      pdg2pdf_timed=pdg2pdf_lhapdf6_value(ih,ipdg,ibeam,x,xmu)
      call cpu_time(time_after)
      tr_pdf=tr_pdf+(time_after-time_before)
      end


      double precision function pdg2pdf(ih,ipdg,ibeam,x,xmu)
      use pdg2pdf_lhapdf6_module, only: pdg2pdf_lhapdf6_value
      implicit none
      integer ih,ipdg,ibeam
      double precision x,xmu

      pdg2pdf=pdg2pdf_lhapdf6_value(ih,ipdg,ibeam,x,xmu)
      end
