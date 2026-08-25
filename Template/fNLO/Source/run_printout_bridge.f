c Python-generated PDF and run-card state boundary for run_printout_module.

      subroutine run_printout
      use run_printout_module, only: write_run_summary
      implicit none
      include 'PDF/pdf.inc'
      include 'run.inc'
      include 'alfas.inc'

      call write_run_summary(lpp,ebeam,pdlabel,asmz,nloop,
     $     fixed_ren_scale,scale,fixed_fac_scale,q2fact)
      end
