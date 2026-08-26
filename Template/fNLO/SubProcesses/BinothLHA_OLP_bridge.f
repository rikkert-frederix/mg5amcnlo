      subroutine BinothLHA(pin,born_wgt,virt_wgt)
      use binoth_lha_olp_backend, only: binoth_lha_eval
      use fnlo_process_common, only: nexternal
      implicit none
      include 'coupl.inc'
      include 'Binoth_proc.inc'
      double precision pin(0:3,nexternal-1)
      double precision born_wgt,virt_wgt
      double precision pmass(nexternal),zero
      parameter (zero=0d0)
      include 'pmass.inc'

      call binoth_lha_eval(pin,born_wgt,virt_wgt,proc_label,pmass)
      return
      end


      subroutine binoth_lha_update_couplings(mu_r_value,alpha_s)
      use fnlo_process_common, only: qes2
      implicit none
      include 'coupl.inc'
      double precision mu_r_value,alpha_s,pi
      parameter (pi=3.1415926535897932385d0)
      mu_r=sqrt(qes2)
      call update_as_param()
      alpha_s=g**2/(4d0*pi)
      mu_r_value=mu_r
      return
      end
