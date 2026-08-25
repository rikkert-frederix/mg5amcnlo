      subroutine BinothLHA(p,born_wgt,virt_wgt)
      use binoth_lha_user_backend, only: binoth_lha_eval
      implicit none
      include 'nexternal.inc'
      include 'orders.inc'
      double precision p(0:3,nexternal-1)
      double precision born_wgt,virt_wgt
      double precision amp_split_finite(amp_split_size)
      common /to_amp_split_finite/ amp_split_finite

      call binoth_lha_eval(p,born_wgt,virt_wgt,amp_split,
     &     amp_split_finite)
      return
      end


      subroutine binoth_lha_update_couplings(mu_r_value,alpha_s)
      implicit none
      include 'coupl.inc'
      double precision mu_r_value,alpha_s,pi
      parameter (pi=3.1415926535897932385d0)
      double precision qes2
      common /coupl_es/ qes2

      mu_r=sqrt(qes2)
      call update_as_param()
      alpha_s=g**2/(4d0*pi)
      mu_r_value=mu_r
      return
      end
