      subroutine BinothLHA(p,born_wgt,virt_wgt)
      use binoth_lha_user_backend, only: binoth_lha_eval
      implicit none
      include 'nexternal.inc'
      include 'orders.inc'
      double precision p(0:3,nexternal-1)
      double precision born_wgt,virt_wgt
      double precision amp_split_finite(amp_split_size)
      common /to_amp_split_finite/ amp_split_finite

      call binoth_lha_eval(virt_wgt,amp_split,
     &     amp_split_finite)
      return
      end
