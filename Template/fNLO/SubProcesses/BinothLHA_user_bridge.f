      subroutine BinothLHA(p,born_wgt,virt_wgt)
      use binoth_lha_user_backend, only: binoth_lha_eval
      use fnlo_process_common, only: nexternal,amp_split,
     $     amp_split_finite=>amp_split_finite_ml
      implicit none
      double precision p(0:3,nexternal-1)
      double precision born_wgt,virt_wgt
      call binoth_lha_eval(virt_wgt,amp_split,
     &     amp_split_finite)
      return
      end

      subroutine BinothLHA_factorized(contribution,event_slot,
     $     born_wgt,virt_wgt)
      implicit none
      integer contribution,event_slot
      double precision born_wgt,virt_wgt
      write(*,*) 'Factorized density matrices require MadLoop virtuals'
      stop 1
      end
