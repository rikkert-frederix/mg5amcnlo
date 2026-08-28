c Legacy OLP ABI for the module-owned dummy implementation.

      subroutine BinothLHA(p_born,born_wgt,virt_wgt)
      use process_dimensions, only: nexternal
      use binoth_lha_dummy, only: module_BinothLHA => BinothLHA
      implicit none
      double precision p_born(0:3,nexternal-1)
      double precision born_wgt,virt_wgt

      call module_BinothLHA(virt_wgt)
      return
      end

      subroutine BinothLHA_factorized(contribution,event_slot,
     $     born_wgt,virt_wgt)
      use binoth_lha_dummy, only: module_BinothLHA => BinothLHA
      implicit none
      integer contribution,event_slot
      double precision born_wgt,virt_wgt
      call module_BinothLHA(virt_wgt)
      return
      end
