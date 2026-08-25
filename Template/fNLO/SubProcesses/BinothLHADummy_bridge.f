c Legacy OLP ABI for the module-owned dummy implementation.

      subroutine BinothLHA(p_born,born_wgt,virt_wgt)
      use process_dimensions, only: nexternal
      use binoth_lha_dummy, only: module_BinothLHA => BinothLHA
      implicit none
      double precision p_born(0:3,nexternal-1)
      double precision born_wgt,virt_wgt

      call module_BinothLHA(p_born,born_wgt,virt_wgt)
      return
      end


      subroutine BinothLHAInit(filename)
      use binoth_lha_dummy, only:
     &     module_BinothLHAInit => BinothLHAInit
      implicit none
      character*(*) filename

      call module_BinothLHAInit(filename)
      return
      end


      subroutine ctsstatistics(n_mp,n_disc)
      use binoth_lha_dummy, only:
     &     module_ctsstatistics => ctsstatistics
      implicit none
      integer n_mp,n_disc

      call module_ctsstatistics(n_mp,n_disc)
      return
      end


      subroutine sloopmatrix(p_born,virt_wgts)
      use process_dimensions, only: nexternal
      use binoth_lha_dummy, only:
     &     module_sloopmatrix => sloopmatrix
      implicit none
      double precision p_born(0:3,nexternal-1),virt_wgts(3)

      call module_sloopmatrix(p_born,virt_wgts)
      return
      end
