c     Generated-include and COMMON boundary for open_output_files.f90.
      subroutine HwU_write_file
      use open_output_files_module, only: HwU_write_file_impl
      implicit none
      include "reweight_pineappl.inc"
      include "pineappl_common.inc"
      logical pineappl
      common /for_pineappl/ pineappl

      call HwU_write_file_impl(pineappl,nh_obs,appl_obs_num)
      return
      end
