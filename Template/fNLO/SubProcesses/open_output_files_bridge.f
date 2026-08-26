c     Generated-include and shared-state boundary for open_output_files.f90.
      subroutine HwU_write_file
      use open_output_files_module, only: HwU_write_file_impl
      use run_state, only: pineappl
      use fnlo_process_common, only: nh_obs,appl_obs_num
      implicit none

      call HwU_write_file_impl(pineappl,nh_obs,appl_obs_num)
      return
      end
