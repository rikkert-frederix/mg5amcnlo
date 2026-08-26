      program driver
      use check_poles_module, only: run_check_poles
      use fnlo_process_common, only: p_born
      implicit none

      call run_check_poles(p_born)
      end


      subroutine init_check_poles_data_bridge()
      use check_poles_module, only: initialize_check_poles_data
      use fnlo_process_common, only: nexternal
      implicit none
      include 'coupl.inc'
      double precision zero,pmass(nexternal)
      parameter (zero=0d0)
      include 'pmass.inc'

      call initialize_check_poles_data(pmass)
      return
      end


      subroutine check_poles_set_model_scale(scale_value)
      implicit none
      double precision scale_value
      include 'coupl.inc'

      mu_r=scale_value
      return
      end
