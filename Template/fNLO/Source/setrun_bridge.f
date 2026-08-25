      subroutine setrun
c-----------------------------------------------------------------------
c     Generated run-card assignments remain in fixed form.  All static
c     post-processing lives in the free-form setrun_module.
c-----------------------------------------------------------------------
      use extra_weights
      use run_state
      use setrun_module, only: complete_setrun
      implicit none

      call reset_run_state()
      call initialize_extra_weights()
      include 'run_card.inc'
      call complete_setrun()
      end


      subroutine setrun_model_strong_coupling(value)
c-----------------------------------------------------------------------
c     Isolate the Python-generated model declaration from the static
c     Fortran 90 implementation.
c-----------------------------------------------------------------------
      implicit none
      double precision value
      include 'MODEL/coupl.inc'

      value=G
      end
