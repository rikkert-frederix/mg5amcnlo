c Python-generated process and run-card boundary for fixed_order_user_hooks.

      logical function dummy_cuts(p,istatus,ipdg)
      use fixed_order_user_hooks, only: accept_dummy_cuts
      implicit none
      include 'nexternal.inc'
      double precision p(0:4,nexternal)
      integer istatus(nexternal),ipdg(nexternal)

      dummy_cuts=accept_dummy_cuts()
      end


      double precision function user_dynamical_scale(p)
      use fixed_order_user_hooks, only: fixed_user_scale
      implicit none
      include 'nexternal.inc'
      include 'run.inc'
      double precision p(0:3,nexternal)
      character*80 temp_scale_id
      common/ctemp_scale_id/temp_scale_id

      user_dynamical_scale=fixed_user_scale(muR_ref_fixed,temp_scale_id)
      end
