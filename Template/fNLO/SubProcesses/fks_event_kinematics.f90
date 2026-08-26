module fks_event_kinematics
  use run_state, only: ebeam, xbk
  use fnlo_process_common, only: ybst_til_tolab, ybst_til_tocm, &
                                 sqrtshat, shat, xbjrk_ev, xbjrk_cnt, &
                                 sqrtshat_ev, shat_ev, sqrtshat_cnt, &
                                 shat_cnt, ycm_ev, ycm_cnt
  implicit none
  private

  public :: set_cms_stuff

contains

  subroutine set_cms_stuff(icountevts)
    implicit none
    integer icountevts

! rapidity of boost from \tilde{k}_1+\tilde{k}_2 c.m. frame to lab frame --
! same for event and counterevents
! This is the rapidity that enters in the arguments of the sinh() and
! cosh() of the boost, in such a way that
!       y(k)_lab = y(k)_tilde - ybst_til_tolab
! where y(k)_lab and y(k)_tilde are the rapidities computed with a generic
! four-momentum k, in the lab frame and in the \tilde{k}_1+\tilde{k}_2
! c.m. frame respectively
    ybst_til_tolab = -ycm_cnt(0) - 0.5d0*log(ebeam(1)/ebeam(2))
    if (icountevts .eq. -100) then
! set Bjorken x's in run.inc for the computation of PDFs in auto_dsig
      xbk(1) = xbjrk_ev(1)
      xbk(2) = xbjrk_ev(2)
! shat=2*k1.k2 -- consistency of this assignment with momenta checked
! in phspncheck_nocms
      shat = shat_ev
      sqrtshat = sqrtshat_ev
! rapidity of boost from \tilde{k}_1+\tilde{k}_2 c.m. frame to
! k_1+k_2 c.m. frame
      ybst_til_tocm = ycm_ev - ycm_cnt(0)
    else
! do the same as above for the counterevents
      xbk(1) = xbjrk_cnt(1, icountevts)
      xbk(2) = xbjrk_cnt(2, icountevts)
      shat = shat_cnt(icountevts)
      sqrtshat = sqrtshat_cnt(icountevts)
      ybst_til_tocm = ycm_cnt(icountevts) - ycm_cnt(0)
    end if
    return
  end subroutine set_cms_stuff

end module fks_event_kinematics
