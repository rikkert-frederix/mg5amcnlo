module fks_event_kinematics
  use run_state, only: ebeam, xbk
  use fnlo_process_common, only: ybst_til_tolab, ybst_til_tocm, &
                                 sqrtshat, shat, soft_counterevent, &
                                 real_event, event_bjorken_x, &
                                 event_sqrt_shat, event_shat, event_ycm
  implicit none
  private

  public :: set_cms_stuff

contains

  subroutine set_cms_stuff(icountevts)
    implicit none
    integer, intent(in) :: icountevts

    if (icountevts < soft_counterevent .or. &
        icountevts > real_event) then
      write (*, *) 'Invalid event index in set_cms_stuff:', icountevts
      stop 1
    end if

! rapidity of boost from \tilde{k}_1+\tilde{k}_2 c.m. frame to lab frame --
! same for event and counterevents
! This is the rapidity that enters in the arguments of the sinh() and
! cosh() of the boost, in such a way that
!       y(k)_lab = y(k)_tilde - ybst_til_tolab
! where y(k)_lab and y(k)_tilde are the rapidities computed with a generic
! four-momentum k, in the lab frame and in the \tilde{k}_1+\tilde{k}_2
! c.m. frame respectively
    ybst_til_tolab = -event_ycm(soft_counterevent) - &
                     0.5d0*log(ebeam(1)/ebeam(2))
! set Bjorken x's in run.inc for the computation of PDFs in auto_dsig
    xbk(1) = event_bjorken_x(1, icountevts)
    xbk(2) = event_bjorken_x(2, icountevts)
! shat=2*k1.k2 -- consistency of this assignment with momenta checked
! in phspncheck_nocms
    shat = event_shat(icountevts)
    sqrtshat = event_sqrt_shat(icountevts)
! rapidity of boost from \tilde{k}_1+\tilde{k}_2 c.m. frame to
! k_1+k_2 c.m. frame
    ybst_til_tocm = event_ycm(icountevts) - &
                    event_ycm(soft_counterevent)
    return
  end subroutine set_cms_stuff

end module fks_event_kinematics
