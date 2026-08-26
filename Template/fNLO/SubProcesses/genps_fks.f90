module genps_fks
  use boostwdir2_module, only: boostwdir2, boostwdir2_in_place
  use process_dimensions, only: validate_process_dimensions
  use run_state
  use timing_state, only: tGenPS
  use kin_functions_module, only: dot => dot_impl, rho => rho_impl
  use phase_space_kinematics, only: rotate_invar, phspncheck_nocms, &
                                    phase_space_lambda
  use fks_diagnostics, only: xmom_compare
  use genps_born, only: born_phase_space, generate_born_phase_space, &
                        invalidate_born_phase_space
  ! Generated parameters keep the phase-space normalization bit-identical
  ! to the NLO template's compile-time arithmetic.
  use fnlo_process_common, only: nexternal, nincoming, max_particles, &
                                 max_branch, nocntevents, use_evpr, nbody, &
                                 xi_i_hat_ev, xij_aor, &
                                 i_fks, j_fks, ybst_til_tolab, &
                                 ybst_til_tocm, sqrtshat, shat, &
                                 xi_i_fks_ev, y_ij_fks_ev, p_i_fks_ev, &
                                 p_i_fks_cnt, xi_i_fks_cnt, xiimax_ev, &
                                 xiimax_cnt, xbjrk_ev, xbjrk_cnt, &
                                 sqrtshat_ev, shat_ev, sqrtshat_cnt, &
                                 shat_cnt, ycm_ev, ycm_cnt, &
                                 xinorm_ev, xinorm_cnt, &
                                 veckn_ev, veckbarn_ev, xp0jfks, &
                                 softtest, colltest, xi_i_fks_fix, &
                                 y_ij_fks_fix, tau_lower_bound
  implicit none
  private
  public :: generate_momenta
  public :: initialize_genps_fks_state

! fnlo_process_common owns the process-sized storage. The fixed-form bridge
! gives the radiation module access only to event and counterevent state.
  double precision, pointer :: cnt_momenta(:, :, :) => null()
  double precision, pointer :: cnt_weight(:) => null()
  double precision, pointer :: cnt_psweight(:) => null()
  double precision, pointer :: cnt_jacobian(:) => null()

  double precision, pointer :: born_lab_momenta(:, :) => null()
  double precision, pointer :: born_coll_momenta(:, :) => null()
  double precision, pointer :: born_norad_momenta(:, :) => null()
  double precision, pointer :: event_momenta(:, :) => null()
  double precision, pointer :: particle_masses(:) => null()

  double precision :: massive_xjac_cache = 1d0
  logical :: genps_state_initialized = .false.

contains

  subroutine initialize_genps_fks_state(cnt_momenta_in, cnt_weight_in, &
                                        cnt_psweight_in, cnt_jacobian_in, &
                                        born_lab_momenta_in, born_coll_momenta_in, &
                                        born_norad_momenta_in, event_momenta_in, &
                                        particle_masses_in)
    implicit none
    double precision, target, intent(inout) :: cnt_momenta_in(0:, 1:, 0:)
    double precision, target, intent(inout) :: cnt_weight_in(0:)
    double precision, target, intent(inout) :: cnt_psweight_in(0:)
    double precision, target, intent(inout) :: cnt_jacobian_in(0:)
    double precision, target, intent(inout) :: born_lab_momenta_in(0:, 1:)
    double precision, target, intent(inout) :: born_coll_momenta_in(0:, 1:)
    double precision, target, intent(inout) :: born_norad_momenta_in(0:, 1:)
    double precision, target, intent(inout) :: event_momenta_in(0:, 1:)
    double precision, target, intent(inout) :: particle_masses_in(1:)

    call validate_process_dimensions()
    cnt_momenta => cnt_momenta_in
    cnt_weight => cnt_weight_in
    cnt_psweight => cnt_psweight_in
    cnt_jacobian => cnt_jacobian_in
    born_lab_momenta => born_lab_momenta_in
    born_coll_momenta => born_coll_momenta_in
    born_norad_momenta => born_norad_momenta_in
    event_momenta => event_momenta_in
    particle_masses => particle_masses_in

    call validate_bound_genps_state()
    genps_state_initialized = .true.
  end subroutine initialize_genps_fks_state

  subroutine generate_momenta(ndim, iconfig, wgt, x, p)
    implicit none
    integer, intent(in) :: ndim, iconfig
    double precision, intent(inout) :: wgt, x(99)
    double precision, intent(out) :: p(0:3, nexternal)

    type(born_phase_space) :: born
    real :: tBefore, tAfter
    double precision :: jac
    integer :: i
    logical :: pass
    call require_genps_state()
    call cpu_time(tBefore)

! fNLO supports only the standard event-projection mapping.
    use_evpr = .true.
    call generate_born_phase_space(ndim, iconfig, x, born)

    pass = born%valid
    if (pass) then
      call generate_FKS_kinematics(x, ndim, born%xjac, born%xpswgt, &
                                   born%stot, born%shat, born%sqrtshat, &
                                   born%tau, born%ycm, born%ycmhat, born%xbjrk, &
                                   born%masses, born%external_masses, jac, p, pass)
    end if

    if (.not. pass .or. born%xjac < 0d0) then
      jac = -222d0
      cnt_jacobian(0) = -222d0
      cnt_jacobian(1) = -222d0
      cnt_jacobian(2) = -222d0
      p(0, 1) = -99d0
      do i = 0, 2
        cnt_momenta(0, 1, i) = -99d0
      end do
      call invalidate_born_phase_space()
      nocntevents = .true.
    end if

! Apply any incoming integration weight to the event and counterevents.
    do i = 0, 2
      cnt_jacobian(i) = cnt_jacobian(i)*wgt
    end do
    wgt = wgt*jac

    call cpu_time(tAfter)
    tGenPS = tGenPS + (tAfter - tBefore)
  end subroutine generate_momenta

  subroutine generate_FKS_kinematics(x, ndim, xjac0, xpswgt0, &
    & stot, shat_born, sqrtshat_born, tau_born, ycm_born, ycmhat, &
    & xbjrk_born, m, m_born, jac, p, pass)
    implicit none

    double precision xjac0, xpswgt0, x(99), p(0:3, nexternal), &
      & stot, shat_born, sqrtshat_born, tau_born, ycm_born, ycmhat, jac
    double precision xbjrk_born(2)
    double precision M(-max_branch:max_particles), m_born(nexternal - 1)
    integer ndim
    logical pass

    integer icountevts
    integer ixEi, ixyij, ixpi, imother
    double precision xmrec2, m_j_fks, phi_i_fks, rat_xi, tau, &
      & xi_i_fks, y_ij_fks, xi_i_hat, xiimax, xinorm, xjac, xpswgt, &
      & ycm, xp(0:3, nexternal), xbjrk(2), p_i_fks(0:3)
    integer i

    double precision pi
    parameter(pi=3.1415926535897932d0)

    integer isolsign

! check that the starting PS point is meaningful
    if (xjac0 .lt. 0d0) then
      pass = .false.
      return
    end if
!
! Here we start with the FKS Stuff
!
! icountevts=-100 is the event; 0, 1, and 2 are the counterevents.
! The value 5 skips counterevents for the second massive-emission fold;
! it is assigned only after the event has been stored and never indexes
! a counterevent array.
    icountevts = -100
! Event/counterevent values stay negative until their ordinary context is
! generated.  This also forces a failure if a skipped context is consumed.
    p_i_fks_ev(0) = -1.d0
    xiimax_ev = -1.d0
    do i = 0, 2
      p_i_fks_cnt(0, i) = -1.d0
      xiimax_cnt(i) = -1.d0
      cnt_jacobian(i) = -1.d0
    end do
! set cm stuff to values to make the program crash if not set elsewhere
    ybst_til_tolab = 1.d14
    ybst_til_tocm = 1.d14
    sqrtshat = 0.d0
    shat = 0.d0
! If the collinear counterevent is not generated, this stays zero.
    xij_aor = (0.d0, 0.d0)
!
! These will correspond to the vegas x's for the FKS variables xi_i,
! y_ij and phi_i (changing this also requires changing folding parameters)
    ixEi = ndim - 2
    ixyij = ndim - 1
    ixpi = ndim
!
    imother = min(j_fks, i_fks)
    m_j_fks = particle_masses(j_fks)
!
! For final state j_fks, compute the recoil invariant mass
    if (j_fks .gt. nincoming) then
      call get_recoil(born_lab_momenta, imother, shat_born, xmrec2, pass)
      if (.not. pass) then
        xjac0 = -44
        return
      end if
    end if

! Here is the beginning of the loop over the momenta for the event and
! counter-events. This will fill the xp momenta with the event and
! counter-event momenta.
111 continue
    xjac = xjac0
    xpswgt = xpswgt0
!
! Put the Born momenta in the xp momenta, making sure that the mapping
! is correct; put i_fks momenta equal to zero.
    if (i_fks .gt. 1) then
      xp(:, 1:i_fks - 1) = born_lab_momenta(:, 1:i_fks - 1)
      m(1:i_fks - 1) = m_born(1:i_fks - 1)
    end if
    xp(:, i_fks) = 0d0
    m(i_fks) = 0d0
    if (i_fks .lt. nexternal) then
      xp(:, i_fks + 1:nexternal) = born_lab_momenta(:, i_fks:nexternal - 1)
      m(i_fks + 1:nexternal) = m_born(i_fks:nexternal - 1)
    end if
!
! set-up phi_i_fks
!
    phi_i_fks = 2d0*pi*x(ixpi)
    xjac = xjac*2d0*pi
! To keep track of the special phase-space region with massive j_fks
    isolsign = 0
!
! consider the three cases:
! case 1: j_fks is massless final state
! case 2: j_fks is massive final state
! case 3: j_fks is initial state
    if (j_fks .gt. nincoming) then
      shat = shat_born
      sqrtshat = sqrtshat_born
      tau = tau_born
      ycm = ycm_born
      xbjrk(1) = xbjrk_born(1)
      xbjrk(2) = xbjrk_born(2)
      if (m_j_fks .eq. 0d0) then
        isolsign = 1
        call generate_momenta_massless_final(icountevts, i_fks, j_fks &
          & , born_lab_momenta(0:3, imother), shat, sqrtshat, x(ixEi), xmrec2, xp &
          & , phi_i_fks, xiimax, xinorm, xi_i_fks, y_ij_fks, xi_i_hat &
          & , p_i_fks, xjac, xpswgt, pass)
        if (.not. pass) goto 112
      elseif (m_j_fks .gt. 0d0) then
        call generate_momenta_massive_final(icountevts, isolsign &
          & , rat_xi, i_fks, j_fks, born_lab_momenta(0:3, imother) &
          & , shat, sqrtshat, m_j_fks, x(ixEi), xmrec2, xp, phi_i_fks &
          & , xiimax, xinorm, xi_i_fks, y_ij_fks, xi_i_hat, p_i_fks, xjac &
          & , xpswgt, pass)
        if (.not. pass) goto 112
      end if
    elseif (j_fks .le. nincoming) then
      isolsign = 1
      call generate_momenta_initial(icountevts, i_fks, j_fks, xbjrk_born &
        & , tau_born, ycm_born, ycmhat, shat_born, phi_i_fks, xp, x(ixEi) &
        & , shat, stot, sqrtshat, tau, ycm, xbjrk, p_i_fks, xiimax, xinorm &
        & , xi_i_fks, y_ij_fks, xi_i_hat, xpswgt, xjac, pass)
      if (.not. pass) goto 112
    else
      write (*, *) 'Error #2 in genps_fks.f', j_fks
      stop
    end if
! At this point, the phase space lacks a factor xi_i_fks, which need be
! excluded in an NLO computation according to FKS, being taken into
! account elsewhere
!$$$      xpswgt=xpswgt*xi_i_fks
!
! All done, so check four-momentum conservation
    if (xjac .gt. 0.d0) then
      call phspncheck_nocms(nexternal, sqrtshat, m, xp, pass)
      if (.not. pass) then
        xjac = -199
        goto 112
      end if
    end if

    call compute_flux(shat, sqrtshat, m(1), m(2), xpswgt, xjac)
!
112 continue

    call fill_FKS_commons(icountevts, tau, ycm, ycm_born, shat, sqrtshat, xbjrk, &
      & xiimax, xinorm, xi_i_fks, xi_i_hat, p_i_fks, y_ij_fks, xp, p, xjac, jac)
!
    if (icountevts .eq. -100) then
      icountevts = 0
! skips counterevents when integrating over second fold for massive
! j_fks
      if (isolsign .eq. -1) icountevts = 5
    else
      icountevts = icountevts + 1
    end if
    if ((icountevts .le. 2 .and. m_j_fks .eq. 0.d0 .and. (.not. nbody)) .or. &
      & (icountevts .eq. 0 .and. m_j_fks .eq. 0.d0 .and. nbody) .or. &
      & (icountevts .eq. 0 .and. m_j_fks .ne. 0.d0)) then
      goto 111 ! back to the top of the loop

    elseif (icountevts .eq. 5) then
! icountevts=5 only when integrating over the second fold with j_fks
! massive. The counterevents have been skipped, so make sure their
! momenta are unphysical. Born are physical if event was generated, and
! must stay so for the computation of enhancement factors.
      do i = 0, 2
        cnt_jacobian(i) = -299
        cnt_momenta(0, 1, i) = -99
      end do
    end if

    nocntevents = (cnt_jacobian(0) .le. 0.d0) .and. &
      & (cnt_jacobian(1) .le. 0.d0) .and. &
      & (cnt_jacobian(2) .le. 0.d0)
    call xmom_compare(i_fks, j_fks, cnt_jacobian, cnt_momenta, &
                      particle_masses, pass)
!
    return
  end subroutine generate_FKS_kinematics

  subroutine compute_flux(shat, sqrtshat, m1, m2, xpswgt, xjac)
    implicit none
    double precision shat, sqrtshat, m1, m2, xpswgt, xjac

    double precision pwgt, flux

    double precision pi
    parameter(pi=3.1415926535897932d0)

    if (nincoming .eq. 2) then
      flux = 1d0/(2.d0*sqrt(phase_space_lambda(shat, m1**2, m2**2)))
    else                      ! Decays
      flux = 1d0/(2d0*sqrtshat)
    end if
! The pi-dependent factor inserted below is due to the fact that the
! weight computed above is relevant to R_n, as defined in Kajantie's
! book, eq.(III.3.1), while we need the full n-body phase space
    flux = flux/(2d0*pi)**(3*(nexternal - nincoming) - 4)
! This extra pi-dependent factor is due to the fact that the phase-space
! part relevant to i_fks and j_fks does contain all the pi's needed for
! the correct normalization of the phase space
    flux = flux*(2d0*pi)**3
    pwgt = max(xjac*xpswgt, 1d-99)
    xjac = pwgt*flux
!
    return
  end subroutine compute_flux

  subroutine fill_FKS_commons(icountevts, tau, ycm, ycm_born, shat, sqrtshat, xbjrk, &
    & xiimax, xinorm, xi_i_fks, xi_i_hat, p_i_fks, y_ij_fks, xp, p, xjac, jac)

    implicit none
    integer icountevts
    double precision tau, ycm, ycm_born, shat, sqrtshat, xbjrk(2), xiimax, xinorm, &
      & xi_i_fks, xi_i_hat, p_i_fks(0:3), y_ij_fks, xp(0:3, nexternal), p(0:3, nexternal), &
      & xjac, jac

    integer i, j

! Catch the points for which there is no viable phase-space generation
! (still fill the shared state with some information that is needed
! (e.g. ycm_cnt)).
    if (xjac .le. 0d0) then
      xp(0, 1) = -99d0
    end if
!
! Fill shared FKS state
    if (icountevts .eq. -100) then
      ycm_ev = ycm
      shat_ev = shat
      sqrtshat_ev = sqrtshat
      xbjrk_ev(1) = xbjrk(1)
      xbjrk_ev(2) = xbjrk(2)
      xiimax_ev = xiimax
      xinorm_ev = xinorm
      xi_i_fks_ev = xi_i_fks
      xi_i_hat_ev = xi_i_hat
      do i = 0, 3
        p_i_fks_ev(i) = p_i_fks(i)
      end do
      y_ij_fks_ev = y_ij_fks
      do i = 1, nexternal
        do j = 0, 3
          p(j, i) = xp(j, i)
          event_momenta(j, i) = xp(j, i)
        end do
      end do
      jac = xjac
    elseif (icountevts .ge. 0 .and. icountevts .le. 2) then
! Special fix in the case the soft counter-events are not generated but
! the Born and real are. (This can happen if ptj>0 in the
! run_card). This fix is needed for set_cms_stuff to work properly.
      if (icountevts .eq. 0) then
        ycm = ycm_born
      end if
      ycm_cnt(icountevts) = ycm
      shat_cnt(icountevts) = shat
      sqrtshat_cnt(icountevts) = sqrtshat
      xbjrk_cnt(1, icountevts) = xbjrk(1)
      xbjrk_cnt(2, icountevts) = xbjrk(2)
      xiimax_cnt(icountevts) = xiimax
      xinorm_cnt(icountevts) = xinorm
      xi_i_fks_cnt(icountevts) = xi_i_fks
      do i = 0, 3
        p_i_fks_cnt(i, icountevts) = p_i_fks(i)
      end do
      do i = 1, nexternal
        do j = 0, 3
          cnt_momenta(j, i, icountevts) = xp(j, i)
        end do
      end do
      cnt_jacobian(icountevts) = xjac
! The following two are obsolete, but remain part of the generated ABI:
! so give some non-physical values
      cnt_weight(icountevts) = -1d99
      cnt_psweight(icountevts) = -1d99
    else
      write (*, *) 'Invalid counterevent index in fill_FKS_commons:', icountevts
      stop 1
    end if

    return
  end subroutine fill_FKS_commons

  subroutine generate_momenta_massless_final(icountevts, i_fks, j_fks &
    & , p_born_imother, shat, sqrtshat, x, xmrec2, xp, phi_i_fks, xiimax &
    & , xinorm, xi_i_fks, y_ij_fks, xi_i_hat, p_i_fks, xjac, xpswgt &
    & , pass)
    implicit none
! arguments
    integer icountevts, i_fks, j_fks
    double precision shat, sqrtshat, x(2), xmrec2, xp(0:3, nexternal) &
      & , y_ij_fks, p_born_imother(0:3), phi_i_fks, xi_i_hat
    double precision xiimax, xinorm, xi_i_fks, p_i_fks(0:3), xjac, xpswgt
    logical pass
! local
    integer i, j
    double precision E_i_fks, x3len_i_fks, x3len_j_fks, x3len_fks_mother &
      & , costh_i_fks, sinth_i_fks, xpifksred(0:3), th_mother_fks &
      & , costh_mother_fks, sinth_mother_fks, phi_mother_fks &
      & , cosphi_mother_fks, sinphi_mother_fks, recoil(0:3), sumrec &
      & , sumrec2, betabst, gammabst, shybst, chybst, chybstmo, xdir(3) &
      & , veckn, veckbarn, xp_mother(0:3), cosphi_i_fks &
      & , sinphi_i_fks
    complex(kind=kind(0d0)) resAoR0
! external
! parameters
    double precision pi
    parameter(pi=3.1415926535897932d0)
    double precision stiny, sstiny, qtiny, ctiny, cctiny
    complex(kind=kind(0d0)) ximag
    parameter(stiny=1d-6)
    parameter(qtiny=1d-7)
    parameter(ctiny=5d-7)
    parameter(ximag=(0d0, 1d0))
!
    pass = .true.
    if (softtest) then
      sstiny = 0.d0
    else
      sstiny = stiny
    end if
    if (colltest) then
      cctiny = 0.d0
    else
      cctiny = ctiny
    end if
!
! set-up y_ij_fks
!
    if ((icountevts .eq. -100 .or. icountevts .eq. 0) .and. &
      & ((.not. softtest) .or. &
      & (softtest .and. y_ij_fks_fix .eq. -2.d0)) .and. &
      & (.not. colltest)) then
! importance sampling towards collinear singularity
! insert here further importance sampling towards y_ij_fks->1
      y_ij_fks = -2d0*(cctiny + (1 - cctiny)*x(2)**2) + 1d0
    elseif ((icountevts .eq. -100 .or. icountevts .eq. 0) .and. &
      & ((softtest .and. y_ij_fks_fix .ne. -2.d0) .or. &
      & colltest)) then
      y_ij_fks = y_ij_fks_fix
    elseif (icountevts .eq. 1 .or. icountevts .eq. 2) then
      y_ij_fks = 1d0
    else
      write (*, *) 'Error #3 in genps_fks.f', icountevts
      stop
    end if
! importance sampling towards collinear singularity
    xjac = xjac*2d0*x(2)*2d0

    call getangles(p_born_imother, &
      & th_mother_fks, costh_mother_fks, sinth_mother_fks, &
      & phi_mother_fks, cosphi_mother_fks, sinphi_mother_fks)
!
! Compute maximum allowed xi_i_fks
    xiimax = 1 - xmrec2/shat
    xinorm = xiimax
!
! Define xi_i_fks
!
    if ((icountevts .eq. -100 .or. icountevts .eq. 1) .and. &
      & ((.not. colltest) .or. &
      & (colltest .and. xi_i_fks_fix .eq. -2.d0)) .and. &
      & (.not. softtest)) then
      if (icountevts .eq. -100) then
! importance sampling towards soft singularity
! insert here further importance sampling towards xi_i_hat->0
        xi_i_hat = sstiny + (1 - sstiny)*x(1)**2
      end if
! in the case of counter events, xi_i_hat is an input to this function
      xi_i_fks = xi_i_hat*xiimax
    elseif ((icountevts .eq. -100 .or. icountevts .eq. 1) .and. &
      & (colltest .and. xi_i_fks_fix .ne. -2.d0) .and. &
      & (.not. softtest)) then
! This is to keep xi_i_hat, rather than xi_i, fixed in the tests.
      if (xi_i_fks_fix .lt. xiimax) then
        xi_i_fks = xi_i_fks_fix*xiimax
      else
        xi_i_fks = xi_i_fks_fix*xiimax
      end if
    elseif ((icountevts .eq. -100 .or. icountevts .eq. 1) .and. &
      & softtest) then
      if (xi_i_fks_fix .lt. 1d0) then
        xi_i_fks = xi_i_fks_fix*xiimax
      else
        xjac = -102
        pass = .false.
        return
      end if
    elseif (icountevts .eq. 2 .or. icountevts .eq. 0) then
      xi_i_fks = 0d0
    else
      write (*, *) 'Error #4 in genps_fks.f', icountevts
      stop
    end if
! remove the following if no importance sampling towards soft
! singularity is performed when integrating over xi_i_hat
    xjac = xjac*2d0*x(1)

! Check that xii is in the allowed range
    if (icountevts .eq. -100 .or. icountevts .eq. 1) then
      if (xi_i_fks .gt. (1 - xmrec2/shat)) then
        xjac = -101
        pass = .false.
        return
      end if
    elseif (icountevts .eq. 0 .or. icountevts .eq. 2) then
! May insert here a check on whether xii<xicut, rather than doing it
! in the cross sections
      continue
    end if
!
! Compute costh_i_fks from xi_i_fks et al.
!
    E_i_fks = xi_i_fks*sqrtshat/2d0
    x3len_i_fks = E_i_fks
    x3len_j_fks = (shat - xmrec2 - 2*sqrtshat*x3len_i_fks)/ &
      & (2*(sqrtshat - x3len_i_fks*(1 - y_ij_fks)))
    x3len_fks_mother = sqrt(x3len_i_fks**2 + x3len_j_fks**2 + &
      & 2*x3len_i_fks*x3len_j_fks*y_ij_fks)
    if (xi_i_fks .lt. qtiny) then
      costh_i_fks = y_ij_fks + shat*(1 - y_ij_fks**2)*xi_i_fks/ &
        & (shat - xmrec2)
      if (abs(costh_i_fks) .gt. 1.d0) costh_i_fks = y_ij_fks
    elseif (1 - y_ij_fks .lt. qtiny) then
      costh_i_fks = 1 - (shat*(1 - xi_i_fks) - xmrec2)**2*(1 - y_ij_fks)/ &
        & (shat - xmrec2)**2
      if (abs(costh_i_fks) .gt. 1.d0) costh_i_fks = 1.d0
    else
      costh_i_fks = (x3len_fks_mother**2 - x3len_j_fks**2 + x3len_i_fks**2) &
        & /(2*x3len_fks_mother*x3len_i_fks)
      if (abs(costh_i_fks) .gt. 1.d0) then
        if (abs(costh_i_fks) .le. (1.d0 + 1.d-5)) then
          costh_i_fks = sign(1.d0, costh_i_fks)
        else
          write (*, *) 'Fatal error #5 in one_tree', &
            & costh_i_fks, xi_i_fks, y_ij_fks, xmrec2
          stop
        end if
      end if
    end if
    sinth_i_fks = sqrt(1 - costh_i_fks**2)
    cosphi_i_fks = cos(phi_i_fks)
    sinphi_i_fks = sin(phi_i_fks)
    xpifksred(1) = sinth_i_fks*cosphi_i_fks
    xpifksred(2) = sinth_i_fks*sinphi_i_fks
    xpifksred(3) = costh_i_fks
!
! The momentum if i_fks and j_fks
!
    xp(0, i_fks) = E_i_fks
    xp(0, j_fks) = sqrt(x3len_j_fks**2)
    p_i_fks(0) = sqrtshat/2d0
    do j = 1, 3
      p_i_fks(j) = sqrtshat/2d0*xpifksred(j)
      xp(j, i_fks) = E_i_fks*xpifksred(j)
      if (j .ne. 3) then
        xp(j, j_fks) = -xp(j, i_fks)
      else
        xp(j, j_fks) = x3len_fks_mother - xp(j, i_fks)
      end if
    end do
!
    call rotate_invar(xp(0, i_fks), xp(0, i_fks), &
      & costh_mother_fks, sinth_mother_fks, &
      & cosphi_mother_fks, sinphi_mother_fks)
    call rotate_invar(xp(0, j_fks), xp(0, j_fks), &
      & costh_mother_fks, sinth_mother_fks, &
      & cosphi_mother_fks, sinphi_mother_fks)
    call rotate_invar(p_i_fks, p_i_fks, &
      & costh_mother_fks, sinth_mother_fks, &
      & cosphi_mother_fks, sinphi_mother_fks)
!
! Now the xp four vectors of all partons except i_fks and j_fks will be
! boosted along the direction of the mother; start by redefining the
! mother four momenta
    do i = 0, 3
      xp_mother(i) = xp(i, i_fks) + xp(i, j_fks)
      if (nincoming .eq. 2) then
        recoil(i) = xp(i, 1) + xp(i, 2) - xp_mother(i)
      else
        recoil(i) = xp(i, 1) - xp_mother(i)
      end if
    end do
    sumrec = recoil(0) + rho(recoil)
    sumrec2 = sumrec**2
    betabst = -(shat - sumrec2)/(shat + sumrec2)
    gammabst = 1/sqrt(1 - betabst**2)
    shybst = -(shat - sumrec2)/(2*sumrec*sqrtshat)
    chybst = (shat + sumrec2)/(2*sumrec*sqrtshat)
! cosh(y) is very often close to one, so define cosh(y)-1 as well
    chybstmo = (sqrtshat - sumrec)**2/(2*sumrec*sqrtshat)
    do j = 1, 3
      xdir(j) = xp_mother(j)/x3len_fks_mother
    end do
! Perform the boost here
    do i = nincoming + 1, nexternal
      if (i .ne. i_fks .and. i .ne. j_fks .and. shybst .ne. 0.d0) &
        & call boostwdir2_in_place(chybst, shybst, chybstmo, xdir, &
        & xp(0, i))
    end do
!
! Collinear limit of <ij>/[ij]. See innerp3.m.
    if ((icountevts .eq. -100 .or. &
      & (icountevts .eq. 1 .and. xij_aor .eq. 0))) then
      resAoR0 = -exp(2*ximag*(phi_mother_fks + phi_i_fks))
! The term O(srt(1-y)) is formally correct but may be numerically large
! Set it to zero
!$$$          resAoR5=-ximag*sqrt(2.d0)*
!$$$       &          sinphi_i_fks*tan(th_mother_fks/2.d0)*
!$$$       &          exp( 2*ximag*(phi_mother_fks+phi_i_fks) )
!$$$          xij_aor=resAoR0+resAoR5*sqrt(1-y_ij_fks)
      xij_aor = resAoR0
    end if
!
! Phase-space factor for (xii,yij,phii)
    veckn = rho(xp(0, j_fks))
    veckbarn = rho(p_born_imother)
!
! Store event-kinematics quantities.
    if (icountevts .eq. -100) then
      veckn_ev = veckn
      veckbarn_ev = veckbarn
      xp0jfks = xp(0, j_fks)
    end if
!
    xpswgt = xpswgt*2*shat/(4*pi)**3*veckn/veckbarn/ &
      & (2 - xi_i_fks*(1 - xp(0, j_fks)/veckn*y_ij_fks))
    xpswgt = abs(xpswgt)
    return
  end subroutine generate_momenta_massless_final

  subroutine generate_momenta_massive_final(icountevts, isolsign &
    & , rat_xi, i_fks, j_fks, p_born_imother, shat &
    & , sqrtshat, m_j_fks, x, xmrec2, xp, phi_i_fks, xiimax, xinorm &
    & , xi_i_fks, y_ij_fks, xi_i_hat, p_i_fks, xjac, xpswgt, pass)
    implicit none
! arguments
    integer icountevts, i_fks, j_fks, isolsign
    double precision shat, sqrtshat, x(2), xmrec2, xp(0:3, nexternal) &
      & , y_ij_fks, p_born_imother(0:3), m_j_fks, phi_i_fks, xi_i_hat
    double precision xiimax, xinorm, xi_i_fks, p_i_fks(0:3), xjac, xpswgt
    logical pass
! local
    integer i, j
    double precision xmj, xmj2, xmjhat, xmhat, xim, cffA2, cffB2, cffC2 &
      & , cffDEL2, xiBm, ximax, xirplus, xirminus, rat_xi, xitmp1 &
      & , E_i_fks, x3len_i_fks, b2m4ac, x3len_j_fks_num, x3len_j_fks_den &
      & , x3len_j_fks, x3len_fks_mother, costh_i_fks, sinth_i_fks &
      & , xpifksred(0:3), recoil(0:3), xp_mother(0:3), sumrec, expybst &
      & , shybst, chybst, chybstmo, xdir(3), veckn, veckbarn, cosphi_i_fks &
      & , sinphi_i_fks, cosphi_mother_fks, costh_mother_fks &
      & , phi_mother_fks, sinphi_mother_fks, th_mother_fks, xitmp2 &
      & , sinth_mother_fks
! external
! parameters
    double precision pi
    parameter(pi=3.1415926535897932d0)
    double precision stiny, sstiny, qtiny, ctiny, cctiny
    parameter(stiny=1d-6)
    parameter(qtiny=1d-7)
    parameter(ctiny=5d-7)
!
    if (colltest .or. &
      & icountevts .eq. 1 .or. icountevts .eq. 2) then
      write (*, *) 'Error #5 in genps_fks.f:'
      write (*, *) &
        & 'This parametrization cannot be used in FS coll limits'
      stop
    end if
!
    pass = .true.
    if (softtest) then
      sstiny = 0.d0
    else
      sstiny = stiny
    end if
    if (colltest) then
      cctiny = 0.d0
    else
      cctiny = ctiny
    end if
!
! set-up y_ij_fks
!
    if ((icountevts .eq. -100 .or. icountevts .eq. 0) .and. &
      & ((.not. softtest) .or. &
      & (softtest .and. y_ij_fks_fix .eq. -2.d0)) .and. &
      & (.not. colltest)) then
! importance sampling towards collinear singularity
! insert here further importance sampling towards y_ij_fks->1
      y_ij_fks = -2d0*(cctiny + (1 - cctiny)*x(2)**2) + 1d0
    elseif ((icountevts .eq. -100 .or. icountevts .eq. 0) .and. &
      & ((softtest .and. y_ij_fks_fix .ne. -2.d0) .or. &
      & colltest)) then
      y_ij_fks = y_ij_fks_fix
    else
      write (*, *) 'Error #6 in genps_fks.f', icountevts
      stop
    end if
! importance sampling towards collinear singularity
    xjac = xjac*2d0*x(2)*2d0

    call getangles(p_born_imother, &
      & th_mother_fks, costh_mother_fks, sinth_mother_fks, &
      & phi_mother_fks, cosphi_mother_fks, sinphi_mother_fks)
!
! Compute the maximum allowed xi_i_fks
!
    xmj = m_j_fks
    xmj2 = xmj**2
    xmjhat = xmj/sqrtshat
    xmhat = sqrt(xmrec2)/sqrtshat
    xim = (1 - xmhat**2 - 2*xmjhat + xmjhat**2)/(1 - xmjhat)
    cffA2 = 1 - xmjhat**2*(1 - y_ij_fks**2)
    cffB2 = -2*(1 - xmhat**2 - xmjhat**2)
    cffC2 = (1 - (xmhat - xmjhat)**2)*(1 - (xmhat + xmjhat)**2)
    cffDEL2 = cffB2**2 - 4*cffA2*cffC2
    xiBm = (-cffB2 - sqrt(cffDEL2))/(2*cffA2)
    ximax = 1 - (xmhat + xmjhat)**2
    if (xiBm .lt. (xim - 1.d-8) .or. xim .lt. 0.d0 .or. xiBm .lt. 0.d0 .or. &
      & xiBm .gt. (ximax + 1.d-8) .or. ximax .gt. 1 .or. ximax .lt. 0.d0) then
      write (*, *) 'WARNING #4 in one_tree', xim, xiBm, ximax
      xjac = -104d0
      pass = .false.
      return
    end if
    if (y_ij_fks .ge. 0.d0) then
      xirplus = xim
      xirminus = 0.d0
    else
      xirplus = xiBm
      xirminus = xiBm - xim
    end if
    xiimax = xirplus
    xinorm = xirplus + xirminus
    rat_xi = xiimax/xinorm
!
! Generate xi_i_fks
!
    if (icountevts .eq. -100 .and. &
      & ((.not. colltest) .or. &
      & (colltest .and. xi_i_fks_fix .eq. -2.d0)) .and. &
      & (.not. softtest)) then
      massive_xjac_cache = 1.d0
      xitmp1 = x(1)
! Map regions (0,A) and (A,1) in xitmp1 onto regions (0,rat_xi) and (rat_xi,1)
! in xi_i_hat respectively. The parameter A is free, but it appears to be
! convenient to choose A=rat_xi
      if (xitmp1 .le. rat_xi) then
        xitmp1 = xitmp1/rat_xi
        massive_xjac_cache = massive_xjac_cache/rat_xi
! importance sampling towards soft singularity
! insert here further importance samplings
        xitmp2 = sstiny + (1 - sstiny)*xitmp1**2
        massive_xjac_cache = massive_xjac_cache*2*xitmp1
        xi_i_hat = xitmp2*rat_xi
        massive_xjac_cache = massive_xjac_cache*rat_xi
        xi_i_fks = xinorm*xi_i_hat
        isolsign = 1
      else
! insert here further importance samplings
        xi_i_hat = xitmp1
        xi_i_fks = -xinorm*xi_i_hat + 2*xiimax
        isolsign = -1
      end if
    elseif (icountevts .eq. -100 .and. &
      & (colltest .and. xi_i_fks_fix .ne. -2.d0) .and. &
      & (.not. softtest)) then
      massive_xjac_cache = 1.d0
      if (xi_i_fks_fix .lt. xiimax) then
        xi_i_fks = xi_i_fks_fix
      else
        xi_i_fks = xi_i_fks_fix*xiimax
      end if
      isolsign = 1
    elseif ((icountevts .eq. -100) .and. &
      & softtest) then
      massive_xjac_cache = 1.d0
      if (xi_i_fks_fix .lt. xiimax) then
        xi_i_fks = xi_i_fks_fix
      else
        xjac = -102
        pass = .false.
        return
      end if
      isolsign = 1
    elseif (icountevts .eq. 0) then
! Keep the event Jacobian cache here for the matching counterevent.
! used for the (real-emission) event
      xi_i_fks = 0d0
      isolsign = 1
    else
      write (*, *) 'Error #7 in genps_fks.f', icountevts
      stop
    end if
    xjac = xjac*massive_xjac_cache
!
    if (isolsign .eq. 0) then
      write (*, *) 'Fatal error #11 in one_tree', isolsign
      stop
    end if
!
! Compute costh_i_fks
!
    E_i_fks = xi_i_fks*sqrtshat/2d0
    x3len_i_fks = E_i_fks
    b2m4ac = xi_i_fks**2*cffA2 + xi_i_fks*cffB2 + cffC2
    if (b2m4ac .le. 0.d0) then
      if (abs(b2m4ac) .lt. 1.d-3) then
        b2m4ac = 0.d0
      else
        write (*, *) 'Fatal error #6 in one_tree'
        write (*, *) b2m4ac, xi_i_fks, cffA2, cffB2, cffC2
        write (*, *) y_ij_fks, xim, xiBm
        stop
      end if
    end if
    x3len_j_fks_num = -xi_i_fks*y_ij_fks* &
      & (1 - xmhat**2 + xmjhat**2 - xi_i_fks) + &
      & (2 - xi_i_fks)*sqrt(b2m4ac)*isolsign
    x3len_j_fks_den = (2 - xi_i_fks*(1 - y_ij_fks))* &
      & (2 - xi_i_fks*(1 + y_ij_fks))
    x3len_j_fks = sqrtshat*x3len_j_fks_num/x3len_j_fks_den
    if (x3len_j_fks .lt. 0.d0) then
      write (*, *) 'WARNING #7 in one_tree', &
        & x3len_j_fks_num, x3len_j_fks_den, xi_i_fks, y_ij_fks
      xjac = -107d0
      pass = .false.
      return
    end if
    x3len_fks_mother = sqrt(x3len_i_fks**2 + x3len_j_fks**2 + &
      & 2*x3len_i_fks*x3len_j_fks*y_ij_fks)
    if (xi_i_fks .lt. qtiny) then
      costh_i_fks = y_ij_fks + (1 - y_ij_fks**2)*xi_i_fks/sqrt(cffC2)
      if (abs(costh_i_fks) .gt. 1.d0) costh_i_fks = y_ij_fks
    else
      costh_i_fks = (x3len_fks_mother**2 - x3len_j_fks**2 + x3len_i_fks**2) &
        & /(2*x3len_fks_mother*x3len_i_fks)
      if (abs(costh_i_fks) .gt. 1.d0 + qtiny) then
        write (*, *) 'Fatal error #8 in one_tree', &
          & costh_i_fks, xi_i_fks, y_ij_fks, xmrec2
        stop
      elseif (abs(costh_i_fks) .gt. 1.d0) then
        costh_i_fks = sign(1d0, costh_i_fks)
      end if
    end if
    sinth_i_fks = sqrt(1 - costh_i_fks**2)
    cosphi_i_fks = cos(phi_i_fks)
    sinphi_i_fks = sin(phi_i_fks)
    xpifksred(1) = sinth_i_fks*cosphi_i_fks
    xpifksred(2) = sinth_i_fks*sinphi_i_fks
    xpifksred(3) = costh_i_fks
!
! Generate momenta for j_fks and i_fks
!
    xp(0, i_fks) = E_i_fks
    xp(0, j_fks) = sqrt(x3len_j_fks**2 + m_j_fks**2)
    p_i_fks(0) = sqrtshat/2d0
    do j = 1, 3
      p_i_fks(j) = sqrtshat/2d0*xpifksred(j)
      xp(j, i_fks) = E_i_fks*xpifksred(j)
      if (j .ne. 3) then
        xp(j, j_fks) = -xp(j, i_fks)
      else
        xp(j, j_fks) = x3len_fks_mother - xp(j, i_fks)
      end if
    end do
!
    call rotate_invar(xp(0, i_fks), xp(0, i_fks), &
      & costh_mother_fks, sinth_mother_fks, &
      & cosphi_mother_fks, sinphi_mother_fks)
    call rotate_invar(xp(0, j_fks), xp(0, j_fks), &
      & costh_mother_fks, sinth_mother_fks, &
      & cosphi_mother_fks, sinphi_mother_fks)
    call rotate_invar(p_i_fks, p_i_fks, &
      & costh_mother_fks, sinth_mother_fks, &
      & cosphi_mother_fks, sinphi_mother_fks)
!
! Now the xp four vectors of all partons except i_fks and j_fks will be
! boosted along the direction of the mother; start by redefining the
! mother four momenta
    do i = 0, 3
      xp_mother(i) = xp(i, i_fks) + xp(i, j_fks)
      if (nincoming .eq. 2) then
        recoil(i) = xp(i, 1) + xp(i, 2) - xp_mother(i)
      else
        recoil(i) = xp(i, 1) - xp_mother(i)
      end if
    end do
!
    sumrec = recoil(0) + rho(recoil)
    if (xmrec2 .lt. 1.d-16*shat) then
      expybst = sqrtshat*sumrec/(shat - xmj2)* &
        & (1 + xmj2*xmrec2/(shat - xmj2)**2)
    else
      expybst = sumrec/(2*sqrtshat*xmrec2)* &
        & (shat + xmrec2 - xmj2 - shat*sqrt(cffC2))
    end if
    if (expybst .le. 0.d0) then
      write (*, *) 'Fatal error #10 in one_tree', expybst
      stop
    end if
    shybst = (expybst - 1/expybst)/2.d0
    chybst = (expybst + 1/expybst)/2.d0
    chybstmo = chybst - 1.d0
!
    do j = 1, 3
      xdir(j) = xp_mother(j)/x3len_fks_mother
    end do
! Boost the momenta
    do i = nincoming + 1, nexternal
      if (i .ne. i_fks .and. i .ne. j_fks .and. shybst .ne. 0.d0) &
        & call boostwdir2_in_place(chybst, shybst, chybstmo, xdir, &
        & xp(0, i))
    end do
!
! Phase-space factor for (xii,yij,phii)
    veckn = rho(xp(0, j_fks))
    veckbarn = rho(p_born_imother)
!
! Store event-kinematics quantities.
    if (icountevts .eq. -100) then
      veckn_ev = veckn
      veckbarn_ev = veckbarn
      xp0jfks = xp(0, j_fks)
    end if
!
    xpswgt = xpswgt*2*shat/(4*pi)**3*veckn/veckbarn/ &
      & (2 - xi_i_fks*(1 - xp(0, j_fks)/veckn*y_ij_fks))
    xpswgt = abs(xpswgt)
    return
  end subroutine generate_momenta_massive_final

  subroutine generate_momenta_initial(icountevts, i_fks, j_fks, &
    & xbjrk_born, tau_born, ycm_born, ycmhat, shat_born, phi_i_fks, xp, x &
    & , shat, stot, sqrtshat, tau, ycm, xbjrk, p_i_fks, xiimax, xinorm &
    & , xi_i_fks, y_ij_fks, xi_i_hat, xpswgt, xjac, pass)
    implicit none
! arguments
    integer icountevts, i_fks, j_fks
    double precision xbjrk_born(2), tau_born, ycm_born, ycmhat, shat_born &
      & , phi_i_fks, xpswgt, xjac, xiimax, xinorm, xp(0:3, nexternal), stot &
      & , x(2), y_ij_fks, xi_i_hat
    double precision shat, sqrtshat, tau, ycm, xbjrk(2), p_i_fks(0:3)
    logical pass
! local
    integer i, j, idir
    double precision yijdir, costh_i_fks, x1bar2, x2bar2, yij_sol, xi1, xi2 &
      & , ximaxtmp, omega, bstfact, shy_tbst, chy_tbst, chy_tbstmo &
      & , xdir_t(3), cosphi_i_fks, sinphi_i_fks, shy_lbst, chy_lbst &
      & , encmso2, E_i_fks, sinth_i_fks, xpifksred(0:3), xi_i_fks &
      & , xiimin, yij_upp, yij_low, y_ij_fks_upp, y_ij_fks_low
    complex(kind=kind(0d0)) resAoR0

    double precision omx1bar2, omx2bar2
    double precision ltau_born, e2ycm_born, em2ycm_born
! external
!
! parameters
    double precision pi
    parameter(pi=3.1415926535897932d0)
    complex(kind=kind(0d0)) ximag
    parameter(ximag=(0d0, 1d0))
    double precision stiny, sstiny, qtiny, zero, ctiny, cctiny
    parameter(stiny=1d-6)
    parameter(qtiny=1d-7)
    parameter(zero=0d0)
    parameter(ctiny=5d-7)
!
    pass = .true.
    if (softtest) then
      sstiny = 0.d0
    else
      sstiny = stiny
    end if
!
! FKS for left or right incoming parton
!
    if (j_fks .eq. 1) then
      idir = 1
    elseif (j_fks .eq. 2) then
      idir = -1
    else
      write (*, *) 'Invalid initial-state FKS sister', j_fks
      stop 1
    end if

    if (1d0 - tau_born .gt. stiny) then
      ltau_born = log(tau_born)
    else
      ltau_born = tau_born - 1d0
    end if
    if (abs(ycm_born) .gt. stiny) then
      e2ycm_born = exp(2*ycm_born)
      em2ycm_born = exp(-2*ycm_born)
    else
      e2ycm_born = 1d0 + 2*ycm_born + 2*ycm_born**2
      em2ycm_born = 1d0 - 2*ycm_born + 2*ycm_born**2
    end if

!
! set-up lower and upper bounds on y_ij_fks
!
    if (tau_born .le. tau_lower_bound .and. ycm_born .gt. &
      & (0.5d0*ltau_born - log(tau_lower_bound))) then
      yij_upp = (tau_lower_bound + tau_born)* &
        & (1 - e2ycm_born*tau_lower_bound)/ &
        & ((tau_lower_bound - tau_born)* &
        & (1 + e2ycm_born*tau_lower_bound))
    else
      yij_upp = 1.d0
    end if
    if (tau_born .le. tau_lower_bound .and. ycm_born .lt. &
      & (-0.5d0*ltau_born + log(tau_lower_bound))) then
      yij_low = -(tau_lower_bound + tau_born)* &
        & (1 - em2ycm_born*tau_lower_bound)/ &
        & ((tau_lower_bound - tau_born)* &
        & (1 + em2ycm_born*tau_lower_bound))
    else
      yij_low = -1.d0
    end if
!
    if (idir .eq. 1) then
      y_ij_fks_upp = yij_upp
      y_ij_fks_low = yij_low
    elseif (idir .eq. -1) then
      y_ij_fks_upp = -yij_low
      y_ij_fks_low = -yij_upp
    end if

!
! set-up y_ij_fks
!
    if (colltest) then
      cctiny = 0.d0
    else
      cctiny = ctiny
    end if
    if ((icountevts .eq. -100 .or. icountevts .eq. 0) .and. &
      & ((.not. softtest) .or. &
      & (softtest .and. y_ij_fks_fix .eq. -2.d0)) .and. &
      & (.not. colltest)) then
! importance sampling towards collinear singularity
! insert here further importance sampling towards y_ij_fks->1
      y_ij_fks = y_ij_fks_upp - &
        & (y_ij_fks_upp - y_ij_fks_low)*(cctiny + (1 - cctiny)*x(2)**2)
    elseif ((icountevts .eq. -100 .or. icountevts .eq. 0) .and. &
      & ((softtest .and. y_ij_fks_fix .ne. -2.d0) .or. &
      & colltest)) then
      y_ij_fks = y_ij_fks_fix
      if (y_ij_fks .gt. y_ij_fks_upp + 1d-12 .or. &
        & y_ij_fks .lt. y_ij_fks_low - 1d-12) then
        xjac = -33d0
        pass = .false.
        return
      end if
    elseif (icountevts .eq. 2 .or. icountevts .eq. 1) then
      y_ij_fks = 1d0
! Check that y_ij_fks is in the allowed range. If not, counter events
! cannot be generated
      if (y_ij_fks .gt. y_ij_fks_upp + 1d-12 .or. &
        & y_ij_fks .lt. y_ij_fks_low - 1d-12) then
        xjac = -33d0
        pass = .false.
        return
      end if
    else
      write (*, *) 'Error #8 in genps_fks.f', icountevts
      stop
    end if
! importance sampling towards collinear singularity
    xjac = xjac*(y_ij_fks_upp - y_ij_fks_low)*x(2)*2d0
!
! Compute costh_i_fks
!
    yijdir = idir*y_ij_fks
    costh_i_fks = yijdir
!
! Compute maximal xi_i_fks
!
    x1bar2 = xbjrk_born(1)**2
    omx1bar2 = 1d0 - x1bar2
    x2bar2 = xbjrk_born(2)**2
    omx2bar2 = 1d0 - x2bar2

    if (1 - tau_born .gt. 1.d-5) then
      yij_sol = -sinh(ycm_born)*(1 + tau_born)/ &
        & (cosh(ycm_born)*(1 - tau_born))
    else
      yij_sol = -ycmhat
    end if
    if (abs(yij_sol) .gt. 1.d0) then
      if (abs(yij_sol) .lt. 1d0 + qtiny) then
        yij_sol = sign(1d0, yij_sol)
      else
        write (*, *) 'Error #9 in genps_fks.f', yij_sol, icountevts
        write (*, *) xbjrk_born(1), xbjrk_born(2), yijdir
      end if
    end if

    if (yijdir .eq. yij_sol) then
      ximaxtmp = 1 - xbjrk_born(1)*xbjrk_born(2)
    elseif (yijdir .ge. yij_sol) then
!this is an expansion when both yij->-1 and x1->1
! in this case there may be precision loosses
! from the argument in the sqrt
      if (abs(yijdir + 1d0) .lt. ctiny .and. omx1bar2 .lt. ctiny) then
        xi1 = (4*x1bar2 + yijdir + 11*x1bar2*yijdir - 5*x1bar2**2*yijdir + &
          & x1bar2**3*yijdir + 4*yijdir**2)/(2*(1 + yijdir)**2)
        ximaxtmp = 1 - xi1
      else if (omx1bar2 .lt. ctiny) then
! compute directly ximaxtmp
        ximaxtmp = omx1bar2/(1 + yijdir)
      else
        xi1 = 2*(1 + yijdir)*x1bar2/( &
          & sqrt(((1 + x1bar2)*(1 - yijdir))**2 + 16*yijdir*x1bar2) + &
          & (1 - yijdir)*(omx1bar2))
        ximaxtmp = 1 - xi1
      end if
    elseif (yijdir .lt. yij_sol) then
!this is an expansion when both yij->+1 and x1->1
! in this case there may be precision loosses
! from the argument in the sqrt
      if (abs(yijdir - 1d0) .lt. ctiny .and. omx2bar2 .lt. ctiny) then
        xi2 = (4*x2bar2 - yijdir - 11*x2bar2*yijdir + 5*x2bar2**2*yijdir - &
          & x2bar2**3*yijdir + 4*yijdir**2)/(4*(-1 + yijdir)**2)
        ximaxtmp = 1 - xi2
      else if (omx2bar2 .lt. ctiny) then
! compute directly ximaxtmp
        ximaxtmp = omx2bar2/(1 - yijdir)
      else
        xi2 = 2*(1 - yijdir)*x2bar2/( &
          & sqrt(((1 + x2bar2)*(1 + yijdir))**2 - 16*yijdir*x2bar2) + &
          & (1 + yijdir)*(omx2bar2))
        ximaxtmp = 1 - xi2
      end if
    else
      write (*, *) 'Fatal error #14 in one_tree: unknown option'
      write (*, *) y_ij_fks, yij_sol, idir
      stop
    end if
    xiimax = ximaxtmp
!
! Lower bound on xi_i_fks
!
    if (tau_born .lt. tau_lower_bound) then
      xiimin = 1d0 - tau_born/tau_lower_bound
    else
      xiimin = 0d0
    end if
    if (xiimax .lt. xiimin) then
      write (*, *) 'WARNING #10 in genps_fks.f', icountevts, xiimax &
        & , xiimin
      xjac = -342d0
      pass = .false.
      return
    end if

    xinorm = xiimax - xiimin
    if (icountevts .ge. 1 .and. &
      & ((idir .eq. 1 .and. &
      & abs(ximaxtmp - (1 - xbjrk_born(1))) .gt. 1.d-5) .or. &
      & (idir .eq. -1 .and. &
      & abs(ximaxtmp - (1 - xbjrk_born(2))) .gt. 1.d-5))) then
      write (*, *) 'Fatal error #15 in one_tree'
      write (*, *) ximaxtmp, xbjrk_born(1), xbjrk_born(2), idir
      stop
    end if
!
! Define xi_i_fks
!
    if ((icountevts .eq. -100 .or. icountevts .eq. 1) .and. &
      & ((.not. colltest) .or. &
      & (colltest .and. xi_i_fks_fix .eq. -2.d0)) .and. &
      & (.not. softtest)) then
      if (icountevts .eq. -100) then
! importance sampling towards soft singularity
! insert here further importance sampling towards xi_i_hat->0
        xi_i_hat = sstiny + (1 - sstiny)*x(1)**2
      end if
      xi_i_fks = xiimin + (xiimax - xiimin)*xi_i_hat
    elseif ((icountevts .eq. -100 .or. icountevts .eq. 1) .and. &
      & (colltest .and. xi_i_fks_fix .ne. -2.d0) .and. &
      & (.not. softtest)) then
      if (xi_i_fks_fix .lt. xiimax) then
        xi_i_fks = xi_i_fks_fix
      else
        xi_i_fks = xi_i_fks_fix*xiimax
      end if
    elseif ((icountevts .eq. -100 .or. icountevts .eq. 1) .and. &
      & softtest) then
      if (xi_i_fks_fix .lt. xiimax) then
        xi_i_fks = xi_i_fks_fix
      else
        xjac = -102
        pass = .false.
        return
      end if
    elseif (icountevts .eq. 2 .or. icountevts .eq. 0) then
      xi_i_fks = 0d0
! Check that xi_i_fks is in the allowed range. If not, counter events
! cannot be generated
      if (xi_i_fks .gt. xiimax + 1d-12 .or. &
        & xi_i_fks .lt. xiimin - 1d-12) then
        xjac = -34d0
        pass = .false.
        return
      end if
    else
      write (*, *) 'Error #11 in genps_fks.f', icountevts
      stop
    end if
! remove the following if no importance sampling towards soft
! singularity is performed when integrating over xi_i_hat
    xjac = xjac*2d0*x(1)
!
! Initial state variables are different for events and counterevents. Update them here.
!
    omega = sqrt((2 - xi_i_fks*(1 + yijdir))/ &
      & (2 - xi_i_fks*(1 - yijdir)))
    if (icountevts .ne. 0) then
      tau = tau_born/(1 - xi_i_fks)
      ycm = ycm_born - log(omega)
      shat = tau*stot
      sqrtshat = sqrt(shat)
      xbjrk(1) = xbjrk_born(1)/(sqrt(1 - xi_i_fks)*omega)
      xbjrk(2) = xbjrk_born(2)*omega/sqrt(1 - xi_i_fks)
    else
      tau = tau_born
      ycm = ycm_born
      shat = shat_born
      sqrtshat = sqrt(shat)
      xbjrk(1) = xbjrk_born(1)
      xbjrk(2) = xbjrk_born(2)
    end if
!
! Define the boost factor here
!
    bstfact = sqrt((2 - xi_i_fks*(1 - yijdir))*(2 - xi_i_fks*(1 + yijdir)))
    shy_tbst = -xi_i_fks*sqrt(1 - yijdir**2)/(2*sqrt(1 - xi_i_fks))
    chy_tbst = bstfact/(2*sqrt(1 - xi_i_fks))
    chy_tbstmo = chy_tbst - 1.d0
    cosphi_i_fks = cos(phi_i_fks)
    sinphi_i_fks = sin(phi_i_fks)
    xdir_t(1) = -cosphi_i_fks
    xdir_t(2) = -sinphi_i_fks
    xdir_t(3) = zero
!
    shy_lbst = -xi_i_fks*yijdir/bstfact
    chy_lbst = (2 - xi_i_fks)/bstfact
! Boost the momenta
    do i = 3, nexternal
      if (i .ne. i_fks .and. shy_tbst .ne. 0.d0) &
        & call boostwdir2_in_place(chy_tbst, shy_tbst, chy_tbstmo, &
        & xdir_t, xp(0, i))
    end do
!
    encmso2 = sqrtshat/2.d0
    p_i_fks(0) = encmso2
    E_i_fks = xi_i_fks*encmso2
    sinth_i_fks = sqrt(1 - costh_i_fks**2)
!
    xp(0, 1) = encmso2*(chy_lbst - shy_lbst)
    xp(1, 1) = 0.d0
    xp(2, 1) = 0.d0
    xp(3, 1) = xp(0, 1)
!
    xp(0, 2) = encmso2*(chy_lbst + shy_lbst)
    xp(1, 2) = 0.d0
    xp(2, 2) = 0.d0
    xp(3, 2) = -xp(0, 2)
!
    xp(0, i_fks) = E_i_fks*(chy_lbst - shy_lbst*yijdir)
    p_i_fks(0) = p_i_fks(0)*(chy_lbst - shy_lbst*yijdir)
    xpifksred(1) = sinth_i_fks*cosphi_i_fks
    xpifksred(2) = sinth_i_fks*sinphi_i_fks
    xpifksred(3) = chy_lbst*yijdir - shy_lbst
!
    do j = 1, 3
      xp(j, i_fks) = E_i_fks*xpifksred(j)
      p_i_fks(j) = encmso2*xpifksred(j)
    end do
!
! Collinear limit of <ij>/[ij]. See innerpin.m.
    if (icountevts .eq. -100 .or. &
      & (icountevts .eq. 1 .and. xij_aor .eq. 0)) then
      resAoR0 = -exp(2*idir*ximag*phi_i_fks)
      xij_aor = resAoR0
    end if
!
! Phase-space factor for (xii,yij,phii) * (tau,ycm)
    xpswgt = xpswgt*shat
    xpswgt = xpswgt/(4*pi)**3/(1 - xi_i_fks)
    xpswgt = abs(xpswgt)
!
    return
  end subroutine generate_momenta_initial

  subroutine getangles(pin, th, cth, sth, phi, cphi, sphi)
    implicit none
    double precision pin(0:3), th, cth, sth, phi, cphi, sphi, xlength
!
    xlength = pin(1)**2 + pin(2)**2 + pin(3)**2
    if (xlength .eq. 0) then
      th = 0.d0
      cth = 1.d0
      sth = 0.d0
      phi = 0.d0
      cphi = 1.d0
      sphi = 0.d0
    else
      xlength = sqrt(xlength)
      cth = pin(3)/xlength
      th = acos(cth)
      if (cth .ne. 1.d0) then
        sth = sqrt(1 - cth**2)
        phi = atan2(pin(2), pin(1))
        cphi = cos(phi)
        sphi = sin(phi)
      else
        sth = 0.d0
        phi = 0.d0
        cphi = 1.d0
        sphi = 0.d0
      end if
    end if
    return
  end subroutine getangles

  subroutine get_recoil(p_born, imother, shat_born, xmrec2, pass)
    implicit none
    double precision p_born(0:3, nexternal - 1), xmrec2, shat_born
    logical pass
    integer imother, i
    double precision recoilbar(0:3)
    pass = .true.
    do i = 0, 3
      if (nincoming .eq. 2) then
        recoilbar(i) = p_born(i, 1) + p_born(i, 2) - p_born(i, imother)
      else
        recoilbar(i) = p_born(i, 1) - p_born(i, imother)
      end if
    end do
    xmrec2 = dot(recoilbar, recoilbar)
    if (xmrec2 .lt. 0.d0) then
      if (abs(xmrec2) .gt. (1.d-4*shat_born)) then
        write (*, *) 'Fatal error #14 in genps_fks.f', xmrec2, imother
        stop
      else
        write (*, *) 'Error #15 in genps_fks.f', xmrec2, imother
        pass = .false.
        return
      end if
    end if
    if (xmrec2 .ne. xmrec2) then
      write (*, *) 'Error #16 in setting up event in genps_fks.f,'// &
        & ' skipping event'
      pass = .false.
      return
    end if
    return
  end subroutine get_recoil

  subroutine validate_bound_genps_state()
    implicit none

    if (size(cnt_momenta, 1) /= 4 .or. &
        size(cnt_momenta, 2) /= nexternal .or. &
        size(cnt_momenta, 3) /= 3) then
      call fail_genps_state('counterevent momenta have inconsistent bounds')
    end if
    if (size(cnt_weight) /= 3 .or. size(cnt_psweight) /= 3 .or. &
        size(cnt_jacobian) /= 3) then
      call fail_genps_state('counterevent weights have inconsistent bounds')
    end if
    if (size(born_lab_momenta, 1) /= 4 .or. &
        size(born_lab_momenta, 2) /= nexternal - 1 .or. &
        any(shape(born_coll_momenta) /= shape(born_lab_momenta)) .or. &
        any(shape(born_norad_momenta) /= shape(born_lab_momenta))) then
      call fail_genps_state('Born radiation momenta have inconsistent bounds')
    end if
    if (size(event_momenta, 1) /= 4 .or. size(event_momenta, 2) /= nexternal) then
      call fail_genps_state('event momenta have inconsistent bounds')
    end if
    if (size(particle_masses) /= nexternal) then
      call fail_genps_state('particle masses have inconsistent bounds')
    end if
  end subroutine validate_bound_genps_state

  subroutine require_genps_state()
    implicit none

    if (.not. genps_state_initialized) then
      call fail_genps_state('module state has not been initialized')
    end if
    if (.not. associated(cnt_momenta) .or. &
        .not. associated(cnt_weight) .or. &
        .not. associated(cnt_psweight) .or. &
        .not. associated(cnt_jacobian) .or. &
        .not. associated(born_lab_momenta) .or. &
        .not. associated(born_coll_momenta) .or. &
        .not. associated(born_norad_momenta) .or. &
        .not. associated(event_momenta) .or. &
        .not. associated(particle_masses)) then
      call fail_genps_state('module state is incomplete')
    end if
  end subroutine require_genps_state

  subroutine fail_genps_state(message)
    implicit none
    character(len=*), intent(in) :: message

    write (*, *) 'genps_fks: ', trim(message)
    stop 1
  end subroutine fail_genps_state

end module genps_fks
