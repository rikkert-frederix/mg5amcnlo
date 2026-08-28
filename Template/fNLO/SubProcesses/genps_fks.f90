module genps_fks
  use boostwdir2_module, only: boostwdir2_in_place
  use process_dimensions, only: validate_process_dimensions
  use run_state
  use timing_state, only: tGenPS
  use kin_functions_module, only: dot => dot_impl, rho => rho_impl
  use phase_space_kinematics, only: rotate_invar, phspncheck_nocms, &
                                    phase_space_lambda
  use fks_diagnostics, only: xmom_compare
  use genps_born, only: born_phase_space, generate_born_phase_space, &
                        invalidate_born_phase_space
  use factorized_phase_space, only: factorized_radiation_state, &
       factorized_measure_state, &
       reset_factorized_phase_space, store_factorized_kernel_momenta, &
       store_factorized_radiation_state, &
       scale_factorized_radiation_jacobians, &
       store_factorized_event_measure, &
       multiply_factorized_event_measure, &
       compose_factorized_event_measure, &
       scale_factorized_event_measures
  use decay_chain_metadata, only: has_decay_chains, context_for_fks, &
                                  real_phase_space_dimension
  use decay_chain_kinematics, only: get_core_born_momenta, &
       get_core_mass_buffer, active_core_count, expand_real_decay_momenta, &
       store_core_event_momenta, fks_leg_mass
  use nlo_decay_metadata, only: has_nlo_decay, nlo_decay_corrected_node
  use nlo_decay_kinematics, only: generate_nlo_decay_fks_event, &
                                  nlo_decay_fks_sister_mass, &
                                  nlo_decay_parent_mass, &
                                  fill_nlo_decay_event_masses
  ! Generated parameters keep the phase-space normalization bit-identical
  ! to the NLO template's compile-time arithmetic.
  use fnlo_process_common, only: nexternal, nincoming, max_particles, &
                                 max_branch, nfksprocess, nocntevents, nbody, &
                                 soft_counterevent, collinear_counterevent, &
                                 soft_collinear_counterevent, real_event, &
                                 first_counterevent, last_counterevent, &
                                 xij_aor, i_fks, j_fks, &
                                 event_xi, event_y, event_xi_hat, &
                                 event_fks_momentum, event_xi_max, &
                                 event_xi_norm, event_bjorken_x, &
                                 event_sqrt_shat, event_shat, &
                                 ybst_til_tolab, ybst_til_tocm, &
                                 softtest, colltest, xi_i_fks_fix, &
                                 y_ij_fks_fix, tau_lower_bound, &
                                 stored_event_momenta => event_momenta, &
                                 stored_event_jacobian => event_jacobian, &
                                 born_lab_momenta => p_born_l, &
                                 particle_masses
  implicit none
  private
  public :: generate_momenta

  double precision :: massive_xjac_cache = 1d0

contains

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
    call validate_process_dimensions()
    call cpu_time(tBefore)

    ! Block-local momenta belong to this one event/counterevent family.  Each
    ! block is generated in its own rest frame and boosted independently;
    ! matrix-element providers consume that cache, not the visible event.
    call reset_factorized_phase_space()
    call generate_born_phase_space(ndim, iconfig, x, born)

    pass = born%valid
    if (pass) then
      if (has_nlo_decay()) then
        call generate_nlo_decay_FKS_kinematics(x, ndim, born%shat, &
             born%sqrtshat, born%ycm, born%xbjrk, &
             jac, p, pass)
      else
        call generate_FKS_kinematics(x, ndim, born%xjac, born%xpswgt, &
                                     born%stot, born%shat, born%sqrtshat, &
                                     born%tau, born%ycm, born%ycmhat, &
                                     born%xbjrk, born%masses, &
                                     born%external_masses, jac, p, pass)
      end if
    end if

    if (.not. pass .or. born%xjac < 0d0) then
      jac = -222d0
      stored_event_jacobian = -222d0
      p(0, 1) = -99d0
      do i = first_counterevent, real_event
        stored_event_momenta(0, 1, i) = -99d0
      end do
      call invalidate_born_phase_space()
      nocntevents = .true.
    end if

! Apply any incoming integration weight to the event and counterevents.
    do i = first_counterevent, real_event
      stored_event_jacobian(i) = stored_event_jacobian(i)*wgt
    end do
    if (has_nlo_decay() .or. has_decay_chains()) then
      call scale_factorized_event_measures(wgt)
    end if
    call scale_factorized_radiation_jacobians(wgt)
    wgt = stored_event_jacobian(real_event)

    call cpu_time(tAfter)
    tGenPS = tGenPS + (tAfter - tBefore)
  end subroutine generate_momenta


  subroutine generate_nlo_decay_FKS_kinematics(x, ndim, &
       shat_born, sqrtshat_born, ycm_born, xbjrk_born, jac, p, pass)
    implicit none
    integer, intent(in) :: ndim
    double precision, intent(in) :: x(99)
    double precision, intent(in) :: shat_born, sqrtshat_born, ycm_born
    double precision, intent(in) :: xbjrk_born(2)
    double precision, intent(out) :: jac, p(0:3, nexternal)
    logical, intent(out) :: pass

    integer :: event_position, event_slot, solution_sign
    integer :: event_generation_order(4)
    double precision :: visible(0:3, nexternal), event_masses(-max_branch:max_particles)
    double precision :: xiimax, xinorm, xi_i, xi_hat, y_ij
    double precision :: p_i_hat(0:3), kernel_p_i_hat(0:3), xjac, xpswgt
    double precision :: radiation_jacobian, radiation_weight
    double precision :: decay_mass
    logical :: event_pass, real_pass, massive_sister, skip_counterevents

    event_generation_order = (/real_event, soft_counterevent, &
         collinear_counterevent, soft_collinear_counterevent/)
    pass = .false.
    real_pass = .false.
    skip_counterevents = .false.
    massive_sister = nlo_decay_fks_sister_mass(nfksprocess) > 0d0
    xi_hat = 0d0
    xij_aor = (0d0, 0d0)
    ybst_til_tolab = -ycm_born - 0.5d0*log(ebeam(1)/ebeam(2))
    ybst_til_tocm = 0d0
    event_masses = 0d0
    decay_mass = nlo_decay_parent_mass()
    call fill_nlo_decay_event_masses(nfksprocess, &
                                     event_masses(1:nexternal))
    do event_slot = first_counterevent, real_event
      event_fks_momentum(0, event_slot) = -1d0
      event_xi_max(event_slot) = -1d0
      stored_event_jacobian(event_slot) = -1d0
      stored_event_momenta(:, :, event_slot) = 0d0
      stored_event_momenta(0, 1, event_slot) = -99d0
      event_shat(event_slot) = shat_born
      event_sqrt_shat(event_slot) = sqrtshat_born
      event_bjorken_x(:, event_slot) = xbjrk_born
    end do

    do event_position = 1, size(event_generation_order)
      event_slot = event_generation_order(event_position)
      if (skip_counterevents) exit
      if (massive_sister .and. &
          (event_slot == collinear_counterevent .or. &
           event_slot == soft_collinear_counterevent)) cycle

      xjac = 1d0
      xpswgt = 1d0
      xiimax = -1d0
      xinorm = -1d0
      xi_i = -1d0
      y_ij = -2d0
      p_i_hat = -1d0
      kernel_p_i_hat = -1d0
      solution_sign = 1
      call generate_nlo_decay_fks_event(nfksprocess, event_slot, x, ndim, &
           xjac, xpswgt, visible, xiimax, xinorm, xi_i, xi_hat, y_ij, &
           p_i_hat, kernel_p_i_hat, solution_sign, event_pass)
      if (event_pass .and. xjac > 0d0) then
        call phspncheck_nocms(nexternal, sqrtshat_born, event_masses, &
                              visible, event_pass)
        if (.not. event_pass) xjac = -197d0
      end if
      if (event_pass .and. xjac > 0d0) then
        radiation_jacobian = xjac
        radiation_weight = xpswgt
        call finalize_factorized_event_measure( &
             event_slot, nlo_decay_corrected_node(), &
             radiation_jacobian, radiation_weight, shat_born, &
             sqrtshat_born, event_masses(1), event_masses(2), xjac)
      else
        xjac = -196d0
        radiation_jacobian = -1d0
        radiation_weight = -1d0
      end if
      call store_FKS_event(event_slot, xiimax, xinorm, xi_i, xi_hat, &
                           p_i_hat, y_ij, visible, p, xjac, jac, &
                           nlo_decay_corrected_node(), kernel_p_i_hat, &
                           decay_mass**2, decay_mass, 0d0, &
                           radiation_jacobian, radiation_weight)

      if (event_slot == real_event) then
        real_pass = event_pass .and. xjac > 0d0
        if (solution_sign == -1) skip_counterevents = .true.
      end if
    end do

    if (skip_counterevents) then
      do event_slot = first_counterevent, last_counterevent
        stored_event_jacobian(event_slot) = -299d0
        stored_event_momenta(0, 1, event_slot) = -99d0
      end do
    end if
    nocntevents = all(stored_event_jacobian( &
         first_counterevent:last_counterevent) <= 0d0)
    pass = real_pass
  end subroutine generate_nlo_decay_FKS_kinematics

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

    integer icountevts, event_position
    integer ixEi, ixpi, imother
    double precision xmrec2, m_j_fks, phi_i_fks, rat_xi, &
      & xi_i_fks, y_ij_fks, xi_i_hat, xiimax, xinorm, xjac, xpswgt, &
      & xp(0:3, nexternal), visible_xp(0:3, nexternal), p_i_fks(0:3)
    double precision :: radiation_jacobian, radiation_weight
    double precision :: core_born(0:3, nexternal - 1)
    double precision :: core_mass_values(nexternal)
    integer i
    integer :: core_particle_count
    integer, parameter :: event_generation_order(4) = (/ &
      real_event, soft_counterevent, collinear_counterevent, &
      soft_collinear_counterevent /)

    double precision pi
    parameter(pi=3.1415926535897932d0)

    integer isolsign
    logical skip_counterevents

! check that the starting PS point is meaningful
    if (xjac0 .lt. 0d0) then
      pass = .false.
      return
    end if
!
! Here we start with the FKS Stuff
!
! The real event occupies slot 3; slots 0, 1, and 2 are its soft,
! collinear, and soft-collinear counterevents. They are generated in
! that explicit order because the real point determines whether the
! massive-emission second solution has counterevents.
    skip_counterevents = .false.
! Event/counterevent values stay negative until their ordinary context is
! generated.  This also forces a failure if a skipped context is consumed.
    do i = first_counterevent, real_event
      event_fks_momentum(0, i) = -1.d0
      event_xi_max(i) = -1.d0
      stored_event_jacobian(i) = -1.d0
    end do
! All event momenta use the same tilde-to-laboratory boost. The
! tilde-to-partonic-CM boost depends on the event and is filled below.
    ybst_til_tolab = -ycm_born - 0.5d0*log(ebeam(1)/ebeam(2))
    ybst_til_tocm = 1d14
! The soft counterevent is always at the Born rapidity, including when
! cuts prevent its ordinary construction.
    ybst_til_tocm(soft_counterevent) = 0d0
! If the collinear counterevent is not generated, this stays zero.
    xij_aor = (0.d0, 0.d0)
!
! These will correspond to the vegas x's for the FKS variables xi_i,
! y_ij and phi_i (changing this also requires changing folding parameters)
    ixEi = ndim - 2
    ixpi = ndim
!
    imother = min(j_fks, i_fks)
    if (has_decay_chains()) then
      core_particle_count = active_core_count(nfksprocess)
      call get_core_born_momenta(core_born)
      call get_core_mass_buffer(context_for_fks(nfksprocess), &
                                core_mass_values)
      m_j_fks = fks_leg_mass(nfksprocess, j_fks)
    else
      core_particle_count = nexternal
      core_born = born_lab_momenta
      m_j_fks = particle_masses(j_fks)
    end if
!
! For final state j_fks, compute the recoil invariant mass
    if (j_fks .gt. nincoming) then
      call get_recoil(core_born, imother, shat_born, xmrec2, pass)
      if (.not. pass) then
        xjac0 = -44
        return
      end if
    end if

! Here is the beginning of the loop over the momenta for the event and
! counter-events. This will fill the xp momenta with the event and
! counter-event momenta.
    counterevent_loop: do event_position = 1, size(event_generation_order)
      icountevts = event_generation_order(event_position)
! The second solution of the massive mapping has no counterevents.
      if (skip_counterevents) exit counterevent_loop
! Massive emitters and n-body integrations need only the soft
! counterevent in addition to the real point.
      if (icountevts .ne. real_event .and. &
          icountevts .gt. soft_counterevent .and. &
          (m_j_fks .ne. 0d0 .or. nbody)) cycle counterevent_loop

      if (has_decay_chains()) then
        xjac = 1d0
        xpswgt = 1d0
      else
        xjac = xjac0
        xpswgt = xpswgt0
      end if
!
! Put the Born momenta in the xp momenta, making sure that the mapping
! is correct; put i_fks momenta equal to zero.
      xp = 0d0
      if (i_fks .gt. 1) xp(:, 1:i_fks - 1) = core_born(:, 1:i_fks - 1)
      xp(:, i_fks) = 0d0
      if (i_fks .lt. core_particle_count) &
        xp(:, i_fks + 1:core_particle_count) = &
             core_born(:, i_fks:core_particle_count - 1)
      if (has_decay_chains()) then
        m = 0d0
        m(1:core_particle_count) = core_mass_values(1:core_particle_count)
      else
        m(i_fks) = 0d0
        if (i_fks .gt. 1) m(1:i_fks - 1) = m_born(1:i_fks - 1)
        if (i_fks .lt. core_particle_count) &
          m(i_fks + 1:core_particle_count) = &
               m_born(i_fks:core_particle_count - 1)
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
        event_shat(icountevts) = shat_born
        event_sqrt_shat(icountevts) = sqrtshat_born
        ybst_til_tocm(icountevts) = 0d0
        event_bjorken_x(:, icountevts) = xbjrk_born
        if (m_j_fks .eq. 0d0) then
          isolsign = 1
          call generate_momenta_massless_final(icountevts, i_fks, j_fks &
            & , core_born(0:3, imother), event_shat(icountevts) &
            & , event_sqrt_shat(icountevts), x(ixEi), xmrec2, xp &
            & , phi_i_fks, xiimax, xinorm, xi_i_fks, y_ij_fks, xi_i_hat &
            & , p_i_fks, xjac, xpswgt, core_particle_count, pass)
        elseif (m_j_fks .gt. 0d0) then
          call generate_momenta_massive_final(icountevts, isolsign &
            & , rat_xi, i_fks, j_fks, core_born(0:3, imother) &
            & , event_shat(icountevts), event_sqrt_shat(icountevts) &
            & , m_j_fks, x(ixEi), xmrec2, xp, phi_i_fks &
            & , xiimax, xinorm, xi_i_fks, y_ij_fks, xi_i_hat, p_i_fks, xjac &
            & , xpswgt, core_particle_count, pass)
        end if
      elseif (j_fks .le. nincoming) then
        isolsign = 1
        call generate_momenta_initial(icountevts, i_fks, j_fks, xbjrk_born &
          & , tau_born, ycm_born, ycmhat, shat_born, phi_i_fks, xp, x(ixEi) &
          & , event_shat(icountevts), stot, event_sqrt_shat(icountevts) &
          & , ybst_til_tocm(icountevts), event_bjorken_x(:, icountevts) &
          & , p_i_fks, xiimax, xinorm &
          & , xi_i_fks, y_ij_fks, xi_i_hat, xpswgt, xjac, &
          & core_particle_count, pass)
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
      if (pass .and. xjac .gt. 0.d0) then
        call phspncheck_nocms(core_particle_count, &
                              event_sqrt_shat(icountevts), m, xp, pass)
        if (.not. pass) xjac = -199
      end if
      if (pass) then
        radiation_jacobian = xjac
        radiation_weight = xpswgt
        if (has_decay_chains()) then
          call finalize_factorized_event_measure( &
               icountevts, 0, radiation_jacobian, radiation_weight, &
               event_shat(icountevts), event_sqrt_shat(icountevts), &
               m(1), m(2), xjac)
        else
          call compute_flux(event_shat(icountevts), &
                            event_sqrt_shat(icountevts), &
                            m(1), m(2), xpswgt, xjac)
        end if
      else
        radiation_jacobian = -1d0
        radiation_weight = -1d0
      end if
!
      if (has_decay_chains() .and. pass .and. xjac > 0d0) then
        call store_core_event_momenta(icountevts, xp, core_particle_count)
        call expand_real_decay_momenta(nfksprocess, icountevts, x, xp, &
                                       visible_xp, pass)
        if (.not. pass) xjac = -198d0
      else
        if (pass .and. xjac > 0d0) then
          call store_factorized_kernel_momenta( &
               icountevts, 0, core_particle_count, xp)
        end if
        visible_xp = xp
      end if
      call store_FKS_event(icountevts, xiimax, xinorm, &
        & xi_i_fks, xi_i_hat, p_i_fks, y_ij_fks, visible_xp, p, xjac, &
        & jac, 0, p_i_fks, event_shat(icountevts), &
        & event_sqrt_shat(icountevts), ybst_til_tocm(icountevts), &
        & radiation_jacobian, radiation_weight)
      if (icountevts .eq. real_event .and. isolsign .eq. -1) &
        skip_counterevents = .true.
    end do counterevent_loop

    if (skip_counterevents) then
! When integrating over the second fold with massive j_fks, the
! counterevents have been skipped, so make sure their
! momenta are unphysical. Born are physical if event was generated, and
! must stay so for the computation of enhancement factors.
      do i = first_counterevent, last_counterevent
        stored_event_jacobian(i) = -299
        stored_event_momenta(0, 1, i) = -99
      end do
    end if

    nocntevents = all(stored_event_jacobian( &
                        first_counterevent:last_counterevent) .le. 0.d0)
    if (.not. has_decay_chains()) then
      call xmom_compare(i_fks, j_fks, &
                        stored_event_jacobian(first_counterevent:last_counterevent), &
                        stored_event_momenta(:, :, first_counterevent:last_counterevent), &
                        particle_masses, pass)
    end if
!
    return
  end subroutine generate_FKS_kinematics

  subroutine compute_flux(shat, sqrtshat, m1, m2, xpswgt, xjac)
    implicit none
    double precision shat, sqrtshat, m1, m2, xpswgt, xjac

    double precision pwgt, flux

    flux = phase_space_flux_factor(shat, sqrtshat, m1, m2)
    pwgt = max(xjac*xpswgt, 1d-99)
    xjac = pwgt*flux
!
    return
  end subroutine compute_flux


  double precision function phase_space_flux_factor( &
       shat, sqrtshat, m1, m2)
    double precision, intent(in) :: shat, sqrtshat, m1, m2
    double precision, parameter :: pi = 3.1415926535897932d0

    if (nincoming == 2) then
      phase_space_flux_factor = &
           1d0/(2d0*sqrt(phase_space_lambda(shat, m1**2, m2**2)))
    else
      phase_space_flux_factor = 1d0/(2d0*sqrtshat)
    end if
    ! Convert Kajantie's R_n normalization to the complete phase space,
    ! retaining the three pi factors already present in the FKS map.
    phase_space_flux_factor = phase_space_flux_factor/ &
         (2d0*pi)**real_phase_space_dimension()*(2d0*pi)**3
  end function phase_space_flux_factor


  subroutine finalize_factorized_event_measure(event_slot, radiation_block, &
       radiation_jacobian, radiation_weight, shat, sqrtshat, m1, m2, xjac)
    integer, intent(in) :: event_slot, radiation_block
    double precision, intent(in) :: radiation_jacobian, radiation_weight
    double precision, intent(in) :: shat, sqrtshat, m1, m2
    double precision, intent(out) :: xjac
    type(factorized_measure_state) :: radiation_measure, flux_measure
    double precision :: composed_jacobian, composed_weight, flux
    logical :: available

    radiation_measure%jacobian = radiation_jacobian
    radiation_measure%phase_space_weight = radiation_weight
    call store_factorized_event_measure( &
         event_slot, radiation_block, radiation_measure)
    call compose_factorized_event_measure( &
         event_slot, composed_jacobian, composed_weight, available)
    if (.not. available) then
      write (*, '(a)') &
           'ERROR in genps_fks: a factorized base measure is unavailable'
      stop 1
    end if

    flux = phase_space_flux_factor(shat, sqrtshat, m1, m2)
    xjac = max(composed_jacobian*composed_weight, 1d-99)*flux
    flux_measure = factorized_measure_state()
    flux_measure%phase_space_weight = flux
    call multiply_factorized_event_measure(event_slot, 0, flux_measure)
  end subroutine finalize_factorized_event_measure

  subroutine store_FKS_event(icountevts, xiimax, xinorm, &
    & xi_i_fks, xi_i_hat, p_i_fks, y_ij_fks, xp, p, xjac, jac, &
    & radiation_block, kernel_p_i_fks, kernel_shat, kernel_sqrt_shat, &
    & kernel_y_to_cm, radiation_jacobian, radiation_weight)

    implicit none
    integer icountevts, radiation_block
    double precision xiimax, xinorm, &
      & xi_i_fks, xi_i_hat, p_i_fks(0:3), y_ij_fks, xp(0:3, nexternal), p(0:3, nexternal), &
      & xjac, jac, kernel_p_i_fks(0:3), kernel_shat, kernel_sqrt_shat, &
      & kernel_y_to_cm, radiation_jacobian, radiation_weight

    integer i, j
    type(factorized_radiation_state) :: radiation

! Catch the points for which there is no viable phase-space generation
! (still fill the shared state with some information that is needed
! (e.g. the counterevent rapidity)).
    if (xjac .le. 0d0) then
      xp(0, 1) = -99d0
    end if
!
! Fill shared FKS state.
    if (icountevts < soft_counterevent .or. icountevts > real_event) then
      write (*, *) 'Invalid event index in store_FKS_event:', icountevts
      stop 1
    end if
    event_xi_max(icountevts) = xiimax
    event_xi_norm(icountevts) = xinorm
    event_xi(icountevts) = xi_i_fks
    event_xi_hat(icountevts) = xi_i_hat
    event_y(icountevts) = y_ij_fks
    do i = 0, 3
      event_fks_momentum(i, icountevts) = p_i_fks(i)
    end do
    do i = 1, nexternal
      do j = 0, 3
        stored_event_momenta(j, i, icountevts) = xp(j, i)
      end do
    end do
    stored_event_jacobian(icountevts) = xjac

    radiation%jacobian = xjac
    radiation%radiation_jacobian = radiation_jacobian
    radiation%radiation_weight = radiation_weight
    radiation%xi = xi_i_fks
    radiation%y = y_ij_fks
    radiation%xi_hat = xi_i_hat
    radiation%xi_max = xiimax
    radiation%xi_norm = xinorm
    radiation%fks_momentum = kernel_p_i_fks
    radiation%shat = kernel_shat
    radiation%sqrt_shat = kernel_sqrt_shat
    radiation%y_to_cm = kernel_y_to_cm
    call store_factorized_radiation_state( &
         icountevts, radiation_block, radiation)

    if (icountevts .eq. real_event) then
      do i = 1, nexternal
        do j = 0, 3
          p(j, i) = xp(j, i)
        end do
      end do
      jac = xjac
    end if

    return
  end subroutine store_FKS_event

  subroutine generate_momenta_massless_final(icountevts, i_fks, j_fks &
    & , p_born_imother, shat, sqrtshat, x, xmrec2, xp, phi_i_fks, xiimax &
    & , xinorm, xi_i_fks, y_ij_fks, xi_i_hat, p_i_fks, xjac, xpswgt &
    & , particle_count, pass)
    implicit none
! arguments
    integer icountevts, i_fks, j_fks, particle_count
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
    if ((icountevts .eq. real_event .or. &
         icountevts .eq. soft_counterevent) .and. &
      & ((.not. softtest) .or. &
      & (softtest .and. y_ij_fks_fix .eq. -2.d0)) .and. &
      & (.not. colltest)) then
! importance sampling towards collinear singularity
! insert here further importance sampling towards y_ij_fks->1
      y_ij_fks = -2d0*(cctiny + (1 - cctiny)*x(2)**2) + 1d0
    elseif ((icountevts .eq. real_event .or. &
             icountevts .eq. soft_counterevent) .and. &
      & ((softtest .and. y_ij_fks_fix .ne. -2.d0) .or. &
      & colltest)) then
      y_ij_fks = y_ij_fks_fix
    elseif (icountevts .eq. collinear_counterevent .or. &
            icountevts .eq. soft_collinear_counterevent) then
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
    if ((icountevts .eq. real_event .or. &
         icountevts .eq. collinear_counterevent) .and. &
      & ((.not. colltest) .or. &
      & (colltest .and. xi_i_fks_fix .eq. -2.d0)) .and. &
      & (.not. softtest)) then
      if (icountevts .eq. real_event) then
! importance sampling towards soft singularity
! insert here further importance sampling towards xi_i_hat->0
        xi_i_hat = sstiny + (1 - sstiny)*x(1)**2
      end if
! in the case of counter events, xi_i_hat is an input to this function
      xi_i_fks = xi_i_hat*xiimax
    elseif ((icountevts .eq. real_event .or. &
             icountevts .eq. collinear_counterevent) .and. &
      & (colltest .and. xi_i_fks_fix .ne. -2.d0) .and. &
      & (.not. softtest)) then
! This is to keep xi_i_hat, rather than xi_i, fixed in the tests.
      if (xi_i_fks_fix .lt. xiimax) then
        xi_i_fks = xi_i_fks_fix*xiimax
      else
        xi_i_fks = xi_i_fks_fix*xiimax
      end if
    elseif ((icountevts .eq. real_event .or. &
             icountevts .eq. collinear_counterevent) .and. &
      & softtest) then
      if (xi_i_fks_fix .lt. 1d0) then
        xi_i_fks = xi_i_fks_fix*xiimax
      else
        xjac = -102
        pass = .false.
        return
      end if
    elseif (icountevts .eq. soft_collinear_counterevent .or. &
            icountevts .eq. soft_counterevent) then
      xi_i_fks = 0d0
    else
      write (*, *) 'Error #4 in genps_fks.f', icountevts
      stop
    end if
! remove the following if no importance sampling towards soft
! singularity is performed when integrating over xi_i_hat
    xjac = xjac*2d0*x(1)

! Check that xii is in the allowed range
    if (icountevts .eq. real_event .or. &
        icountevts .eq. collinear_counterevent) then
      if (xi_i_fks .gt. (1 - xmrec2/shat)) then
        xjac = -101
        pass = .false.
        return
      end if
    elseif (icountevts .eq. soft_counterevent .or. &
            icountevts .eq. soft_collinear_counterevent) then
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
    do i = nincoming + 1, particle_count
      if (i .ne. i_fks .and. i .ne. j_fks .and. shybst .ne. 0.d0) &
        & call boostwdir2_in_place(chybst, shybst, chybstmo, xdir, &
        & xp(0, i))
    end do
!
! Collinear limit of <ij>/[ij]. See innerp3.m.
    if ((icountevts .eq. real_event .or. &
      & (icountevts .eq. collinear_counterevent .and. &
         xij_aor .eq. 0))) then
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
    xpswgt = xpswgt*2*shat/(4*pi)**3*veckn/veckbarn/ &
      & (2 - xi_i_fks*(1 - xp(0, j_fks)/veckn*y_ij_fks))
    xpswgt = abs(xpswgt)
    return
  end subroutine generate_momenta_massless_final

  subroutine generate_momenta_massive_final(icountevts, isolsign &
    & , rat_xi, i_fks, j_fks, p_born_imother, shat &
    & , sqrtshat, m_j_fks, x, xmrec2, xp, phi_i_fks, xiimax, xinorm &
    & , xi_i_fks, y_ij_fks, xi_i_hat, p_i_fks, xjac, xpswgt, &
    & particle_count, pass)
    implicit none
! arguments
    integer icountevts, i_fks, j_fks, isolsign, particle_count
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
      & icountevts .eq. collinear_counterevent .or. &
      & icountevts .eq. soft_collinear_counterevent) then
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
    if ((icountevts .eq. real_event .or. &
         icountevts .eq. soft_counterevent) .and. &
      & ((.not. softtest) .or. &
      & (softtest .and. y_ij_fks_fix .eq. -2.d0)) .and. &
      & (.not. colltest)) then
! importance sampling towards collinear singularity
! insert here further importance sampling towards y_ij_fks->1
      y_ij_fks = -2d0*(cctiny + (1 - cctiny)*x(2)**2) + 1d0
    elseif ((icountevts .eq. real_event .or. &
             icountevts .eq. soft_counterevent) .and. &
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
    if (icountevts .eq. real_event .and. &
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
    elseif (icountevts .eq. real_event .and. &
      & (colltest .and. xi_i_fks_fix .ne. -2.d0) .and. &
      & (.not. softtest)) then
      massive_xjac_cache = 1.d0
      if (xi_i_fks_fix .lt. xiimax) then
        xi_i_fks = xi_i_fks_fix
      else
        xi_i_fks = xi_i_fks_fix*xiimax
      end if
      isolsign = 1
    elseif ((icountevts .eq. real_event) .and. &
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
    elseif (icountevts .eq. soft_counterevent) then
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
    do i = nincoming + 1, particle_count
      if (i .ne. i_fks .and. i .ne. j_fks .and. shybst .ne. 0.d0) &
        & call boostwdir2_in_place(chybst, shybst, chybstmo, xdir, &
        & xp(0, i))
    end do
!
! Phase-space factor for (xii,yij,phii)
    veckn = rho(xp(0, j_fks))
    veckbarn = rho(p_born_imother)
!
    xpswgt = xpswgt*2*shat/(4*pi)**3*veckn/veckbarn/ &
      & (2 - xi_i_fks*(1 - xp(0, j_fks)/veckn*y_ij_fks))
    xpswgt = abs(xpswgt)
    return
  end subroutine generate_momenta_massive_final

  subroutine generate_momenta_initial(icountevts, i_fks, j_fks, &
    & xbjrk_born, tau_born, ycm_born, ycmhat, shat_born, phi_i_fks, xp, x &
    & , shat, stot, sqrtshat, cm_boost, xbjrk, p_i_fks, xiimax, xinorm &
    & , xi_i_fks, y_ij_fks, xi_i_hat, xpswgt, xjac, particle_count, pass)
    implicit none
! arguments
    integer icountevts, i_fks, j_fks, particle_count
    double precision xbjrk_born(2), tau_born, ycm_born, ycmhat, shat_born &
      & , phi_i_fks, xpswgt, xjac, xiimax, xinorm, xp(0:3, nexternal), stot &
      & , x(2), y_ij_fks, xi_i_hat
    double precision shat, sqrtshat, cm_boost, xbjrk(2), p_i_fks(0:3)
    logical pass
! local
    integer i, j, idir
    double precision yijdir, costh_i_fks, x1bar2, x2bar2, yij_sol, xi1, xi2 &
      & , ximaxtmp, omega, bstfact, shy_tbst, chy_tbst, chy_tbstmo &
      & , xdir_t(3), cosphi_i_fks, sinphi_i_fks, shy_lbst, chy_lbst &
      & , encmso2, E_i_fks, sinth_i_fks, xpifksred(0:3), xi_i_fks &
      & , xiimin, yij_upp, yij_low, y_ij_fks_upp, y_ij_fks_low, tau
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
    if ((icountevts .eq. real_event .or. &
         icountevts .eq. soft_counterevent) .and. &
      & ((.not. softtest) .or. &
      & (softtest .and. y_ij_fks_fix .eq. -2.d0)) .and. &
      & (.not. colltest)) then
! importance sampling towards collinear singularity
! insert here further importance sampling towards y_ij_fks->1
      y_ij_fks = y_ij_fks_upp - &
        & (y_ij_fks_upp - y_ij_fks_low)*(cctiny + (1 - cctiny)*x(2)**2)
    elseif ((icountevts .eq. real_event .or. &
             icountevts .eq. soft_counterevent) .and. &
      & ((softtest .and. y_ij_fks_fix .ne. -2.d0) .or. &
      & colltest)) then
      y_ij_fks = y_ij_fks_fix
      if (y_ij_fks .gt. y_ij_fks_upp + 1d-12 .or. &
        & y_ij_fks .lt. y_ij_fks_low - 1d-12) then
        xjac = -33d0
        pass = .false.
        return
      end if
    elseif (icountevts .eq. soft_collinear_counterevent .or. &
            icountevts .eq. collinear_counterevent) then
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
    if (icountevts .ge. collinear_counterevent .and. &
        icountevts .le. last_counterevent .and. &
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
    if ((icountevts .eq. real_event .or. &
         icountevts .eq. collinear_counterevent) .and. &
      & ((.not. colltest) .or. &
      & (colltest .and. xi_i_fks_fix .eq. -2.d0)) .and. &
      & (.not. softtest)) then
      if (icountevts .eq. real_event) then
! importance sampling towards soft singularity
! insert here further importance sampling towards xi_i_hat->0
        xi_i_hat = sstiny + (1 - sstiny)*x(1)**2
      end if
      xi_i_fks = xiimin + (xiimax - xiimin)*xi_i_hat
    elseif ((icountevts .eq. real_event .or. &
             icountevts .eq. collinear_counterevent) .and. &
      & (colltest .and. xi_i_fks_fix .ne. -2.d0) .and. &
      & (.not. softtest)) then
      if (xi_i_fks_fix .lt. xiimax) then
        xi_i_fks = xi_i_fks_fix
      else
        xi_i_fks = xi_i_fks_fix*xiimax
      end if
    elseif ((icountevts .eq. real_event .or. &
             icountevts .eq. collinear_counterevent) .and. &
      & softtest) then
      if (xi_i_fks_fix .lt. xiimax) then
        xi_i_fks = xi_i_fks_fix
      else
        xjac = -102
        pass = .false.
        return
      end if
    elseif (icountevts .eq. soft_collinear_counterevent .or. &
            icountevts .eq. soft_counterevent) then
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
    if (icountevts .ne. soft_counterevent) then
      tau = tau_born/(1 - xi_i_fks)
      cm_boost = (ycm_born - log(omega)) - ycm_born
      shat = tau*stot
      sqrtshat = sqrt(shat)
      xbjrk(1) = xbjrk_born(1)/(sqrt(1 - xi_i_fks)*omega)
      xbjrk(2) = xbjrk_born(2)*omega/sqrt(1 - xi_i_fks)
    else
      tau = tau_born
      cm_boost = 0d0
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
    do i = 3, particle_count
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
    if (icountevts .eq. real_event .or. &
      & (icountevts .eq. collinear_counterevent .and. &
         xij_aor .eq. 0)) then
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

end module genps_fks
