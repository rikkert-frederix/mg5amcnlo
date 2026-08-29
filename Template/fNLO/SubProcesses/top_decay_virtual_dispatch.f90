module top_decay_virtual_dispatch
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use top_decay_virtual_cdr, only: tdv_virtual_rho_top, &
       tdv_virtual_rho_top_w, tdv_virtual_rho_top_3body
  use nlo_contribution_bundle, only: nlo_contribution_count
  implicit none
  private

  integer, parameter, public :: tdv_required_validation_points = 3
  integer, parameter :: top_pdg = 6
  integer, parameter :: top_spin_size = 2
  integer, parameter :: top_w_spin_size = 6
  integer, parameter :: decay_momentum_count = 4
  double precision, parameter :: validation_tolerance = 1d-8
  integer, allocatable, save :: validated_point_count(:)
  double precision, allocatable, save :: validated_momenta(:, :, :, :)

  public :: tdv_evaluate_two_body_top
  public :: tdv_evaluate_two_body_top_w
  public :: tdv_evaluate_three_body_top
  public :: tdv_madloop_required
  public :: tdv_validate_against_madloop

contains

  subroutine tdv_evaluate_two_body_top(parent_pdg, pt, pb, pw, mt, mb, &
       mw, mur, alphas, gc11, rho, analytic_available)
    integer, intent(in) :: parent_pdg
    double precision, intent(in) :: pt(0:3), pb(0:3), pw(0:3)
    double precision, intent(in) :: mt, mb, mw, mur, alphas
    complex(kind=8), intent(in) :: gc11
    complex(kind=8), intent(out) :: rho(3, top_spin_size, top_spin_size)
    logical, intent(out) :: analytic_available
    double precision :: parity_pt(0:3), parity_pb(0:3), parity_pw(0:3)
    complex(kind=8) :: top_rho(3, top_spin_size, top_spin_size)

    call check_parent(parent_pdg)
    analytic_available = .true.
    if (parent_pdg == top_pdg) then
      call tdv_virtual_rho_top(pt, pb, pw, mt, mb, mw, mur, alphas, &
           gc11, rho)
      return
    end if

    call parity_reflect(pt, parity_pt)
    call parity_reflect(pb, parity_pb)
    call parity_reflect(pw, parity_pw)
    call tdv_virtual_rho_top(parity_pt, parity_pb, parity_pw, mt, mb, &
         mw, mur, alphas, gc11, top_rho)
    call transform_top_only_antitop_density(pt, top_rho, rho)
  end subroutine tdv_evaluate_two_body_top


  subroutine tdv_evaluate_two_body_top_w(parent_pdg, pt, pb, pw, mt, &
       mb, mw, mur, alphas, gc11, rho, analytic_available)
    integer, intent(in) :: parent_pdg
    double precision, intent(in) :: pt(0:3), pb(0:3), pw(0:3)
    double precision, intent(in) :: mt, mb, mw, mur, alphas
    complex(kind=8), intent(in) :: gc11
    complex(kind=8), intent(out) :: rho(3, top_w_spin_size, top_w_spin_size)
    logical, intent(out) :: analytic_available
    double precision :: parity_pt(0:3), parity_pb(0:3), parity_pw(0:3)
    complex(kind=8) :: top_rho(3, top_w_spin_size, top_w_spin_size)

    call check_parent(parent_pdg)
    analytic_available = joint_density_supported(parent_pdg, pt)
    if (.not. analytic_available) then
      rho = (0d0, 0d0)
      return
    end if
    if (parent_pdg == top_pdg) then
      call tdv_virtual_rho_top_w(pt, pb, pw, mt, mb, mw, mur, alphas, &
           gc11, rho)
      return
    end if

    call parity_reflect(pt, parity_pt)
    call parity_reflect(pb, parity_pb)
    call parity_reflect(pw, parity_pw)
    call tdv_virtual_rho_top_w(parity_pt, parity_pb, parity_pw, mt, mb, &
         mw, mur, alphas, gc11, top_rho)
    call transform_antitop_density(pt, top_rho, rho)
  end subroutine tdv_evaluate_two_body_top_w


  subroutine tdv_evaluate_three_body_top(parent_pdg, pt, pb, pl, pnu, &
       mt, mb, mw, ww, mur, alphas, gc11, rho, analytic_available)
    integer, intent(in) :: parent_pdg
    double precision, intent(in) :: pt(0:3), pb(0:3), pl(0:3), pnu(0:3)
    double precision, intent(in) :: mt, mb, mw, ww, mur, alphas
    complex(kind=8), intent(in) :: gc11
    complex(kind=8), intent(out) :: rho(3, top_spin_size, top_spin_size)
    logical, intent(out) :: analytic_available
    double precision :: parity_pt(0:3), parity_pb(0:3)
    double precision :: parity_pl(0:3), parity_pnu(0:3)
    double precision :: q(0:3), q2, pole_distance, width_factor
    complex(kind=8) :: top_rho(3, top_spin_size, top_spin_size)

    call check_parent(parent_pdg)
    q = pl + pnu
    q2 = q(0)**2 - sum(q(1:3)**2)
    pole_distance = q2 - mw**2
    analytic_available = three_body_kinematics_supported( &
         pt, pb, pl, pnu, mt, mb, mw, q2, pole_distance)
    if (.not. analytic_available) then
      rho = (0d0, 0d0)
      return
    end if

    if (parent_pdg == top_pdg) then
      call tdv_virtual_rho_top_3body(pt, pb, pl, pnu, mt, mb, mw, mur, &
           alphas, gc11, rho)
    else
      call parity_reflect(pt, parity_pt)
      call parity_reflect(pb, parity_pb)
      call parity_reflect(pl, parity_pl)
      call parity_reflect(pnu, parity_pnu)
      call tdv_virtual_rho_top_3body(parity_pt, parity_pb, parity_pl, &
           parity_pnu, mt, mb, mw, mur, alphas, gc11, top_rho)
      call transform_top_only_antitop_density(pt, top_rho, rho)
    end if

    ! The packaged kernel has a zero-width W propagator.  For the fixed-width
    ! propagator used by the generated tree and loop amplitudes, only its
    ! common modulus squared changes because the massless leptonic current is
    ! transverse.  Points too close to the zero-width pole use MadLoop above.
    width_factor = pole_distance**2/(pole_distance**2 + (mw*ww)**2)
    rho = width_factor*rho
    if (.not. density_is_finite(rho)) then
      analytic_available = .false.
      rho = (0d0, 0d0)
    end if
  end subroutine tdv_evaluate_three_body_top


  logical function three_body_kinematics_supported( &
       pt, pb, pl, pnu, mt, mb, mw, q2, pole_distance)
    double precision, intent(in) :: pt(0:3), pb(0:3), pl(0:3), pnu(0:3)
    double precision, intent(in) :: mt, mb, mw, q2, pole_distance
    double precision :: scale

    scale = max(mt**2, mw**2, 1d0)
    three_body_kinematics_supported = &
         all(ieee_is_finite(pt)) .and. all(ieee_is_finite(pb)) .and. &
         all(ieee_is_finite(pl)) .and. all(ieee_is_finite(pnu)) .and. &
         ieee_is_finite(q2) .and. ieee_is_finite(pole_distance) .and. &
         q2 > 0d0 .and. q2 < (mt-mb)**2 .and. &
         abs(pole_distance) > 1d-8*scale
  end function three_body_kinematics_supported


  pure logical function joint_density_supported(parent_pdg, pt)
    integer, intent(in) :: parent_pdg
    double precision, intent(in) :: pt(0:3)

    joint_density_supported = .not. (parent_pdg == -top_pdg .and. &
         effectively_at_rest(pt))
  end function joint_density_supported


  logical function tdv_madloop_required(contribution, analytic_available)
    integer, intent(in) :: contribution
    logical, intent(in) :: analytic_available

    if (.not. analytic_available) then
      tdv_madloop_required = .true.
      return
    end if
    call ensure_validation_state()
    call check_contribution(contribution)
    tdv_madloop_required = validated_point_count(contribution) < &
         tdv_required_validation_points
  end function tdv_madloop_required


  subroutine tdv_validate_against_madloop(contribution, momenta, &
       analytic_density, madloop_density, madloop_precision, &
       madloop_return_code)
    integer, intent(in) :: contribution, madloop_return_code
    double precision, intent(in) :: momenta(0:3, decay_momentum_count)
    double precision, intent(in) :: madloop_precision
    complex(kind=8), intent(in) :: analytic_density(:, :, :)
    complex(kind=8), intent(in) :: madloop_density(:, :, :)
    double precision :: difference, relative_difference, scale
    integer :: location(3), point
    logical :: distinct

    call ensure_validation_state()
    call check_contribution(contribution)
    if (size(analytic_density, 1) /= 3 .or. &
        any(shape(analytic_density) /= shape(madloop_density)) .or. &
        size(analytic_density, 2) /= size(analytic_density, 3)) then
      call fail_dispatch('the validation densities have incompatible shapes')
    end if
    if (.not. all(ieee_is_finite(momenta))) then
      call fail_dispatch('the validation point has non-finite momenta')
    end if
    if (.not. density_is_finite(analytic_density)) then
      call fail_dispatch('the analytic validation density is non-finite')
    end if
    if (.not. density_is_finite(madloop_density)) then
      call fail_dispatch('the MadLoop validation density is non-finite')
    end if

    difference = maxval(abs(analytic_density - madloop_density))
    scale = max(maxval(abs(analytic_density)), &
                maxval(abs(madloop_density)), tiny(1d0))
    relative_difference = difference/scale
    if (relative_difference > validation_tolerance) then
      location = maxloc(abs(analytic_density - madloop_density))
      write (*,*) 'ERROR: analytic top-decay virtual validation failed'
      write (*,*) 'Contribution:', contribution
      write (*,*) 'Density index:', location
      write (*,*) 'Relative maximum difference:', relative_difference
      write (*,*) 'Analytic value:', &
           analytic_density(location(1), location(2), location(3))
      write (*,*) 'MadLoop value:', &
           madloop_density(location(1), location(2), location(3))
      write (*,*) 'MadLoop precision and return code:', &
           madloop_precision, madloop_return_code
      stop 1
    end if

    distinct = .true.
    do point = 1, validated_point_count(contribution)
      if (same_momenta(momenta, &
          validated_momenta(:, :, point, contribution))) then
        distinct = .false.
        exit
      end if
    end do
    if (.not. distinct) return

    point = validated_point_count(contribution) + 1
    validated_momenta(:, :, point, contribution) = momenta
    validated_point_count(contribution) = point
    write (*,'(A,I0,A,I0,A,I0,A,ES12.4)') &
         'INFO: analytic top-decay virtual validation ', point, '/', &
         tdv_required_validation_points, ' for contribution ', &
         contribution, '; relative difference ', relative_difference
    if (point == tdv_required_validation_points) then
      write (*,'(A,I0,A)') &
           'INFO: analytic top-decay virtual provider for contribution ', &
           contribution, &
           ' validated; subsequent phase-space points skip MadLoop.'
    end if
  end subroutine tdv_validate_against_madloop


  subroutine ensure_validation_state()
    integer :: contribution_count

    if (allocated(validated_point_count)) return
    contribution_count = nlo_contribution_count()
    if (contribution_count < 1) then
      call fail_dispatch('the NLO contribution count is invalid')
    end if
    allocate(validated_point_count(contribution_count))
    allocate(validated_momenta(0:3, decay_momentum_count, &
         tdv_required_validation_points, contribution_count))
    validated_point_count = 0
  end subroutine ensure_validation_state


  logical function same_momenta(first, second)
    double precision, intent(in) :: first(0:3, decay_momentum_count)
    double precision, intent(in) :: second(0:3, decay_momentum_count)
    double precision :: scale

    scale = max(maxval(abs(first)), maxval(abs(second)), 1d0)
    same_momenta = maxval(abs(first - second)) <= &
         128d0*epsilon(1d0)*scale
  end function same_momenta


  logical function density_is_finite(density)
    complex(kind=8), intent(in) :: density(:, :, :)

    density_is_finite = all(ieee_is_finite(real(density, kind=8))) .and. &
         all(ieee_is_finite(aimag(density)))
  end function density_is_finite


  pure logical function effectively_at_rest(momentum)
    double precision, intent(in) :: momentum(0:3)

    effectively_at_rest = sum(momentum(1:3)**2) <= &
         64d0*epsilon(1d0)*max(momentum(0)**2, 1d0)
  end function effectively_at_rest


  subroutine transform_top_only_antitop_density(pt, top_density, &
       antitop_density)
    double precision, intent(in) :: pt(0:3)
    complex(kind=8), intent(in) :: top_density(3, top_spin_size, &
         top_spin_size)
    complex(kind=8), intent(out) :: antitop_density(3, top_spin_size, &
         top_spin_size)

    ! Charge conjugation maps the antitop density to a parity-reflected top
    ! density.  Match the spinor phase convention used by the generated
    ! MadLoop insertion before returning it to the common contraction code.
    if (effectively_at_rest(pt)) then
      antitop_density(:, 1, 1) = top_density(:, 2, 2)
      antitop_density(:, 1, 2) = top_density(:, 2, 1)
      antitop_density(:, 2, 1) = top_density(:, 1, 2)
      antitop_density(:, 2, 2) = top_density(:, 1, 1)
    else
      call transform_antitop_density(pt, top_density, antitop_density)
    end if
  end subroutine transform_top_only_antitop_density


  subroutine transform_antitop_density(pt, top_density, antitop_density)
    double precision, intent(in) :: pt(0:3)
    complex(kind=8), intent(in) :: top_density(:, :, :)
    complex(kind=8), intent(out) :: antitop_density(:, :, :)
    complex(kind=8) :: phase(6)
    complex(kind=8) :: even_phase
    double precision :: transverse_squared
    integer :: first, second

    if (size(top_density, 1) /= 3 .or. &
        any(shape(top_density) /= shape(antitop_density)) .or. &
        size(top_density, 2) /= size(top_density, 3) .or. &
        (size(top_density, 2) /= 2 .and. &
         size(top_density, 2) /= 6)) then
      call fail_dispatch('the antitop density has an unsupported shape')
    end if

    transverse_squared = pt(1)**2 + pt(2)**2
    if (transverse_squared <= &
        64d0*epsilon(1d0)*max(sum(pt(1:3)**2), 1d0)) then
      even_phase = (-1d0, 0d0)
    else
      even_phase = cmplx( &
           -(pt(1)**2 - pt(2)**2)/transverse_squared, &
           2d0*pt(1)*pt(2)/transverse_squared, kind=8)
    end if
    phase(1:size(top_density, 2)) = (1d0, 0d0)
    phase(2:size(top_density, 2):2) = even_phase
    do first = 1, size(top_density, 2)
      do second = 1, size(top_density, 2)
        antitop_density(:, first, second) = &
             phase(first)*top_density(:, first, second)*conjg(phase(second))
      end do
    end do
  end subroutine transform_antitop_density


  pure subroutine parity_reflect(momentum, reflected)
    double precision, intent(in) :: momentum(0:3)
    double precision, intent(out) :: reflected(0:3)

    reflected(0) = momentum(0)
    reflected(1:3) = -momentum(1:3)
  end subroutine parity_reflect


  subroutine check_parent(parent_pdg)
    integer, intent(in) :: parent_pdg

    if (abs(parent_pdg) /= top_pdg) then
      call fail_dispatch('the analytic provider received a non-top parent')
    end if
  end subroutine check_parent


  subroutine check_contribution(contribution)
    integer, intent(in) :: contribution

    if (contribution < 1 .or. &
        contribution > size(validated_point_count)) then
      call fail_dispatch('the validation contribution is out of range')
    end if
  end subroutine check_contribution


  subroutine fail_dispatch(message)
    character(len=*), intent(in) :: message

    write (*,*) 'ERROR: analytic top-decay virtual dispatch: ', &
         trim(message)
    stop 1
  end subroutine fail_dispatch

end module top_decay_virtual_dispatch
