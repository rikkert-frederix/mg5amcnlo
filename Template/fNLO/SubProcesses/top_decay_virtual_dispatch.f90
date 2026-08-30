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
  integer, parameter :: top_charge_phases(top_spin_size) = (/ 1, -1 /)
  integer, parameter :: top_w_charge_phases(top_w_spin_size) = &
       (/ 1, -1, -1, 1, 1, -1 /)
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
    complex(kind=8) :: top_rho(3, top_spin_size, top_spin_size)

    call check_parent(parent_pdg)
    analytic_available = .true.
    if (parent_pdg == top_pdg) then
      call tdv_virtual_rho_top(pt, pb, pw, mt, mb, mw, mur, alphas, &
           gc11, rho)
      call require_finite_analytic_density(rho)
      return
    end if

    call tdv_virtual_rho_top(pt, pb, pw, mt, mb, mw, mur, alphas, &
         gc11, top_rho)
    call charge_conjugate_density(top_rho, top_charge_phases, rho)
    call require_finite_analytic_density(rho)
  end subroutine tdv_evaluate_two_body_top


  subroutine tdv_evaluate_two_body_top_w(parent_pdg, pt, pb, pw, mt, &
       mb, mw, mur, alphas, gc11, rho, analytic_available)
    integer, intent(in) :: parent_pdg
    double precision, intent(in) :: pt(0:3), pb(0:3), pw(0:3)
    double precision, intent(in) :: mt, mb, mw, mur, alphas
    complex(kind=8), intent(in) :: gc11
    complex(kind=8), intent(out) :: rho(3, top_w_spin_size, top_w_spin_size)
    logical, intent(out) :: analytic_available
    complex(kind=8) :: top_rho(3, top_w_spin_size, top_w_spin_size)

    call check_parent(parent_pdg)
    analytic_available = .true.
    if (parent_pdg == top_pdg) then
      call tdv_virtual_rho_top_w(pt, pb, pw, mt, mb, mw, mur, alphas, &
           gc11, rho)
      call require_finite_analytic_density(rho)
      return
    end if

    call tdv_virtual_rho_top_w(pt, pb, pw, mt, mb, mw, mur, alphas, &
         gc11, top_rho)
    call charge_conjugate_density(top_rho, top_w_charge_phases, rho)
    call require_finite_analytic_density(rho)
  end subroutine tdv_evaluate_two_body_top_w


  subroutine tdv_evaluate_three_body_top(parent_pdg, pt, pb, pl, pnu, &
       mt, mb, mw, ww, mur, alphas, gc11, rho, analytic_available)
    integer, intent(in) :: parent_pdg
    double precision, intent(in) :: pt(0:3), pb(0:3), pl(0:3), pnu(0:3)
    double precision, intent(in) :: mt, mb, mw, ww, mur, alphas
    complex(kind=8), intent(in) :: gc11
    complex(kind=8), intent(out) :: rho(3, top_spin_size, top_spin_size)
    logical, intent(out) :: analytic_available
    complex(kind=8) :: top_rho(3, top_spin_size, top_spin_size)

    call check_parent(parent_pdg)
    analytic_available = .true.

    if (parent_pdg == top_pdg) then
      call tdv_virtual_rho_top_3body(pt, pb, pl, pnu, mt, mb, mw, mur, &
           alphas, gc11, rho, ww)
    else
      ! Keep pl and pnu assigned by physical PDG role.  Charge conjugation
      ! reverses the fermion flow through complex conjugation and the HELAS
      ! top-spin phases below; it does not exchange their momenta.
      call tdv_virtual_rho_top_3body(pt, pb, pl, pnu, mt, mb, mw, mur, &
           alphas, gc11, top_rho, ww)
      call charge_conjugate_density(top_rho, top_charge_phases, rho)
    end if
    call require_finite_analytic_density(rho)
  end subroutine tdv_evaluate_three_body_top


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


  subroutine charge_conjugate_density(top_density, phases, &
       antitop_density)
    complex(kind=8), intent(in) :: top_density(:, :, :)
    integer, intent(in) :: phases(:)
    complex(kind=8), intent(out) :: antitop_density(:, :, :)
    integer :: first, second

    if (size(top_density, 1) /= 3 .or. &
        any(shape(top_density) /= shape(antitop_density)) .or. &
        size(top_density, 2) /= size(top_density, 3) .or. &
        size(phases) /= size(top_density, 2)) then
      call fail_dispatch('the antitop density has an unsupported shape')
    end if

    do first = 1, size(top_density, 2)
      do second = 1, size(top_density, 2)
        antitop_density(:, first, second) = &
             phases(first)*phases(second)* &
             conjg(top_density(:, first, second))
      end do
    end do
  end subroutine charge_conjugate_density


  subroutine require_finite_analytic_density(density)
    complex(kind=8), intent(in) :: density(:, :, :)

    if (.not. density_is_finite(density)) then
      call fail_dispatch('the analytic density is non-finite')
    end if
  end subroutine require_finite_analytic_density


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
