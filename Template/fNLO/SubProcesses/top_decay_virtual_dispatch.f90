module top_decay_virtual_dispatch
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use top_decay_virtual_cdr, only: tdv_virtual_rho_top, &
       tdv_virtual_rho_top_w
  use nlo_contribution_bundle, only: nlo_contribution_count
  implicit none
  private

  integer, parameter, public :: tdv_required_validation_points = 3
  double precision, parameter :: validation_tolerance = 1d-8
  logical, save :: validation_initialized = .false.
  integer, allocatable, save :: validation_count(:)
  double precision, allocatable, save :: validation_momenta(:, :, :, :)

  public :: tdv_evaluate_two_body_top
  public :: tdv_evaluate_two_body_top_w
  public :: tdv_two_body_analytic_supported
  public :: tdv_madloop_validation_required
  public :: tdv_validate_against_madloop

contains

  subroutine tdv_evaluate_two_body_top(parent_pdg, pt, pb, pw, mt, mb, &
       mw, mur, alphas, gc11, rho)
    integer, intent(in) :: parent_pdg
    double precision, intent(in) :: pt(0:3), pb(0:3), pw(0:3)
    double precision, intent(in) :: mt, mb, mw, mur, alphas
    complex(kind=8), intent(in) :: gc11
    complex(kind=8), intent(out) :: rho(3, 2, 2)
    double precision :: parity_pt(0:3), parity_pb(0:3), parity_pw(0:3)
    complex(kind=8) :: top_rho(3, 2, 2)

    call check_parent(parent_pdg)
    if (parent_pdg == 6) then
      call tdv_virtual_rho_top(pt, pb, pw, mt, mb, mw, mur, alphas, &
           gc11, rho)
      return
    end if

    call parity_reflect(pt, parity_pt)
    call parity_reflect(pb, parity_pb)
    call parity_reflect(pw, parity_pw)
    call tdv_virtual_rho_top(parity_pt, parity_pb, parity_pw, mt, mb, &
         mw, mur, alphas, gc11, top_rho)
    ! Charge conjugation maps the antitop density to a parity-reflected
    ! top density.  Match the spinor phase convention used by the generated
    ! MadLoop insertion before returning it to the common contraction code.
    if (sum(pt(1:3)**2) <= &
        64d0*epsilon(1d0)*max(pt(0)**2, 1d0)) then
      rho(:, 1, 1) = top_rho(:, 2, 2)
      rho(:, 1, 2) = top_rho(:, 2, 1)
      rho(:, 2, 1) = top_rho(:, 1, 2)
      rho(:, 2, 2) = top_rho(:, 1, 1)
    else
      call transform_antitop_density(pt, top_rho, rho)
    end if
  end subroutine tdv_evaluate_two_body_top


  subroutine tdv_evaluate_two_body_top_w(parent_pdg, pt, pb, pw, mt, &
       mb, mw, mur, alphas, gc11, rho)
    integer, intent(in) :: parent_pdg
    double precision, intent(in) :: pt(0:3), pb(0:3), pw(0:3)
    double precision, intent(in) :: mt, mb, mw, mur, alphas
    complex(kind=8), intent(in) :: gc11
    complex(kind=8), intent(out) :: rho(3, 6, 6)
    double precision :: parity_pt(0:3), parity_pb(0:3), parity_pw(0:3)
    complex(kind=8) :: top_rho(3, 6, 6)

    call check_parent(parent_pdg)
    if (parent_pdg == 6) then
      call tdv_virtual_rho_top_w(pt, pb, pw, mt, mb, mw, mur, alphas, &
           gc11, rho)
      return
    end if
    if (.not. tdv_two_body_analytic_supported(parent_pdg, pt, 6)) then
      call fail_dispatch('an exactly-resting antitop has an unsupported ' // &
           'joint top-W helicity basis')
    end if

    call parity_reflect(pt, parity_pt)
    call parity_reflect(pb, parity_pb)
    call parity_reflect(pw, parity_pw)
    call tdv_virtual_rho_top_w(parity_pt, parity_pb, parity_pw, mt, mb, &
         mw, mur, alphas, gc11, top_rho)
    call transform_antitop_density(pt, top_rho, rho)
  end subroutine tdv_evaluate_two_body_top_w


  logical function tdv_two_body_analytic_supported(parent_pdg, pt, &
       open_size)
    integer, intent(in) :: parent_pdg, open_size
    double precision, intent(in) :: pt(0:3)
    double precision :: spatial_scale

    call check_parent(parent_pdg)
    if (open_size /= 2 .and. open_size /= 6) then
      tdv_two_body_analytic_supported = .false.
      return
    end if
    spatial_scale = sum(pt(1:3)**2)
    tdv_two_body_analytic_supported = .not. (parent_pdg == -6 .and. &
         open_size == 6 .and. spatial_scale <= &
         64d0*epsilon(1d0)*max(pt(0)**2, 1d0))
  end function tdv_two_body_analytic_supported


  logical function tdv_madloop_validation_required(contribution)
    integer, intent(in) :: contribution

    call ensure_validation_state()
    call check_contribution(contribution)
    tdv_madloop_validation_required = &
         validation_count(contribution) < tdv_required_validation_points
  end function tdv_madloop_validation_required


  subroutine tdv_validate_against_madloop(contribution, momenta, &
       analytic_density, madloop_density, madloop_precision, &
       madloop_return_code)
    integer, intent(in) :: contribution, madloop_return_code
    double precision, intent(in) :: momenta(0:3, 3), madloop_precision
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
    do point = 1, validation_count(contribution)
      if (same_momenta(momenta, &
          validation_momenta(:, :, point, contribution))) then
        distinct = .false.
        exit
      end if
    end do
    if (.not. distinct) return

    point = validation_count(contribution) + 1
    validation_momenta(:, :, point, contribution) = momenta
    validation_count(contribution) = point
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

    if (validation_initialized) return
    contribution_count = nlo_contribution_count()
    if (contribution_count < 1) then
      call fail_dispatch('the NLO contribution count is invalid')
    end if
    allocate(validation_count(contribution_count))
    allocate(validation_momenta(0:3, 3, tdv_required_validation_points, &
                                contribution_count))
    validation_count = 0
    validation_momenta = 0d0
    validation_initialized = .true.
  end subroutine ensure_validation_state


  logical function same_momenta(first, second)
    double precision, intent(in) :: first(0:3, 3), second(0:3, 3)
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

    if (abs(parent_pdg) /= 6) then
      call fail_dispatch('the analytic provider received a non-top parent')
    end if
  end subroutine check_parent


  subroutine check_contribution(contribution)
    integer, intent(in) :: contribution

    if (contribution < 1 .or. contribution > size(validation_count)) then
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
