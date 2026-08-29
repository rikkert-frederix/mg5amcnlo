module multiplicative_runtime
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use process_dimensions, only: nexternal
  use factorized_phase_space, only: factorized_radiation_state, &
       fetch_factorized_radiation_state, compose_factorized_tuple_measure
  use nlo_contribution_bundle, only: multiplicative_event_capacity
  use multiplicative_density_terms, only: block_nlo_distribution, &
       multiplicative_density_tuple
  use multiplicative_density_contraction, only: &
       multiplicative_density_basis, &
       prepare_multiplicative_density_basis, &
       evaluate_multiplicative_density_basis, &
       evaluate_multiplicative_scale_polynomial
  use multiplicative_kinematics, only: &
       realize_factorized_event_tuple, materialize_factorized_event_tuple
  implicit none
  private

  type, public :: multiplicative_event_evaluation
    logical :: available = .false.
    integer :: visible_count = 0
    integer :: nlo_order = 0
    integer :: return_code = 0
    double precision :: precision_found = 0d0
    double precision :: bjorken_x(2) = -1d0
    double precision :: y_to_lab = 0d0
    double precision :: kinematic_weight = 0d0
    complex(kind=8) :: partonic_weight = (0d0, 0d0)
    integer, allocatable :: event_slots(:)
    integer, allocatable :: pdgs(:)
    integer, allocatable :: origin_blocks(:)
    double precision, allocatable :: momenta(:, :)
  end type multiplicative_event_evaluation

  type, public :: multiplicative_partonic_reweight
    logical :: available = .false.
    integer :: nlo_order = 0
    integer :: return_code = 0
    double precision :: precision_found = 0d0
    complex(kind=8) :: partonic_weight = (0d0, 0d0)
  end type multiplicative_partonic_reweight

  public :: evaluate_multiplicative_event_selection
  public :: prepare_multiplicative_event_kinematics
  public :: evaluate_multiplicative_event_density
  public :: evaluate_multiplicative_basis_reweight
  public :: evaluate_multiplicative_basis_projection
  public :: real_multiplicative_weight

contains

  subroutine prepare_multiplicative_event_kinematics( &
       tuple, vegas_weight, evaluation)
    type(multiplicative_density_tuple), intent(in) :: tuple
    double precision, intent(in) :: vegas_weight
    type(multiplicative_event_evaluation), intent(inout) :: evaluation
    type(factorized_radiation_state) :: production_radiation
    double precision :: jacobian, phase_space_weight
    integer :: capacity
    logical :: pass, measure_available, radiation_available

    if (vegas_weight < 0d0) then
      call fail_multiplicative_runtime('the VEGAS weight is negative')
    end if
    capacity = multiplicative_event_capacity()
    call prepare_event_evaluation(evaluation, capacity)
    evaluation%event_slots = tuple%event_slots
    evaluation%nlo_order = tuple%nlo_order

    ! Cuts and observable-event identity need only the already-boosted local
    ! momenta.  Materialize that visible event before loading any matrix
    ! element so rejected fiducial points do not enter a provider.
    call materialize_factorized_event_tuple( &
         evaluation%event_slots, capacity, evaluation%visible_count, &
         evaluation%momenta, evaluation%pdgs, pass, &
         evaluation%origin_blocks)
    if (.not. pass) return
    call compose_factorized_tuple_measure( &
         evaluation%event_slots, jacobian, phase_space_weight, &
         measure_available)
    if (.not. measure_available .or. jacobian <= 0d0 .or. &
        phase_space_weight <= 0d0) return
    call fetch_factorized_radiation_state( &
         evaluation%event_slots(0), 0, production_radiation, &
         radiation_available)
    if (.not. radiation_available .or. &
        any(production_radiation%bjorken_x <= 0d0)) then
      call fail_multiplicative_runtime( &
           'the selected production atom has no Bjorken fractions')
    end if
    evaluation%bjorken_x = production_radiation%bjorken_x
    evaluation%y_to_lab = production_radiation%y_to_lab
    evaluation%kinematic_weight = jacobian*phase_space_weight*vegas_weight
    evaluation%available = .true.
  end subroutine prepare_multiplicative_event_kinematics


  subroutine evaluate_multiplicative_event_density( &
       distributions, tuple, logarithmic_mu2_r, logarithmic_mu2_f, &
       coupling_rescaling, precision_asked, evaluation, density_basis, &
       virtual_sampling_fraction, virtual_sampled)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    type(multiplicative_density_tuple), intent(in) :: tuple
    double precision, intent(in) :: logarithmic_mu2_r(0:)
    double precision, intent(in) :: logarithmic_mu2_f(0:)
    double precision, intent(in) :: coupling_rescaling(0:, 0:)
    double precision, intent(in) :: precision_asked
    type(multiplicative_event_evaluation), intent(inout) :: evaluation
    type(multiplicative_density_basis), intent(inout) :: density_basis
    double precision, intent(in), optional :: virtual_sampling_fraction
    logical, intent(in), optional :: virtual_sampled
    complex(kind=8) :: density_result
    logical :: load_virtual_primitives, sampled_virtual
    double precision :: sampling_fraction

    if (.not. evaluation%available) return
    evaluation%available = .false.
    evaluation%partonic_weight = (0d0, 0d0)
    sampling_fraction = 1d0
    if (present(virtual_sampling_fraction)) &
         sampling_fraction = virtual_sampling_fraction
    sampled_virtual = .true.
    if (present(virtual_sampled)) sampled_virtual = virtual_sampled
    call validate_virtual_sampling( &
         sampling_fraction, sampled_virtual, present(virtual_sampled))
    load_virtual_primitives = &
         sampling_fraction >= 1d0 .or. sampled_virtual

    call prepare_multiplicative_density_basis( &
         distributions, tuple, precision_asked, density_basis, .true., &
         load_virtual_primitives)
    evaluation%nlo_order = density_basis%nlo_order
    evaluation%event_slots = density_basis%event_slots
    evaluation%precision_found = density_basis%precision_found
    evaluation%return_code = density_basis%return_code
    density_result = (0d0, 0d0)
    if (density_basis%prepared) then
      call evaluate_sampled_basis_density( &
           density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
           coupling_rescaling, sampling_fraction, sampled_virtual, &
           density_result)
    end if
    if (.not. usable_multiplicative_density_result( &
        density_result, evaluation%precision_found, &
        evaluation%return_code)) return
    evaluation%partonic_weight = density_result*evaluation%kinematic_weight
    evaluation%available = .true.
  end subroutine evaluate_multiplicative_event_density

  subroutine evaluate_multiplicative_event_selection( &
       distributions, tuple, logarithmic_mu2_r, logarithmic_mu2_f, &
       coupling_rescaling, vegas_weight, precision_asked, evaluation, &
       already_realized, density_basis, virtual_sampling_fraction, &
       virtual_sampled)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    type(multiplicative_density_tuple), intent(in) :: tuple
    double precision, intent(in) :: logarithmic_mu2_r(0:)
    double precision, intent(in) :: logarithmic_mu2_f(0:)
    double precision, intent(in) :: coupling_rescaling(0:, 0:)
    double precision, intent(in) :: vegas_weight, precision_asked
    type(multiplicative_event_evaluation), intent(inout) :: evaluation
    logical, intent(in), optional :: already_realized
    type(multiplicative_density_basis), intent(inout) :: density_basis
    double precision, intent(in), optional :: virtual_sampling_fraction
    logical, intent(in), optional :: virtual_sampled
    logical :: use_realized_kinematics
    logical :: sampled_virtual
    logical :: pass
    double precision :: sampling_fraction

    use_realized_kinematics = .false.
    if (present(already_realized)) then
      use_realized_kinematics = already_realized
    end if
    sampling_fraction = 1d0
    if (present(virtual_sampling_fraction)) &
         sampling_fraction = virtual_sampling_fraction
    sampled_virtual = .true.
    if (present(virtual_sampled)) sampled_virtual = virtual_sampled
    call validate_virtual_sampling( &
         sampling_fraction, sampled_virtual, present(virtual_sampled))
    ! The public one-shot API retains its old contract.  The integration
    ! driver uses the two routines separately and applies cuts between them.
    if (.not. use_realized_kinematics) then
      call realize_factorized_event_tuple(tuple%event_slots, pass)
      if (.not. pass) return
    end if
    call prepare_multiplicative_event_kinematics( &
         tuple, vegas_weight, evaluation)
    if (.not. evaluation%available) return
    call evaluate_multiplicative_event_density( &
         distributions, tuple, logarithmic_mu2_r, logarithmic_mu2_f, &
         coupling_rescaling, precision_asked, evaluation, density_basis, &
         sampling_fraction, sampled_virtual)
  end subroutine evaluate_multiplicative_event_selection


  subroutine evaluate_multiplicative_basis_reweight( &
       density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
       coupling_rescaling, kinematic_weight, reweight, &
       virtual_sampling_fraction, virtual_sampled)
    type(multiplicative_density_basis), intent(inout) :: density_basis
    double precision, intent(in) :: logarithmic_mu2_r(0:)
    double precision, intent(in) :: logarithmic_mu2_f(0:)
    double precision, intent(in) :: coupling_rescaling(0:, 0:)
    double precision, intent(in) :: kinematic_weight
    type(multiplicative_partonic_reweight), intent(out) :: reweight
    double precision, intent(in), optional :: virtual_sampling_fraction
    logical, intent(in), optional :: virtual_sampled
    complex(kind=8) :: density_result
    double precision :: sampling_fraction
    logical :: sampled_virtual

    if (kinematic_weight <= 0d0 .or. &
        .not. ieee_is_finite(kinematic_weight)) return
    if (.not. density_basis%prepared) return
    sampling_fraction = 1d0
    if (present(virtual_sampling_fraction)) &
         sampling_fraction = virtual_sampling_fraction
    sampled_virtual = .true.
    if (present(virtual_sampled)) sampled_virtual = virtual_sampled
    call validate_virtual_sampling( &
         sampling_fraction, sampled_virtual, present(virtual_sampled))
    call evaluate_sampled_basis_scale_polynomial( &
         density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
         coupling_rescaling, sampling_fraction, sampled_virtual, &
         density_result)
    reweight%nlo_order = density_basis%nlo_order
    reweight%precision_found = density_basis%precision_found
    reweight%return_code = density_basis%return_code
    if (.not. usable_multiplicative_density_result( &
        density_result, reweight%precision_found, &
        reweight%return_code)) return
    reweight%partonic_weight = density_result*kinematic_weight
    reweight%available = .true.
  end subroutine evaluate_multiplicative_basis_reweight


  subroutine evaluate_multiplicative_basis_projection( &
       density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
       coupling_rescaling, kinematic_weight, component_orders, reweight, &
       virtual_sampling_fraction, virtual_sampled, radiation_position, &
       radiation_group)
    type(multiplicative_density_basis), intent(inout) :: density_basis
    double precision, intent(in) :: logarithmic_mu2_r(0:)
    double precision, intent(in) :: logarithmic_mu2_f(0:)
    double precision, intent(in) :: coupling_rescaling(0:, 0:)
    double precision, intent(in) :: kinematic_weight
    integer, intent(in) :: component_orders(:)
    type(multiplicative_partonic_reweight), intent(out) :: reweight
    double precision, intent(in), optional :: virtual_sampling_fraction
    logical, intent(in), optional :: virtual_sampled
    integer, intent(in), optional :: radiation_position, radiation_group
    complex(kind=8) :: density_result
    double precision :: sampling_fraction
    logical :: sampled_virtual

    if (kinematic_weight <= 0d0 .or. &
        .not. ieee_is_finite(kinematic_weight)) return
    if (.not. density_basis%prepared) return
    if (size(component_orders) /= density_basis%block_count .or. &
        any(component_orders < 0) .or. any(component_orders > 1)) then
      call fail_multiplicative_runtime( &
           'a projected formal-order vector has the wrong shape')
    end if
    sampling_fraction = 1d0
    if (present(virtual_sampling_fraction)) &
         sampling_fraction = virtual_sampling_fraction
    sampled_virtual = .true.
    if (present(virtual_sampled)) sampled_virtual = virtual_sampled
    call validate_virtual_sampling( &
         sampling_fraction, sampled_virtual, present(virtual_sampled))
    call evaluate_sampled_basis_density( &
         density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
         coupling_rescaling, sampling_fraction, sampled_virtual, &
         density_result, component_orders, radiation_position, &
         radiation_group)
    reweight%nlo_order = sum(component_orders)
    reweight%precision_found = density_basis%precision_found
    reweight%return_code = density_basis%return_code
    if (.not. usable_multiplicative_density_result( &
        density_result, reweight%precision_found, &
        reweight%return_code)) return
    reweight%partonic_weight = density_result*kinematic_weight
    reweight%available = .true.
  end subroutine evaluate_multiplicative_basis_projection


  subroutine evaluate_sampled_basis_density( &
       density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
       coupling_rescaling, sampling_fraction, sampled_virtual, result, &
       component_orders, radiation_position, radiation_group)
    type(multiplicative_density_basis), intent(inout) :: density_basis
    double precision, intent(in) :: logarithmic_mu2_r(0:)
    double precision, intent(in) :: logarithmic_mu2_f(0:)
    double precision, intent(in) :: coupling_rescaling(0:, 0:)
    double precision, intent(in) :: sampling_fraction
    logical, intent(in) :: sampled_virtual
    complex(kind=8), intent(out) :: result
    integer, intent(in), optional :: component_orders(:)
    integer, intent(in), optional :: radiation_position, radiation_group
    complex(kind=8) :: full_result, no_virtual_result

    if (sampling_fraction >= 1d0) then
      call evaluate_multiplicative_density_basis( &
           density_basis, logarithmic_mu2_r, logarithmic_mu2_f, result, &
           coupling_rescaling, component_orders=component_orders, &
           radiation_position=radiation_position, &
           radiation_group=radiation_group)
      return
    end if
    call evaluate_multiplicative_density_basis( &
         density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
         no_virtual_result, coupling_rescaling, .false., &
         component_orders, radiation_position, radiation_group)
    if (.not. sampled_virtual) then
      result = no_virtual_result
      return
    end if
    call evaluate_multiplicative_density_basis( &
         density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
         full_result, coupling_rescaling, .true., component_orders, &
         radiation_position, radiation_group)
    result = no_virtual_result + &
         (full_result - no_virtual_result)/sampling_fraction
  end subroutine evaluate_sampled_basis_density


  subroutine evaluate_sampled_basis_scale_polynomial( &
       density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
       coupling_rescaling, sampling_fraction, sampled_virtual, result)
    type(multiplicative_density_basis), intent(inout) :: density_basis
    double precision, intent(in) :: logarithmic_mu2_r(0:)
    double precision, intent(in) :: logarithmic_mu2_f(0:)
    double precision, intent(in) :: coupling_rescaling(0:, 0:)
    double precision, intent(in) :: sampling_fraction
    logical, intent(in) :: sampled_virtual
    complex(kind=8), intent(out) :: result
    complex(kind=8) :: full_result, no_virtual_result

    if (sampling_fraction >= 1d0) then
      call evaluate_multiplicative_scale_polynomial( &
           density_basis, logarithmic_mu2_r, logarithmic_mu2_f, result, &
           coupling_rescaling)
      return
    end if
    call evaluate_multiplicative_scale_polynomial( &
         density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
         no_virtual_result, coupling_rescaling, .false.)
    if (.not. sampled_virtual) then
      result = no_virtual_result
      return
    end if
    call evaluate_multiplicative_scale_polynomial( &
         density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
         full_result, coupling_rescaling, .true.)
    result = no_virtual_result + &
         (full_result - no_virtual_result)/sampling_fraction
  end subroutine evaluate_sampled_basis_scale_polynomial


  subroutine validate_virtual_sampling(fraction, sampled, mask_present)
    double precision, intent(in) :: fraction
    logical, intent(in) :: sampled, mask_present

    if (.not. ieee_is_finite(fraction) .or. &
        fraction <= 0d0 .or. fraction > 1d0) then
      call fail_multiplicative_runtime( &
           'the virtual-residual sampling fraction is invalid')
    end if
    if (fraction < 1d0 .and. .not. mask_present) then
      call fail_multiplicative_runtime( &
           'stochastic virtual sampling has no common point mask')
    end if
    if (fraction >= 1d0 .and. .not. sampled) then
      call fail_multiplicative_runtime( &
           'a unit virtual sampling fraction cannot skip virtuals')
    end if
  end subroutine validate_virtual_sampling


  subroutine prepare_event_evaluation(evaluation, capacity)
    type(multiplicative_event_evaluation), intent(inout) :: evaluation
    integer, intent(in) :: capacity

    evaluation%available = .false.
    evaluation%visible_count = 0
    evaluation%nlo_order = 0
    evaluation%return_code = 0
    evaluation%precision_found = 0d0
    evaluation%bjorken_x = -1d0
    evaluation%y_to_lab = 0d0
    evaluation%kinematic_weight = 0d0
    evaluation%partonic_weight = (0d0, 0d0)
    if (allocated(evaluation%event_slots)) then
      if (lbound(evaluation%event_slots, 1) /= 0 .or. &
          ubound(evaluation%event_slots, 1) /= nexternal) &
           deallocate(evaluation%event_slots)
    end if
    if (.not. allocated(evaluation%event_slots)) &
         allocate(evaluation%event_slots(0:nexternal))
    if (allocated(evaluation%momenta)) then
      if (size(evaluation%momenta, 1) /= 4 .or. &
          size(evaluation%momenta, 2) /= capacity) &
           deallocate(evaluation%momenta)
    end if
    if (.not. allocated(evaluation%momenta)) &
         allocate(evaluation%momenta(0:3, capacity))
    if (allocated(evaluation%pdgs)) then
      if (size(evaluation%pdgs) /= capacity) deallocate(evaluation%pdgs)
    end if
    if (.not. allocated(evaluation%pdgs)) allocate(evaluation%pdgs(capacity))
    if (allocated(evaluation%origin_blocks)) then
      if (size(evaluation%origin_blocks) /= capacity) &
           deallocate(evaluation%origin_blocks)
    end if
    if (.not. allocated(evaluation%origin_blocks)) &
         allocate(evaluation%origin_blocks(capacity))
    evaluation%event_slots = 0
    evaluation%momenta = 0d0
    evaluation%pdgs = 0
    evaluation%origin_blocks = -1
  end subroutine prepare_event_evaluation


  logical function usable_multiplicative_density_result( &
       density, precision_found, return_code)
    complex(kind=8), intent(in) :: density
    double precision, intent(in) :: precision_found
    integer, intent(in) :: return_code

    usable_multiplicative_density_result = .false.
    if (return_code < 0) return
    if (return_code/100 == 4) return
    if (.not. ieee_is_finite(precision_found) .or. &
        precision_found > 5d-2) return
    if (.not. ieee_is_finite(dble(density)) .or. &
        .not. ieee_is_finite(aimag(density))) return
    usable_multiplicative_density_result = .true.
  end function usable_multiplicative_density_result


  double precision function real_multiplicative_weight(weight)
    complex(kind=8), intent(in) :: weight

    if (.not. ieee_is_finite(dble(weight)) .or. &
        .not. ieee_is_finite(aimag(weight))) then
      call fail_multiplicative_runtime( &
           'a multiplicative density contraction is not finite')
    end if
    if (abs(aimag(weight)) > 1d-8*max(1d0, abs(dble(weight)))) then
      call fail_multiplicative_runtime( &
           'a multiplicative density contraction is not real')
    end if
    real_multiplicative_weight = dble(weight)
  end function real_multiplicative_weight


  subroutine fail_multiplicative_runtime(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') 'ERROR in multiplicative_runtime: '//trim(message)
    stop 1
  end subroutine fail_multiplicative_runtime

end module multiplicative_runtime
