module multiplicative_runtime
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use process_dimensions, only: nexternal
  use factorized_phase_space, only: factorized_radiation_state, &
       fetch_factorized_radiation_state, compose_factorized_tuple_measure
  use nlo_contribution_bundle, only: multiplicative_event_capacity
  use fnlo_process_common, only: soft_counterevent
  use multiplicative_density_terms, only: block_nlo_distribution, &
       multiplicative_density_tuple, decode_density_cartesian_tuple
  use multiplicative_density_contraction, only: &
       multiplicative_density_basis, &
       contract_multiplicative_density_selection, &
       prepare_multiplicative_density_basis, &
       evaluate_multiplicative_density_basis, &
       evaluate_multiplicative_scale_polynomial
  use multiplicative_kinematics, only: &
       materialize_factorized_event_tuple
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

  public :: evaluate_multiplicative_event_tuple
  public :: evaluate_multiplicative_event_selection
  public :: evaluate_multiplicative_partonic_reweight
  public :: evaluate_multiplicative_basis_reweight

contains

  subroutine evaluate_multiplicative_event_tuple( &
       distributions, tuple_index, logarithmic_mu2_r, &
       logarithmic_mu2_f, coupling_rescaling, vegas_weight, &
       precision_asked, evaluation, already_realized)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    integer, intent(in) :: tuple_index
    double precision, intent(in) :: logarithmic_mu2_r(0:)
    double precision, intent(in) :: logarithmic_mu2_f(0:)
    double precision, intent(in) :: coupling_rescaling(0:, 0:)
    double precision, intent(in) :: vegas_weight, precision_asked
    type(multiplicative_event_evaluation), intent(inout) :: evaluation
    logical, intent(in), optional :: already_realized
    type(multiplicative_density_tuple) :: tuple

    call decode_density_cartesian_tuple(distributions, tuple_index, tuple)
    call evaluate_multiplicative_event_selection( &
         distributions, tuple, logarithmic_mu2_r, logarithmic_mu2_f, &
         coupling_rescaling, vegas_weight, precision_asked, evaluation, &
         already_realized)
  end subroutine evaluate_multiplicative_event_tuple


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
    type(multiplicative_density_basis), intent(inout), optional :: &
         density_basis
    double precision, intent(in), optional :: virtual_sampling_fraction
    logical, intent(in), optional :: virtual_sampled
    type(factorized_radiation_state) :: production_radiation
    complex(kind=8) :: density_result
    double precision :: jacobian, phase_space_weight
    integer :: capacity
    logical :: pass, measure_available, radiation_available
    logical :: use_realized_kinematics
    logical :: load_virtual_primitives, sampled_virtual
    double precision :: sampling_fraction

    if (vegas_weight < 0d0) then
      call fail_multiplicative_runtime('the VEGAS weight is negative')
    end if
    capacity = multiplicative_event_capacity()
    call prepare_event_evaluation(evaluation, capacity)
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
    load_virtual_primitives = &
         sampling_fraction >= 1d0 .or. sampled_virtual

    ! CONTRACT_MULTIPLICATIVE_DENSITY_TUPLE recursively boosts every selected
    ! block before calling any provider.  It does not assemble the visible
    ! event.  R/S/C/SC choices therefore remain distinct signed atoms while
    ! every matrix element sees only its boosted block-local momenta.
    if (present(density_basis)) then
      call prepare_multiplicative_density_basis( &
           distributions, tuple, precision_asked, density_basis, &
           use_realized_kinematics, load_virtual_primitives)
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
    else
      if (sampling_fraction < 1d0) then
        call fail_multiplicative_runtime( &
             'virtual sampling requires a reusable density basis')
      end if
      call contract_multiplicative_density_selection( &
           distributions, tuple, logarithmic_mu2_r, &
           logarithmic_mu2_f, precision_asked, density_result, &
           evaluation%nlo_order, evaluation%event_slots, &
           evaluation%precision_found, evaluation%return_code, &
           coupling_rescaling, use_realized_kinematics)
    end if
    ! MadLoop return codes are diagnostic bit fields, not a Boolean success
    ! flag.  Ordinary stable evaluations therefore generally have a nonzero
    ! code.  Reject only an exceptional result, a reported precision worse
    ! than the same five-percent ceiling used by the additive path, or a
    ! non-finite density.  Treating every nonzero code as a failure silently
    ! removed all virtual density matrices from multiplicative runs.
    if (.not. usable_multiplicative_density_result( &
        density_result, evaluation%precision_found, &
        evaluation%return_code)) return

    ! The global event is first made after the complete density contraction.
    ! This is the object on which cuts and observables must subsequently act.
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
    evaluation%partonic_weight = density_result*evaluation%kinematic_weight
    evaluation%available = .true.
  end subroutine evaluate_multiplicative_event_selection


  subroutine evaluate_multiplicative_partonic_reweight( &
       distributions, tuple, logarithmic_mu2_r, logarithmic_mu2_f, &
       coupling_rescaling, kinematic_weight, precision_asked, reweight)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    type(multiplicative_density_tuple), intent(in) :: tuple
    double precision, intent(in) :: logarithmic_mu2_r(0:)
    double precision, intent(in) :: logarithmic_mu2_f(0:)
    double precision, intent(in) :: coupling_rescaling(0:, 0:)
    double precision, intent(in) :: kinematic_weight, precision_asked
    type(multiplicative_partonic_reweight), intent(out) :: reweight
    integer :: event_slots(0:nexternal)
    complex(kind=8) :: density_result

    if (kinematic_weight <= 0d0 .or. &
        .not. ieee_is_finite(kinematic_weight)) return
    call contract_multiplicative_density_selection( &
         distributions, tuple, logarithmic_mu2_r, logarithmic_mu2_f, &
         precision_asked, density_result, reweight%nlo_order, event_slots, &
         reweight%precision_found, reweight%return_code, &
         coupling_rescaling, .true.)
    if (.not. usable_multiplicative_density_result( &
        density_result, reweight%precision_found, &
        reweight%return_code)) return
    reweight%partonic_weight = density_result*kinematic_weight
    reweight%available = .true.
  end subroutine evaluate_multiplicative_partonic_reweight


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


  subroutine evaluate_sampled_basis_density( &
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
      call evaluate_multiplicative_density_basis( &
           density_basis, logarithmic_mu2_r, logarithmic_mu2_f, result, &
           coupling_rescaling)
      return
    end if
    call evaluate_multiplicative_density_basis( &
         density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
         no_virtual_result, coupling_rescaling, .false.)
    if (.not. sampled_virtual) then
      result = no_virtual_result
      return
    end if
    call evaluate_multiplicative_density_basis( &
         density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
         full_result, coupling_rescaling, .true.)
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


  subroutine fail_multiplicative_runtime(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') 'ERROR in multiplicative_runtime: '//trim(message)
    stop 1
  end subroutine fail_multiplicative_runtime

end module multiplicative_runtime
