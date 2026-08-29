module multiplicative_runtime
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use process_dimensions, only: nexternal
  use factorized_phase_space, only: factorized_radiation_state, &
       fetch_factorized_radiation_state, compose_factorized_tuple_measure
  use nlo_contribution_bundle, only: multiplicative_event_capacity
  use fnlo_process_common, only: soft_counterevent
  use multiplicative_density_terms, only: block_nlo_distribution
  use multiplicative_density_contraction, only: &
       contract_multiplicative_density_tuple
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
    complex(kind=8) :: partonic_weight = (0d0, 0d0)
    integer, allocatable :: event_slots(:)
    integer, allocatable :: pdgs(:)
    integer, allocatable :: origin_blocks(:)
    double precision, allocatable :: momenta(:, :)
  end type multiplicative_event_evaluation

  public :: evaluate_multiplicative_event_tuple

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
    type(multiplicative_event_evaluation), intent(out) :: evaluation
    logical, intent(in), optional :: already_realized
    type(factorized_radiation_state) :: production_radiation
    complex(kind=8) :: density_result
    double precision :: jacobian, phase_space_weight
    integer :: capacity
    logical :: pass, measure_available, radiation_available
    logical :: use_realized_kinematics

    if (vegas_weight < 0d0) then
      call fail_multiplicative_runtime('the VEGAS weight is negative')
    end if
    capacity = multiplicative_event_capacity()
    allocate(evaluation%event_slots(0:nexternal))
    allocate(evaluation%momenta(0:3, capacity))
    allocate(evaluation%pdgs(capacity))
    allocate(evaluation%origin_blocks(capacity))
    use_realized_kinematics = .false.
    if (present(already_realized)) then
      use_realized_kinematics = already_realized
    end if

    ! CONTRACT_MULTIPLICATIVE_DENSITY_TUPLE recursively boosts every selected
    ! block before calling any provider.  It does not assemble the visible
    ! event.  R/S/C/SC choices therefore remain distinct signed atoms while
    ! every matrix element sees only its boosted block-local momenta.
    call contract_multiplicative_density_tuple( &
         distributions, tuple_index, logarithmic_mu2_r, &
         logarithmic_mu2_f, precision_asked, density_result, &
         evaluation%nlo_order, evaluation%event_slots, &
         evaluation%precision_found, evaluation%return_code, &
         coupling_rescaling, use_realized_kinematics)
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
    evaluation%partonic_weight = density_result*jacobian* &
         phase_space_weight*vegas_weight
    evaluation%available = .true.
  end subroutine evaluate_multiplicative_event_tuple


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
