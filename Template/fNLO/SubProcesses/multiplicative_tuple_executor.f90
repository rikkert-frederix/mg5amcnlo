module multiplicative_tuple_executor
  use process_dimensions, only: nexternal
  use multiplicative_density_terms, only: block_nlo_distribution, &
       multiplicative_density_tuple
  use multiplicative_density_contraction, only: multiplicative_density_basis
  use multiplicative_kinematics, only: realize_factorized_event_transition
  use multiplicative_scale_state, only: build_multiplicative_scale_tables
  use multiplicative_runtime, only: multiplicative_event_evaluation, &
       prepare_multiplicative_event_kinematics, &
       evaluate_multiplicative_event_density
  implicit none
  private

  public :: execute_multiplicative_tuple
  public :: prepare_multiplicative_tuple
  public :: evaluate_prepared_multiplicative_tuple_density

contains

  subroutine execute_multiplicative_tuple( &
       distributions, tuple, previous_event_slots, previous_available, &
       vegas_weight, precision_asked, density_basis, &
       virtual_sampling_fraction, virtual_sampled, evaluation, &
       logarithmic_mu2_r, logarithmic_mu2_f, coupling_rescaling, &
       production_mu2_r, production_mu2_f, kinematics_available)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    type(multiplicative_density_tuple), intent(in) :: tuple
    integer, intent(inout) :: previous_event_slots(0:nexternal)
    logical, intent(inout) :: previous_available
    double precision, intent(in) :: vegas_weight, precision_asked
    type(multiplicative_density_basis), intent(inout) :: density_basis
    double precision, intent(in) :: virtual_sampling_fraction
    logical, intent(in) :: virtual_sampled
    type(multiplicative_event_evaluation), intent(inout) :: evaluation
    double precision, intent(out) :: logarithmic_mu2_r(0:nexternal)
    double precision, intent(out) :: logarithmic_mu2_f(0:nexternal)
    double precision, intent(out) :: coupling_rescaling(0:nexternal, 0:1)
    double precision, intent(out) :: production_mu2_r, production_mu2_f
    logical, intent(out) :: kinematics_available

    call prepare_multiplicative_tuple( &
         tuple, previous_event_slots, previous_available, vegas_weight, &
         evaluation, logarithmic_mu2_r, logarithmic_mu2_f, &
         coupling_rescaling, production_mu2_r, production_mu2_f, &
         kinematics_available)
    if (.not. kinematics_available .or. .not. evaluation%available) return
    call evaluate_prepared_multiplicative_tuple_density( &
         distributions, tuple, precision_asked, density_basis, &
         virtual_sampling_fraction, virtual_sampled, evaluation, &
         logarithmic_mu2_r, logarithmic_mu2_f, coupling_rescaling)
  end subroutine execute_multiplicative_tuple


  subroutine prepare_multiplicative_tuple( &
       tuple, previous_event_slots, previous_available, vegas_weight, &
       evaluation, logarithmic_mu2_r, logarithmic_mu2_f, &
       coupling_rescaling, production_mu2_r, production_mu2_f, &
       kinematics_available)
    type(multiplicative_density_tuple), intent(in) :: tuple
    integer, intent(inout) :: previous_event_slots(0:nexternal)
    logical, intent(inout) :: previous_available
    double precision, intent(in) :: vegas_weight
    type(multiplicative_event_evaluation), intent(inout) :: evaluation
    double precision, intent(out) :: logarithmic_mu2_r(0:nexternal)
    double precision, intent(out) :: logarithmic_mu2_f(0:nexternal)
    double precision, intent(out) :: coupling_rescaling(0:nexternal, 0:1)
    double precision, intent(out) :: production_mu2_r, production_mu2_f
    logical, intent(out) :: kinematics_available

    evaluation%available = .false.
    call realize_factorized_event_transition( &
         tuple%event_slots, previous_event_slots, previous_available, &
         kinematics_available)
    if (.not. kinematics_available) then
      previous_available = .false.
      return
    end if
    previous_event_slots = tuple%event_slots
    previous_available = .true.
    call build_multiplicative_scale_tables( &
         tuple%event_slots, logarithmic_mu2_r, logarithmic_mu2_f, &
         coupling_rescaling, production_mu2_r, production_mu2_f)
    call prepare_multiplicative_event_kinematics( &
         tuple, vegas_weight, evaluation)
  end subroutine prepare_multiplicative_tuple


  subroutine evaluate_prepared_multiplicative_tuple_density( &
       distributions, tuple, precision_asked, density_basis, &
       virtual_sampling_fraction, virtual_sampled, evaluation, &
       logarithmic_mu2_r, logarithmic_mu2_f, coupling_rescaling)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    type(multiplicative_density_tuple), intent(in) :: tuple
    double precision, intent(in) :: precision_asked
    type(multiplicative_density_basis), intent(inout) :: density_basis
    double precision, intent(in) :: virtual_sampling_fraction
    logical, intent(in) :: virtual_sampled
    type(multiplicative_event_evaluation), intent(inout) :: evaluation
    double precision, intent(in) :: logarithmic_mu2_r(0:nexternal)
    double precision, intent(in) :: logarithmic_mu2_f(0:nexternal)
    double precision, intent(in) :: coupling_rescaling(0:nexternal, 0:1)

    call evaluate_multiplicative_event_density( &
         distributions, tuple, logarithmic_mu2_r, logarithmic_mu2_f, &
         coupling_rescaling, precision_asked, evaluation, density_basis, &
         virtual_sampling_fraction, virtual_sampled)
  end subroutine evaluate_prepared_multiplicative_tuple_density

end module multiplicative_tuple_executor
