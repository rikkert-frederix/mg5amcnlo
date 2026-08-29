module multiplicative_reweighter
  use process_dimensions, only: nexternal
  use multiplicative_density_contraction, only: multiplicative_density_basis
  use multiplicative_runtime, only: multiplicative_partonic_reweight, &
       evaluate_multiplicative_basis_reweight
  use multiplicative_scale_state, only: build_multiplicative_scale_tables
  implicit none
  private

  public :: multiplicative_partonic_reweight
  public :: reweight_multiplicative_scale_point

contains

  subroutine reweight_multiplicative_scale_point( &
       density_basis, event_slots, kinematic_weight, scale_factor_r, &
       scale_factor_f, decay_block_factor_indices, dynamic_scale, &
       virtual_sampling_fraction, virtual_sampled, reweight, &
       logarithmic_mu2_r, logarithmic_mu2_f, coupling_rescaling, &
       production_mu2_r, production_mu2_f)
    type(multiplicative_density_basis), intent(inout) :: density_basis
    integer, intent(in) :: event_slots(0:nexternal)
    double precision, intent(in) :: kinematic_weight
    double precision, intent(in) :: scale_factor_r, scale_factor_f
    integer, intent(in) :: decay_block_factor_indices(0:nexternal)
    integer, intent(in) :: dynamic_scale
    double precision, intent(in) :: virtual_sampling_fraction
    logical, intent(in) :: virtual_sampled
    type(multiplicative_partonic_reweight), intent(out) :: reweight
    double precision, intent(out) :: logarithmic_mu2_r(0:nexternal)
    double precision, intent(out) :: logarithmic_mu2_f(0:nexternal)
    double precision, intent(out) :: coupling_rescaling(0:nexternal, 0:1)
    double precision, intent(out) :: production_mu2_r, production_mu2_f

    call build_multiplicative_scale_tables( &
         event_slots, logarithmic_mu2_r, logarithmic_mu2_f, &
         coupling_rescaling, production_mu2_r, production_mu2_f, &
         scale_factor_r, scale_factor_f, decay_block_factor_indices, &
         dynamic_scale)
    call evaluate_multiplicative_basis_reweight( &
         density_basis, logarithmic_mu2_r, logarithmic_mu2_f, &
         coupling_rescaling, kinematic_weight, reweight, &
         virtual_sampling_fraction, virtual_sampled)
  end subroutine reweight_multiplicative_scale_point

end module multiplicative_reweighter
