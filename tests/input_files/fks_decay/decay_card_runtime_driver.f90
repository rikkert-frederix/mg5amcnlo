program decay_card_runtime_driver
  use decay_chain_parameters
  use decay_chain_scales
  use factorized_phase_space
  use nlo_decay_metadata, only: born_qcd
  use nlo_contribution_bundle, only: active, is_bundle
  implicit none
  character(len=32) :: mode
  double precision :: p(0:3, 3), scales(8)
  integer :: factors(1)
  call get_command_argument(1, mode)
  if (mode == 'qcd_born') born_qcd = 1
  if (mode == 'standalone') then
    is_bundle = .false.
    active = 2
  end if
  call initialize_decay_chain_parameters()
  if (mode == 'orders' .or. mode == 'all_lo') then
    if (mode == 'all_lo') call set_decay_run_order(.true.)
    write(*, '(a,4(1x,l1))') 'ORDERS', nlo_correction_enabled(1), &
         nlo_correction_enabled(2), decay_species_nlo_enabled(6), &
         multiplicative_nlo_enabled()
    write(*, '(a,es25.16)') 'COUNTERTERM ', decay_width_expansion_coefficient()
    stop
  end if
  p(:, 1) = [100d0, 0d0, 0d0, 0d0]
  p(:, 2) = [50d0, 30d0, 0d0, 40d0]
  p(:, 3) = [50d0, -30d0, 0d0, -40d0]
  call store_factorized_block_momenta(0, 1, 3, p)
  p(:, 2) = [50d0, 0d0, 40d0, 30d0]
  p(:, 3) = [50d0, 0d0, -40d0, -30d0]
  call store_factorized_block_momenta(0, 2, 3, p)
  call decay_event_scales(0, scales)
  factors = decay_scale_factor_count()
  if (.not. decay_scale_variation_enabled()) factors = 1
  write(*, '(a,2es25.16)') 'SCALES ', scales(1:2)
  write(*, '(a,2es25.16)') 'WIDTHS ', &
       decay_nlo_width(6, 1, scales(1)), &
       decay_nlo_width(-6, factors(1), scales(2))
  write(*, '(a,es25.16)') 'COUNTERTERM ', &
       decay_width_expansion_coefficient(factors, scales)
  write(*, '(a,es25.16)') 'DENOMINATOR ', &
       decay_width_denominator_rescaling(factors, scales)
  write(*, '(a,2es25.16)') 'LOCAL_WIDTHS ', &
       decay_node_width_rescaling(6, factors(1), scales(1)), &
       decay_node_width_rescaling(-6, factors(1), scales(2))
  active = 2
  write(*, '(a,es25.16)') 'COUPLING ', &
       decay_qcd_coupling_weight(2 + 4*born_qcd, factors, scales)
end program
