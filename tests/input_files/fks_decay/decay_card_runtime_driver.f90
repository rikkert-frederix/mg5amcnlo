program decay_card_runtime_driver
  use decay_chain_parameters
  use decay_chain_scales
  use factorized_phase_space
  use nlo_decay_metadata, only: born_qcd
  use nlo_contribution_bundle, only: active, is_bundle
  use fnlo_scale_variations
  use run_state, only: do_rwgt_scale, do_rwgt_decay_scale
  use weight_lines, only: decay_scales, correction_scale_node
  use spin_density_weight_lines
  use spin_density_matrix_results, only: spin_density_bornlike_branch
  use setscales_module, only: set_alphas, set_ren_scale, set_fac_scale
  use production_scale_fixture, only: production_p, core_pdgs, core_count
  use run_state, only: mur2_current, muf12_current, muf22_current, qes2_current, &
       mur_over_ref, muf1_over_ref, muf2_over_ref, dynamical_scale_choice
  implicit none
  character(len=32) :: mode
  double precision :: p(0:3, 3), scales(8)
  double precision :: visible(0:3, 8), core(0:3, 8)
  double precision :: saved_p(0:3, 8), mur, muf(2)
  integer :: saved_pdgs(8)
  integer, allocatable :: factors(:)
  integer :: itop, iantitop, point, kr, kf, i
  character(len=80) :: label
  complex(kind=8) :: coefficients(3, 1, 1)
  call get_command_argument(1, mode)
  if (mode == 'no_decay') then
    is_bundle = .false.
    active = 0
    visible = 1d0
    call select_production_ren_scale_momenta(visible, 1, core)
    call decay_event_scales(0, scales)
    if (production_w_system_scale()) stop 3
    write(*, '(a,3(1x,l1))') 'NO_DECAY', all(core == visible), &
         all(scales == 0d0), .not. multiplicative_nlo_enabled()
    stop
  end if
  if (mode == 'qcd_born' .or. mode == 'qcd_split') born_qcd = 1
  if (mode == 'standalone' .or. mode == 'standalone_tbar') then
    is_bundle = .false.
    active = 2
    if (mode == 'standalone_tbar') active = 3
  end if
  call initialize_decay_chain_parameters()
  if (index(mode, 'production') == 1) then
    visible = 0d0
    visible(:, 3) = [200d0, 60d0, 0d0, 80d0]
    visible(:, 4) = [200d0, -60d0, 0d0, -80d0]
    visible(:, 5) = [50d0, 30d0, 0d0, 40d0]
    visible(:, 6) = [50d0, 0d0, 40d0, -30d0]
    select case (trim(mode))
    case ('production_real')
      visible(:, 7) = [10d0, 6d0, 8d0, 0d0]
    case ('production_soft')
      visible(:, 7) = 1d-8*[10d0, 6d0, 8d0, 0d0]
    case ('production_beam')
      visible(:, 7) = [10d0, 0d0, 0d0, 10d0]
    case ('production_virtuality')
      visible(:, 5:6) = 1.5d0*visible(:, 5:6)
    case ('production_explicit_w')
      visible(:, 5) = visible(:, 5) + visible(:, 6)
      visible(:, 6) = 0d0
      core_pdgs(5:8) = [24, 21, 0, 0]
      core_count = 6
    case ('production_swap')
      saved_p = visible
      saved_pdgs = core_pdgs
      visible(:, 3:7) = saved_p(:, [7, 6, 3, 5, 4])
      core_pdgs(3:7) = saved_pdgs([7, 6, 3, 5, 4])
    case ('production_minus')
      core_pdgs(5:6) = [11, -12]
    case ('production_bad_pair')
      core_pdgs(6) = 14
    case ('production_ambiguous')
      core_pdgs(7) = 14
    case ('production_bad_core')
      core_pdgs(4) = 6
    case ('production_bad_choice')
      dynamical_scale_choice = 2
    case ('production_decay_t')
      active = 2
    case ('production_decay_tbar')
      active = 3
    end select
    production_p = visible
    saved_p = visible
    call set_alphas(visible)
    write(*, '(a,4es25.16)') 'PRODUCTION_BASE ', &
         sqrt(mur2_current), sqrt(muf12_current), sqrt(muf22_current), sqrt(qes2_current)
    ! These are the same entry points used for production scale reweighting.
    mur_over_ref = 2d0
    muf1_over_ref = .5d0
    muf2_over_ref = .5d0
    call set_ren_scale(visible, mur)
    call set_fac_scale(visible, muf)
    write(*, '(a,3es25.16)') 'PRODUCTION_VARIED ', mur, muf
    write(*, '(a,1x,l1)') 'INPUT_UNCHANGED', all(visible == saved_p)
    stop
  end if
  if (index(mode, 'grid') == 1) then
    if (mode == 'grid_lo') call set_decay_run_order(.true.)
    do_rwgt_scale = mode /= 'grid_decay_only'
    call configure_fnlo_scale_variations()
    i = 0
    if (do_rwgt_decay_scale) i = decay_scale_species_count()
    allocate(factors(i))
    write(*, '(a,20i5)') 'AXES ', (decay_scale_species(i), i=1,size(factors))
    write(*, '(a,i5)') 'GRID_COUNT ', fnlo_scale_point_count(1)
    do point = 1, fnlo_scale_point_count(1)
      call decode_fnlo_scale_point(1, point, kr, kf, factors)
      call fnlo_scale_point_label(1, point, label)
      write(*, '(a,i0,1x,a)') 'POINT', point, trim(label)
      write(*, '(a,i0,20i5)') 'INDICES', point, kr, kf, factors
    end do
    stop
  end if
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
  if (mode == 'direct_split') p = 2d0*p
  call store_factorized_block_momenta(0, 1, 3, p)
  p(:, 2) = [50d0, 0d0, 40d0, 30d0]
  p(:, 3) = [50d0, 0d0, -40d0, -30d0]
  p(:, 1) = [100d0, 0d0, 0d0, 0d0]
  if (mode == 'direct_split') p = .5d0*p
  call store_factorized_block_momenta(0, 2, 3, p)
  call decay_event_scales(0, scales)
  allocate(factors(decay_scale_species_count()))
  factors = decay_scale_factor_count()
  if (.not. decay_scale_variation_enabled()) factors = 1
  itop = decay_scale_species_index(6)
  iantitop = decay_scale_species_index(-6)
  if (mode == 'split' .or. mode == 'qcd_split' .or. mode == 'standalone_tbar') then
    if (itop == iantitop) stop 2
    factors(itop) = 3
    factors(iantitop) = 2
  end if
  write(*, '(a,2es25.16)') 'SCALES ', scales(1:2)
  write(*, '(a,2es25.16)') 'WIDTHS ', &
       decay_nlo_width(6, 1, scales(1)), &
       decay_nlo_width(-6, factors(iantitop), scales(2))
  write(*, '(a,es25.16)') 'COUNTERTERM ', &
       decay_width_expansion_coefficient(factors, scales)
  write(*, '(a,es25.16)') 'DENOMINATOR ', &
       decay_width_denominator_rescaling(factors, scales)
  write(*, '(a,2es25.16)') 'LOCAL_WIDTHS ', &
       decay_node_width_rescaling(6, factors(itop), scales(1)), &
       decay_node_width_rescaling(-6, factors(iantitop), scales(2))
  write(*, '(a,es25.16)') 'PRODUCT_WIDTH ', &
       decay_multiplicative_width_rescaling(factors, scales)
  allocate(decay_scales(8, 2), correction_scale_node(2))
  decay_scales(:, 1) = scales
  decay_scales(:, 2) = scales
  correction_scale_node = [1, 2]
  coefficients = (1d0, 0d0)
  call record_spin_density_weight_line(1, 1, spin_density_bornlike_branch, &
       coefficients, 2, .false., 6)
  call record_spin_density_weight_line(2, 2, spin_density_bornlike_branch, &
       coefficients, 2, .false., -6)
  write(*, '(a,2es25.16)') 'DENSITY_MULTIPLIERS ', &
       spin_density_weight_line_multiplier(1, 1d0, factors, 1d0), &
       spin_density_weight_line_multiplier(2, 1d0, factors, 1d0)
  active = 2
  write(*, '(a,es25.16)') 'COUPLING ', &
       decay_qcd_coupling_weight(2 + 4*born_qcd, factors, scales)
  active = 3
  write(*, '(a,es25.16)') 'ANTITOP_COUPLING ', &
       decay_qcd_coupling_weight(2 + 4*born_qcd, factors, scales)
end program
