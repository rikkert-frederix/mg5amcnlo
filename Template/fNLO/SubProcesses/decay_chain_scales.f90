module decay_chain_scales
  use process_dimensions, only: nexternal
  use decay_chain_metadata, only: has_decay_chains, decay_node_count, &
       context_for_fks, node_pdg, node_qcd_order
  use decay_chain_parameters, only: decay_renormalization_scale, &
       use_decayed_production_ren_scale_momenta, &
       decay_scale_species_count, decay_scale_species_index, &
       decay_dynamical_scale_choice, decay_scale_factor
  use factorized_phase_space, only: factorized_block_size, &
       fetch_factorized_block_momenta
  use fixed_order_user_hooks, only: fixed_user_decay_scale
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use decay_chain_kinematics, only: contract_visible_momenta
  use nlo_decay_metadata, only: has_nlo_decay, corrected_parent_pdg, &
       nlo_decay_production_born_qcd_order, nlo_decay_born_qcd_order, &
       nlo_decay_node_count, nlo_decay_node_pdg, &
       nlo_decay_node_qcd_order, nlo_decay_corrected_node
  use nlo_decay_kinematics, only: get_nlo_decay_production_momenta
  use alfas_functions_module, only: alphas
  implicit none
  private

  public :: decay_qcd_squared_order, production_qcd_squared_order
  public :: active_block_qcd_squared_order
  public :: decay_qcd_coupling_weight
  public :: decay_qcd_coupling_rescaling
  public :: corrected_born_qcd_squared_order
  public :: select_production_core_momenta
  public :: select_production_ren_scale_momenta
  public :: decay_event_scales

contains

  subroutine decay_event_scales(event_slot, scales)
    integer, intent(in) :: event_slot
    double precision, intent(out) :: scales(:)
    integer :: node, nodes, pdg, choice, particles, i
    double precision :: p(0:3, nexternal), value, spatial2
    logical :: available

    scales = 0d0
    if (has_nlo_decay()) then
      nodes = nlo_decay_node_count()
    else if (has_decay_chains()) then
      nodes = decay_node_count()
    else
      return
    end if
    if (size(scales) < nodes) call fail_scales('too few decay event scales')
    do node = 1, nodes
      if (has_nlo_decay()) then
        pdg = nlo_decay_node_pdg(node)
      else
        pdg = node_pdg(node)
      end if
      value = decay_renormalization_scale(pdg)
      choice = decay_dynamical_scale_choice(pdg)
      if (choice /= 0) then
        particles = factorized_block_size(event_slot, node)
        if (particles < 2) call fail_scales('dynamic decay scale has no block momenta')
        call fetch_factorized_block_momenta( &
             event_slot, node, particles, p(:, 1:particles), available)
        if (.not. available) call fail_scales('dynamic decay scale block is absent')
        if (choice == -1) then
          value = fixed_user_decay_scale(pdg, node, p(:, 1:particles), value)
        else
          value = 0d0
          do i = 2, particles
            if (choice == 1) then
              spatial2 = sum(p(1:3, i)**2)
              if (spatial2 > 0d0) &
                   value = value + p(0, i)*sqrt(sum(p(1:2, i)**2)/spatial2)
            else
              value = value + sqrt(max(0d0, (p(0, i)-p(3, i))*(p(0, i)+p(3, i))))
            end if
          end do
          if (choice == 3) value = value/2d0
        end if
        if (.not. ieee_is_finite(value) .or. value < 0d0) &
             call fail_scales('invalid dynamical decay scale')
        value = max(2d0, value)
      end if
      scales(node) = value
    end do
  end subroutine decay_event_scales

  integer function active_block_qcd_squared_order(total_qcd_order)
    integer, intent(in) :: total_qcd_order

    if (has_nlo_decay()) then
      ! The direct density provider contains only the corrected decay block;
      ! remove production and every spectator decay from the flattened order.
      active_block_qcd_squared_order = total_qcd_order - &
           nlo_decay_production_born_qcd_order() - &
           nlo_decay_static_qcd_order()
    else if (has_decay_chains()) then
      ! In the bundle member without NLO-decay metadata the active block is
      ! production.  All decay orders are spectators of this local density.
      active_block_qcd_squared_order = &
           production_qcd_squared_order(total_qcd_order)
    else
      active_block_qcd_squared_order = total_qcd_order
    end if
    if (active_block_qcd_squared_order < 0) then
      call fail_scales('the active block has a negative QCD order')
    end if
  end function active_block_qcd_squared_order

  integer function decay_qcd_squared_order(total_qcd_order)
    integer, intent(in), optional :: total_qcd_order
    integer :: node

    decay_qcd_squared_order = 0
    if (has_nlo_decay()) then
      if (present(total_qcd_order)) then
        decay_qcd_squared_order = total_qcd_order - &
             nlo_decay_production_born_qcd_order()
      else
        do node = 1, nlo_decay_node_count()
          decay_qcd_squared_order = decay_qcd_squared_order + &
               2*nlo_decay_node_qcd_order(node)
        end do
      end if
      if (decay_qcd_squared_order < nlo_decay_born_total_qcd_order()) then
        call fail_scales('total QCD order is below the NLO-decay Born order')
      end if
      return
    end if
    if (.not. has_decay_chains()) return
    do node = 1, decay_node_count()
      decay_qcd_squared_order = decay_qcd_squared_order + &
           2*node_qcd_order(node)
    end do
  end function decay_qcd_squared_order


  integer function production_qcd_squared_order(total_qcd_order)
    integer, intent(in) :: total_qcd_order

    if (has_nlo_decay()) then
      production_qcd_squared_order = &
           nlo_decay_production_born_qcd_order()
      if (total_qcd_order - production_qcd_squared_order < &
          nlo_decay_born_total_qcd_order()) then
        call fail_scales('NLO-decay QCD orders are inconsistent')
      end if
      return
    end if
    production_qcd_squared_order = &
         total_qcd_order - decay_qcd_squared_order()
    if (production_qcd_squared_order < 0) then
      write (*, '(a)') &
           'ERROR in decay_chain_scales: decay QCD order exceeds total order'
      stop 1
    end if
  end function production_qcd_squared_order


  double precision function decay_qcd_coupling_weight(qcd_power, &
                                                       factor_indices, node_scales)
    integer, intent(in), optional :: qcd_power
    integer, intent(in), optional :: factor_indices(:)
    double precision, intent(in), optional :: node_scales(:)
    integer :: node, qcd_order, total_power, corrected_power
    integer :: factor_index
    double precision :: coupling, scale
    double precision, parameter :: pi = 3.14159265358979323846d0

    decay_qcd_coupling_weight = 1d0
    if (has_nlo_decay()) then
      total_power = decay_qcd_squared_order()
      if (present(qcd_power)) total_power = qcd_power
      corrected_power = total_power - nlo_decay_static_qcd_order()
      if (corrected_power < nlo_decay_born_qcd_order()) then
        call fail_scales('corrected-decay QCD power is below its Born order')
      end if
      do node = 1, nlo_decay_node_count()
        if (node == nlo_decay_corrected_node()) then
          qcd_order = corrected_power
        else
          qcd_order = 2*nlo_decay_node_qcd_order(node)
        end if
        if (qcd_order == 0) cycle
        factor_index = selected_factor_index(&
             nlo_decay_node_pdg(node), factor_indices)
        scale = decay_renormalization_scale(&
             nlo_decay_node_pdg(node), factor_index)
        if (present(node_scales)) &
             scale = node_scales(node)*decay_scale_factor(factor_index)
        coupling = sqrt(4d0*pi*alphas(scale))
        decay_qcd_coupling_weight = decay_qcd_coupling_weight* &
             coupling**qcd_order
      end do
      return
    end if
    if (.not. has_decay_chains()) return
    do node = 1, decay_node_count()
      qcd_order = node_qcd_order(node)
      if (qcd_order == 0) cycle
      factor_index = selected_factor_index(node_pdg(node), factor_indices)
      scale = decay_renormalization_scale(node_pdg(node), factor_index)
      if (present(node_scales)) &
           scale = node_scales(node)*decay_scale_factor(factor_index)
      coupling = sqrt(4d0*pi*alphas(scale))
      decay_qcd_coupling_weight = decay_qcd_coupling_weight* &
           coupling**(2*qcd_order)
    end do
  end function decay_qcd_coupling_weight


  integer function nlo_decay_static_qcd_order()
    integer :: node
    nlo_decay_static_qcd_order = 0
    do node = 1, nlo_decay_node_count()
      if (node == nlo_decay_corrected_node()) cycle
      nlo_decay_static_qcd_order = nlo_decay_static_qcd_order + &
           2*nlo_decay_node_qcd_order(node)
    end do
  end function nlo_decay_static_qcd_order


  integer function nlo_decay_born_total_qcd_order()
    nlo_decay_born_total_qcd_order = nlo_decay_static_qcd_order() + &
         nlo_decay_born_qcd_order()
  end function nlo_decay_born_total_qcd_order


  double precision function decay_qcd_coupling_rescaling(production_g, &
                                                          qcd_power, &
                                                          factor_indices, node_scales)
    double precision, intent(in) :: production_g
    integer, intent(in), optional :: qcd_power
    integer, intent(in), optional :: factor_indices(:)
    double precision, intent(in), optional :: node_scales(:)
    integer :: power

    power = decay_qcd_squared_order()
    if (present(qcd_power)) power = qcd_power
    if (power == 0) then
      decay_qcd_coupling_rescaling = 1d0
      return
    end if
    if (production_g <= 0d0) then
      write (*, '(a)') &
           'ERROR in decay_chain_scales: production coupling is not positive'
      stop 1
    end if
    decay_qcd_coupling_rescaling = &
         decay_qcd_coupling_weight(power, factor_indices, node_scales)/production_g**power
  end function decay_qcd_coupling_rescaling


  integer function selected_factor_index(pdg, factor_indices)
    integer, intent(in) :: pdg
    integer, intent(in), optional :: factor_indices(:)
    integer :: species_index

    selected_factor_index = 1
    if (.not. present(factor_indices)) return
    if (size(factor_indices) == 0) return
    if (size(factor_indices) /= decay_scale_species_count()) then
      call fail_scales('decay factor-index array has the wrong size')
    end if
    species_index = decay_scale_species_index(pdg)
    if (species_index > 0) then
      selected_factor_index = factor_indices(species_index)
    end if
  end function selected_factor_index


  integer function corrected_born_qcd_squared_order(total_nlo_qcd_order)
    integer, intent(in) :: total_nlo_qcd_order
    if (has_nlo_decay()) then
      corrected_born_qcd_squared_order = nlo_decay_born_qcd_order()
    else
      corrected_born_qcd_squared_order = &
           production_qcd_squared_order(total_nlo_qcd_order - 2)
    end if
  end function corrected_born_qcd_squared_order


  subroutine select_production_core_momenta(visible_momenta, configuration, &
                                            core_momenta)
    double precision, intent(in) :: visible_momenta(0:3, nexternal)
    integer, intent(in) :: configuration
    double precision, intent(out) :: core_momenta(0:3, nexternal)

    if (has_nlo_decay()) then
      call get_nlo_decay_production_momenta(core_momenta)
    else if (has_decay_chains()) then
      call contract_visible_momenta(context_for_fks(configuration), &
                                    visible_momenta, core_momenta)
    else
      core_momenta = visible_momenta
    end if
  end subroutine select_production_core_momenta


  subroutine select_production_ren_scale_momenta(visible_momenta, &
                                                 configuration, &
                                                 scale_momenta)
    double precision, intent(in) :: visible_momenta(0:3, nexternal)
    integer, intent(in) :: configuration
    double precision, intent(out) :: scale_momenta(0:3, nexternal)

    if ((has_decay_chains() .or. has_nlo_decay()) .and. &
        use_decayed_production_ren_scale_momenta()) then
      scale_momenta = visible_momenta
    else
      call select_production_core_momenta(visible_momenta, configuration, &
                                          scale_momenta)
    end if
  end subroutine select_production_ren_scale_momenta


  subroutine fail_scales(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in decay_chain_scales: '//trim(message)
    stop 1
  end subroutine fail_scales

end module decay_chain_scales
