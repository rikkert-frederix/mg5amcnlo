module decay_chain_scales
  use process_dimensions, only: nexternal
  use decay_chain_metadata, only: has_decay_chains, decay_node_count, &
       context_for_fks, node_pdg, node_qcd_order
  use decay_chain_parameters, only: decay_renormalization_scale, &
       use_decayed_production_ren_scale_momenta
  use decay_chain_kinematics, only: contract_visible_momenta
  use nlo_decay_metadata, only: has_nlo_decay, corrected_parent_pdg, &
       nlo_decay_production_born_qcd_order, nlo_decay_born_qcd_order
  use nlo_decay_kinematics, only: get_nlo_decay_production_momenta
  use alfas_functions_module, only: alphas
  implicit none
  private

  public :: decay_qcd_squared_order, production_qcd_squared_order
  public :: decay_qcd_coupling_weight
  public :: decay_qcd_coupling_rescaling
  public :: corrected_born_qcd_squared_order
  public :: select_production_core_momenta
  public :: select_production_ren_scale_momenta

contains

  integer function decay_qcd_squared_order(total_qcd_order)
    integer, intent(in), optional :: total_qcd_order
    integer :: node

    decay_qcd_squared_order = 0
    if (has_nlo_decay()) then
      if (present(total_qcd_order)) then
        decay_qcd_squared_order = total_qcd_order - &
             nlo_decay_production_born_qcd_order()
      else
        decay_qcd_squared_order = nlo_decay_born_qcd_order()
      end if
      if (decay_qcd_squared_order < nlo_decay_born_qcd_order()) then
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
          nlo_decay_born_qcd_order()) then
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


  double precision function decay_qcd_coupling_weight(qcd_power)
    integer, intent(in), optional :: qcd_power
    integer :: node, qcd_order
    double precision :: coupling, scale
    double precision, parameter :: pi = 3.14159265358979323846d0

    decay_qcd_coupling_weight = 1d0
    if (has_nlo_decay()) then
      qcd_order = decay_qcd_squared_order()
      if (present(qcd_power)) qcd_order = qcd_power
      if (qcd_order == 0) return
      scale = decay_renormalization_scale(corrected_parent_pdg())
      coupling = sqrt(4d0*pi*alphas(scale))
      decay_qcd_coupling_weight = coupling**qcd_order
      return
    end if
    if (.not. has_decay_chains()) return
    do node = 1, decay_node_count()
      qcd_order = node_qcd_order(node)
      if (qcd_order == 0) cycle
      scale = decay_renormalization_scale(node_pdg(node))
      coupling = sqrt(4d0*pi*alphas(scale))
      decay_qcd_coupling_weight = decay_qcd_coupling_weight* &
           coupling**(2*qcd_order)
    end do
  end function decay_qcd_coupling_weight


  double precision function decay_qcd_coupling_rescaling(production_g, &
                                                          qcd_power)
    double precision, intent(in) :: production_g
    integer, intent(in), optional :: qcd_power
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
         decay_qcd_coupling_weight(power)/production_g**power
  end function decay_qcd_coupling_rescaling


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
