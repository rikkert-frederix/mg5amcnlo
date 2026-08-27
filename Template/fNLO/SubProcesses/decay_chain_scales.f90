module decay_chain_scales
  use process_dimensions, only: nexternal
  use decay_chain_metadata, only: has_decay_chains, decay_node_count, &
       context_for_fks, node_pdg, node_qcd_order
  use decay_chain_parameters, only: decay_renormalization_scale, &
       use_decayed_production_ren_scale_momenta
  use decay_chain_kinematics, only: contract_visible_momenta
  use alfas_functions_module, only: alphas
  implicit none
  private

  public :: decay_qcd_squared_order, production_qcd_squared_order
  public :: decay_qcd_coupling_weight
  public :: decay_qcd_coupling_rescaling
  public :: select_production_core_momenta
  public :: select_production_ren_scale_momenta

contains

  integer function decay_qcd_squared_order()
    integer :: node

    decay_qcd_squared_order = 0
    if (.not. has_decay_chains()) return
    do node = 1, decay_node_count()
      decay_qcd_squared_order = decay_qcd_squared_order + &
           2*node_qcd_order(node)
    end do
  end function decay_qcd_squared_order


  integer function production_qcd_squared_order(total_qcd_order)
    integer, intent(in) :: total_qcd_order

    production_qcd_squared_order = &
         total_qcd_order - decay_qcd_squared_order()
    if (production_qcd_squared_order < 0) then
      write (*, '(a)') &
           'ERROR in decay_chain_scales: decay QCD order exceeds total order'
      stop 1
    end if
  end function production_qcd_squared_order


  double precision function decay_qcd_coupling_weight()
    integer :: node, qcd_order
    double precision :: coupling, scale
    double precision, parameter :: pi = 3.14159265358979323846d0

    decay_qcd_coupling_weight = 1d0
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


  double precision function decay_qcd_coupling_rescaling(production_g)
    double precision, intent(in) :: production_g
    integer :: power

    power = decay_qcd_squared_order()
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
         decay_qcd_coupling_weight()/production_g**power
  end function decay_qcd_coupling_rescaling


  subroutine select_production_core_momenta(visible_momenta, configuration, &
                                            core_momenta)
    double precision, intent(in) :: visible_momenta(0:3, nexternal)
    integer, intent(in) :: configuration
    double precision, intent(out) :: core_momenta(0:3, nexternal)

    if (has_decay_chains()) then
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

    if (has_decay_chains() .and. &
        use_decayed_production_ren_scale_momenta()) then
      scale_momenta = visible_momenta
    else
      call select_production_core_momenta(visible_momenta, configuration, &
                                          scale_momenta)
    end if
  end subroutine select_production_ren_scale_momenta

end module decay_chain_scales
