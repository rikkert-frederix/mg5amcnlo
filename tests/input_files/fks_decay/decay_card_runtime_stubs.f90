! A two-resonance fixture for testing the real fNLO card/scale modules.
module process_dimensions
  implicit none
  integer, parameter :: nexternal = 8
contains
  subroutine validate_process_dimensions()
  end subroutine
end module

module fnlo_process_common
  integer, parameter :: soft_counterevent = 0, real_event = 3
end module

module nlo_contribution_bundle
  logical :: is_bundle = .true.
  integer :: active = 1
contains
  logical function has_nlo_contribution_bundle()
    has_nlo_contribution_bundle = is_bundle
  end function
  logical function bundle_species_is_nlo(pdg)
    integer, intent(in) :: pdg
    bundle_species_is_nlo = abs(pdg) == 6
  end function
  logical function contribution_is_nlo_decay(contribution)
    integer, intent(in) :: contribution
    contribution_is_nlo_decay = contribution > 1
  end function
  integer function contribution_parent_pdg(contribution)
    integer, intent(in) :: contribution
    contribution_parent_pdg = 6
  end function
  integer function active_nlo_contribution()
    active_nlo_contribution = active
  end function
end module

module nlo_decay_metadata
  use nlo_contribution_bundle, only: active
  integer :: born_qcd = 0
contains
  logical function has_nlo_decay()
    has_nlo_decay = active > 1
  end function
  integer function corrected_parent_pdg()
    corrected_parent_pdg = 6
  end function
  integer function nlo_decay_node_count()
    nlo_decay_node_count = 2
  end function
  integer function nlo_decay_node_pdg(node)
    integer, intent(in) :: node
    nlo_decay_node_pdg = 6*(-1)**(node-1)
  end function
  integer function nlo_decay_node_qcd_order(node)
    integer, intent(in) :: node
    nlo_decay_node_qcd_order = born_qcd
  end function
  integer function nlo_decay_production_born_qcd_order()
    nlo_decay_production_born_qcd_order = 4
  end function
  integer function nlo_decay_born_qcd_order()
    nlo_decay_born_qcd_order = 2*born_qcd
  end function
  integer function nlo_decay_corrected_node()
    nlo_decay_corrected_node = 1
  end function
end module

module decay_chain_metadata
  use nlo_decay_metadata, only: born_qcd
  use nlo_contribution_bundle, only: active
contains
  logical function has_decay_chains()
    has_decay_chains = active == 1
  end function
  integer function decay_node_count()
    decay_node_count = 2
  end function
  integer function node_pdg(node)
    integer, intent(in) :: node
    node_pdg = 6*(-1)**(node-1)
  end function
  integer function node_qcd_order(node)
    integer, intent(in) :: node
    node_qcd_order = born_qcd
  end function
  integer function context_for_fks(configuration)
    integer, intent(in) :: configuration
    context_for_fks = 1
  end function
end module

module alfas_functions_module
contains
  double precision function alphas(mu)
    double precision, intent(in) :: mu
    alphas = 0.1d0/(1d0 + 0.1d0*log(mu/100d0))
  end function
end module

module decay_chain_kinematics
contains
  subroutine contract_visible_momenta(context, visible, core)
    integer, intent(in) :: context
    double precision, intent(in) :: visible(0:3, 8)
    double precision, intent(out) :: core(0:3, 8)
    core = visible
  end subroutine
end module

module nlo_decay_kinematics
contains
  subroutine get_nlo_decay_production_momenta(core)
    double precision, intent(out) :: core(0:3, 8)
    core = 0d0
  end subroutine
end module
