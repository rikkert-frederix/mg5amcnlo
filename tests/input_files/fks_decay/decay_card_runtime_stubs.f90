! A two-resonance fixture for testing the real fNLO card/scale modules.
module process_dimensions
  implicit none
  integer, parameter :: nexternal = 8, nincoming = 2
contains
  subroutine validate_process_dimensions()
  end subroutine
end module

module fnlo_process_common
  integer, parameter :: soft_counterevent = 0, real_event = 3
  integer :: nfksprocess = 1
  logical :: calculated_born = .false.
  character(len=80) :: mur_id_str, muf1_id_str, muf2_id_str, qes_id_str, temp_scale_id
end module

module production_scale_fixture
  integer :: core_count = 7, core_pdgs(8) = [2, -1, 6, -6, -11, 12, 21, 0]
  double precision :: production_p(0:3, 8) = 0d0
end module

module timing_state
  double precision :: t_coupl = 0d0
end module

module kin_functions_module
contains
  double precision function et_impl(p)
    double precision, intent(in) :: p(0:3)
    et_impl = sqrt(p(1)**2+p(2)**2)
  end function
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
    if (contribution == 3) contribution_parent_pdg = -6
  end function
  integer function active_nlo_contribution()
    active_nlo_contribution = active
  end function
end module

module nlo_decay_metadata
  use nlo_contribution_bundle, only: active
  use production_scale_fixture
  integer :: born_qcd = 0
contains
  integer function nlo_decay_production_count()
    nlo_decay_production_count = core_count
  end function
  integer function nlo_decay_production_pdg(leg)
    integer, intent(in) :: leg
    nlo_decay_production_pdg = core_pdgs(leg)
  end function
  logical function has_nlo_decay()
    has_nlo_decay = active > 1
  end function
  integer function corrected_parent_pdg()
    corrected_parent_pdg = 6
    if (active == 3) corrected_parent_pdg = -6
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
    if (active == 3) nlo_decay_corrected_node = 2
  end function
end module

module extra_weights
  integer :: dyn_scale(0:1) = [1, 3]
  logical :: lscalevar(1) = .true.
  double precision :: scalevarR(0:3) = [3d0, 1d0, .5d0, 2d0]
  double precision :: scalevarF(0:3) = [3d0, 1d0, .5d0, 2d0]
end module

module run_state
  logical :: do_rwgt_scale = .true., do_rwgt_decay_scale = .false.
  logical :: fixed_ren_scale = .false., fixed_fac_scale = .false., fixed_qes_scale = .false.
  integer :: dynamical_scale_choice = 3
  double precision :: scale, mur2_current, muf12_current, muf22_current, qes2_current, q2fact(2)
  double precision :: mur_over_ref = 1d0, muf1_over_ref = 1d0, muf2_over_ref = 1d0, qes_over_ref = 1d0
  double precision :: mur_ref_fixed = 100d0, muf1_ref_fixed = 100d0, muf2_ref_fixed = 100d0
  double precision :: qes_ref_fixed = 100d0
end module

! The real weight-line module is tested; its unrelated aggregation interfaces
! below are link stubs, not implementations of the physics under test.
module spin_density_matrix_results
  integer, parameter :: spin_density_bornlike_branch = 1, spin_density_real_branch = 2
end module

module multiplicative_nlo_decay
  type multiplicative_nlo_workspace
    integer :: component_count, weight_count
    integer, allocatable :: component_open_sizes(:)
  end type
contains
  subroutine add_multiplicative_block_density(workspace, component, branch, density)
    type(multiplicative_nlo_workspace), intent(inout) :: workspace
    integer, intent(in) :: component, branch
    complex(kind=8), intent(in) :: density(:, :, :)
    stop 'aggregation is not part of this fixture'
  end subroutine
end module

integer function sdm_branch_max_open_size()
  sdm_branch_max_open_size = 1
end function

module decay_chain_metadata
  use nlo_decay_metadata, only: born_qcd
  use nlo_contribution_bundle, only: active
  use production_scale_fixture
contains
  integer function context_core_count(context)
    integer, intent(in) :: context
    context_core_count = core_count
  end function
  integer function core_leg_pdg(context, leg)
    integer, intent(in) :: context, leg
    core_leg_pdg = core_pdgs(leg)
  end function
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
  use production_scale_fixture
contains
  subroutine get_nlo_decay_production_momenta(core)
    double precision, intent(out) :: core(0:3, 8)
    core = production_p
  end subroutine
end module

! Model side effects are not needed to check actual production scale routing.
subroutine set_model_ren_scale_bridge(mur, g_value)
  double precision, intent(in) :: mur, g_value
end subroutine
subroutine set_model_qes_scale_bridge(qes_squared)
  double precision, intent(in) :: qes_squared
end subroutine
subroutine update_model_momenta_bridge(p, reset_momenta)
  double precision, intent(in) :: p(0:3, *)
  logical, intent(in) :: reset_momenta
end subroutine
