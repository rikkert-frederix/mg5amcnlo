module multiplicative_scale_state
  use process_dimensions, only: nexternal, validate_process_dimensions
  use alfas_functions_module, only: alphas
  use factorized_phase_space, only: factorized_radiation_state, &
       fetch_factorized_radiation_state
  use decay_chain_parameters, only: decay_renormalization_scale, &
       use_decayed_production_ren_scale_momenta
  use run_state, only: scale, q2fact, mur2_current, muf12_current, &
       muf22_current, qes2_current, fixed_ren_scale, fixed_fac_scale, &
       mur_over_ref, muf1_over_ref, mur_ref_fixed, muf1_ref_fixed, &
       dynamical_scale_choice
  use fixed_order_user_hooks, only: fixed_user_scale
  use multiplicative_kinematics, only: &
       factorized_production_scale_sums, factorized_visible_scale_sums
  use fnlo_process_common, only: soft_counterevent, real_event
  use multiplicative_generated_metadata, only: &
       multiplicative_block_count, multiplicative_physical_blocks, &
       multiplicative_block_pdgs, multiplicative_born_qcd_powers, &
       multiplicative_component_position
  implicit none
  private

  double precision, parameter :: pi = 3.14159265358979323846d0
  double precision, allocatable, save :: reference_scale_squared(:)
  double precision, allocatable, save :: reference_coupling(:)
  logical, allocatable, save :: reference_is_valid(:)
  integer(kind=8), save :: reference_generation = 0_8
  integer(kind=8), save :: active_coupling_context = 0_8

  public :: initialize_multiplicative_scale_references
  public :: activate_multiplicative_block_reference
  public :: multiplicative_reference_scale_squared
  public :: multiplicative_reference_coupling
  public :: multiplicative_active_coupling_context
  public :: multiplicative_block_scale_logarithm
  public :: multiplicative_block_coupling_rescaling
  public :: build_multiplicative_scale_tables

  interface
    subroutine set_model_loop_ren_scale_bridge(mur, g_value)
      double precision, intent(in) :: mur, g_value
    end subroutine set_model_loop_ren_scale_bridge

    subroutine set_model_qes_scale_bridge(qes_squared)
      double precision, intent(in) :: qes_squared
    end subroutine set_model_qes_scale_bridge

  end interface

contains

  subroutine initialize_multiplicative_scale_references()
    type(factorized_radiation_state) :: radiation
    integer :: position, block, block_count, pdg
    double precision :: reference_scale
    logical :: available

    call validate_process_dimensions()
    if (.not. allocated(reference_scale_squared)) then
      allocate(reference_scale_squared(0:nexternal))
      allocate(reference_coupling(0:nexternal))
      allocate(reference_is_valid(0:nexternal))
    end if
    reference_scale_squared = 0d0
    reference_coupling = 0d0
    reference_is_valid = .false.
    reference_generation = reference_generation + 1_8
    if (reference_generation <= 0_8) then
      call fail_multiplicative_scale_state( &
           'the coupling-reference generation counter overflowed')
    end if
    active_coupling_context = 0_8

    block_count = multiplicative_block_count
    if (block_count < 1) then
      call fail_multiplicative_scale_state( &
           'the generated density graph has no blocks')
    end if
    do position = 1, block_count
      block = multiplicative_physical_blocks(position)
      call validate_block(block)
      if (reference_is_valid(block)) then
        call fail_multiplicative_scale_state( &
             'two generated components own one physical block')
      end if
      if (block == 0) then
        call fetch_production_reference_radiation(radiation, available)
        if (.not. available) then
          call fail_multiplicative_scale_state( &
               'the production block has no positive reference scale')
        end if
        reference_scale = max(2d0, radiation%sqrt_shat)
      else
        pdg = multiplicative_block_pdgs(position)
        if (pdg == 0) then
          call fail_multiplicative_scale_state( &
               'a decay block has no generated parent PDG')
        end if
        reference_scale = decay_renormalization_scale(pdg)
      end if
      if (reference_scale <= 0d0) then
        call fail_multiplicative_scale_state( &
             'a block reference scale is not positive')
      end if
      reference_scale_squared(block) = reference_scale**2
      reference_coupling(block) = &
           sqrt(4d0*pi*alphas(reference_scale))
      if (reference_coupling(block) <= 0d0) then
        call fail_multiplicative_scale_state( &
             'a block reference coupling is not positive')
      end if
      reference_is_valid(block) = .true.
    end do
  end subroutine initialize_multiplicative_scale_references


  subroutine fetch_production_reference_radiation(radiation, available)
    type(factorized_radiation_state), intent(out) :: radiation
    logical, intent(out) :: available
    type(factorized_radiation_state) :: candidate
    integer :: event_slot
    logical :: candidate_available

    radiation = factorized_radiation_state()
    available = .false.
    ! Massive FKS sectors may intentionally omit one or more mapped
    ! counterevents.  The reference is shared by every atom in this block,
    ! so any stored positive production scale is a valid reference; prefer
    ! the Born-like slots and fall back to the real event.
    do event_slot = soft_counterevent, real_event
      call fetch_factorized_radiation_state( &
           event_slot, 0, candidate, candidate_available)
      if (.not. candidate_available .or. candidate%sqrt_shat <= 0d0) cycle
      radiation = candidate
      available = .true.
      return
    end do
  end subroutine fetch_production_reference_radiation


  subroutine activate_multiplicative_block_reference(block)
    integer, intent(in) :: block
    double precision :: reference_scale
    integer(kind=8) :: target_coupling_context

    call require_reference(block)
    target_coupling_context = reference_generation* &
         int(nexternal + 2, kind=8) + int(block + 1, kind=8)
    ! A coalesced kinematic family can contain many primitives from the same
    ! block.  Their immutable coupling reference is identical, so avoid
    ! repeating UPDATE_AS_PARAM (including the loop-only UV/R2 refresh) for
    ! every descriptor in that family.
    if (active_coupling_context == target_coupling_context) return
    reference_scale = sqrt(reference_scale_squared(block))

    ! Every primitive in a block is generated at this immutable reference.
    ! In particular, the loop provider sees MU_R=Q and the same coupling as
    ! the tree/integrated operators in that block.  Desired physical scales
    ! are introduced later by the explicit logarithms and coupling powers.
    scale = reference_scale
    q2fact = reference_scale_squared(block)
    mur2_current = reference_scale_squared(block)
    muf12_current = reference_scale_squared(block)
    muf22_current = reference_scale_squared(block)
    qes2_current = reference_scale_squared(block)
    call set_model_qes_scale_bridge(reference_scale_squared(block))
    call set_model_loop_ren_scale_bridge( &
         reference_scale, reference_coupling(block))
    active_coupling_context = target_coupling_context
  end subroutine activate_multiplicative_block_reference


  integer(kind=8) function multiplicative_active_coupling_context()
    multiplicative_active_coupling_context = active_coupling_context
  end function multiplicative_active_coupling_context


  double precision function multiplicative_reference_scale_squared(block)
    integer, intent(in) :: block

    call require_reference(block)
    multiplicative_reference_scale_squared = &
         reference_scale_squared(block)
  end function multiplicative_reference_scale_squared


  double precision function multiplicative_reference_coupling(block)
    integer, intent(in) :: block

    call require_reference(block)
    multiplicative_reference_coupling = reference_coupling(block)
  end function multiplicative_reference_coupling


  double precision function multiplicative_block_scale_logarithm( &
       block, physical_scale_squared)
    integer, intent(in) :: block
    double precision, intent(in) :: physical_scale_squared

    call require_reference(block)
    if (physical_scale_squared <= 0d0) then
      call fail_multiplicative_scale_state( &
           'a physical scale squared is not positive')
    end if
    multiplicative_block_scale_logarithm = &
         log(physical_scale_squared/reference_scale_squared(block))
  end function multiplicative_block_scale_logarithm


  double precision function multiplicative_block_coupling_rescaling( &
       block, nlo_order, physical_scale)
    integer, intent(in) :: block, nlo_order
    double precision, intent(in) :: physical_scale
    integer :: qcd_power
    double precision :: physical_coupling

    call require_reference(block)
    if (nlo_order < 0 .or. nlo_order > 1) then
      call fail_multiplicative_scale_state( &
           'one block term has an invalid NLO order')
    end if
    if (physical_scale <= 0d0) then
      call fail_multiplicative_scale_state( &
           'a physical renormalization scale is not positive')
    end if
    qcd_power = multiplicative_born_qcd_powers( &
         multiplicative_component_position(block)) + 2*nlo_order
    if (qcd_power < 0) then
      call fail_multiplicative_scale_state( &
           'a block has a negative QCD coupling power')
    end if
    if (qcd_power == 0) then
      multiplicative_block_coupling_rescaling = 1d0
      return
    end if
    physical_coupling = sqrt(4d0*pi*alphas(physical_scale))
    multiplicative_block_coupling_rescaling = &
         (physical_coupling/reference_coupling(block))**qcd_power
  end function multiplicative_block_coupling_rescaling


  subroutine build_multiplicative_scale_tables( &
       event_slots, logarithmic_mu2_r, logarithmic_mu2_f, &
       coupling_rescaling, production_mu2_r, production_mu2_f, &
       production_ren_factor, production_fac_factor, &
       decay_factor_indices, production_dynamic_choice)
    integer, intent(in) :: event_slots(0:)
    double precision, intent(out) :: logarithmic_mu2_r(0:)
    double precision, intent(out) :: logarithmic_mu2_f(0:)
    double precision, intent(out) :: coupling_rescaling(0:, 0:)
    double precision, intent(out) :: production_mu2_r, production_mu2_f
    double precision, intent(in), optional :: production_ren_factor
    double precision, intent(in), optional :: production_fac_factor
    integer, intent(in), optional :: decay_factor_indices(0:)
    integer, intent(in), optional :: production_dynamic_choice
    integer :: block_count, position, block, pdg, factor_index, nlo_order
    double precision :: mu_r, mu_f, ren_factor, fac_factor

    if (ubound(event_slots, 1) < nexternal .or. &
        ubound(logarithmic_mu2_r, 1) < nexternal .or. &
        ubound(logarithmic_mu2_f, 1) < nexternal .or. &
        ubound(coupling_rescaling, 1) < nexternal .or. &
        ubound(coupling_rescaling, 2) < 1) then
      call fail_multiplicative_scale_state( &
           'a multiplicative scale table has the wrong shape')
    end if
    if (present(decay_factor_indices)) then
      if (ubound(decay_factor_indices, 1) < nexternal) then
        call fail_multiplicative_scale_state( &
             'a decay-scale factor table has the wrong shape')
      end if
    end if
    ren_factor = 1d0
    fac_factor = 1d0
    if (present(production_ren_factor)) ren_factor = production_ren_factor
    if (present(production_fac_factor)) fac_factor = production_fac_factor
    if (ren_factor <= 0d0 .or. fac_factor <= 0d0) then
      call fail_multiplicative_scale_state( &
           'a production scale multiplier is not positive')
    end if

    call production_physical_scales( &
         event_slots, mu_r, mu_f, production_dynamic_choice)
    mu_r = mu_r*ren_factor
    mu_f = mu_f*fac_factor
    production_mu2_r = mu_r**2
    production_mu2_f = mu_f**2
    logarithmic_mu2_r = 0d0
    logarithmic_mu2_f = 0d0
    coupling_rescaling = 1d0

    block_count = multiplicative_block_count
    do position = 1, block_count
      block = multiplicative_physical_blocks(position)
      if (block == 0) then
        mu_r = sqrt(production_mu2_r)
        mu_f = sqrt(production_mu2_f)
      else
        factor_index = 1
        if (present(decay_factor_indices)) then
          factor_index = decay_factor_indices(block)
        end if
        pdg = multiplicative_block_pdgs(position)
        mu_r = decay_renormalization_scale(pdg, factor_index)
        ! Decay blocks carry no incoming-PDF scale.  Using their own
        ! renormalization scale for the formally present mu_F slot makes the
        ! zero coefficient explicit and keeps one uniform table interface.
        mu_f = mu_r
      end if
      logarithmic_mu2_r(block) = &
           multiplicative_block_scale_logarithm(block, mu_r**2)
      logarithmic_mu2_f(block) = &
           multiplicative_block_scale_logarithm(block, mu_f**2)
      do nlo_order = 0, 1
        coupling_rescaling(block, nlo_order) = &
             multiplicative_block_coupling_rescaling( &
             block, nlo_order, mu_r)
      end do
    end do
  end subroutine build_multiplicative_scale_tables


  subroutine production_physical_scales( &
       event_slots, mu_r, mu_f, dynamic_choice)
    integer, intent(in) :: event_slots(0:)
    double precision, intent(out) :: mu_r, mu_f
    integer, intent(in), optional :: dynamic_choice
    double precision :: core_sum_et, core_sum_mt
    double precision :: ren_sum_et, ren_sum_mt, reference

    call factorized_production_scale_sums( &
         event_slots, core_sum_et, core_sum_mt)
    ren_sum_et = core_sum_et
    ren_sum_mt = core_sum_mt
    if (use_decayed_production_ren_scale_momenta()) then
      call factorized_visible_scale_sums( &
           event_slots, ren_sum_et, ren_sum_mt)
    end if
    if (fixed_ren_scale) then
      reference = mur_ref_fixed
    else
      reference = max(2d0, production_dynamic_reference( &
           ren_sum_et, ren_sum_mt, dynamic_choice))
    end if
    mu_r = mur_over_ref*reference
    if (fixed_fac_scale) then
      reference = muf1_ref_fixed
    else
      reference = max(2d0, production_dynamic_reference( &
           core_sum_et, core_sum_mt, dynamic_choice))
    end if
    mu_f = muf1_over_ref*reference
    if (mu_r <= 0d0 .or. mu_f <= 0d0) then
      call fail_multiplicative_scale_state( &
           'a physical production scale is not positive')
    end if
  end subroutine production_physical_scales


  double precision function production_dynamic_reference( &
       sum_et, sum_mt, requested_choice)
    double precision, intent(in) :: sum_et, sum_mt
    integer, intent(in), optional :: requested_choice
    character(len=80) :: scale_id
    integer :: choice

    choice = dynamical_scale_choice
    if (present(requested_choice)) choice = requested_choice
    select case (choice)
    case (1)
      production_dynamic_reference = sum_et
    case (2)
      production_dynamic_reference = sum_mt
    case (3, -1)
      production_dynamic_reference = sum_mt/2d0
    case (-2)
      production_dynamic_reference = mur_ref_fixed
    case (10, 0)
      production_dynamic_reference = &
           fixed_user_scale(mur_ref_fixed, scale_id)
    case default
      call fail_multiplicative_scale_state( &
           'the production dynamical-scale choice is unknown')
    end select
  end function production_dynamic_reference


  subroutine require_reference(block)
    integer, intent(in) :: block

    call validate_block(block)
    if (.not. allocated(reference_is_valid) .or. &
        .not. reference_is_valid(block)) then
      call fail_multiplicative_scale_state( &
           'a block scale reference is uninitialized')
    end if
  end subroutine require_reference


  subroutine validate_block(block)
    integer, intent(in) :: block

    if (block < 0 .or. block > nexternal) then
      call fail_multiplicative_scale_state( &
           'a physical block index is out of range')
    end if
  end subroutine validate_block


  subroutine fail_multiplicative_scale_state(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_scale_state: '//trim(message)
    stop 1
  end subroutine fail_multiplicative_scale_state
end module multiplicative_scale_state
