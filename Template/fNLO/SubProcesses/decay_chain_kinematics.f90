module decay_chain_kinematics
  use process_dimensions, only: nexternal, nincoming, validate_process_dimensions
  use fnlo_process_common, only: event_from_decay => from_decay, &
       soft_counterevent
  use factorized_phase_space, only: factorized_measure_state, &
       store_factorized_block_momenta, store_factorized_kernel_momenta, &
       store_factorized_base_measure, compose_factorized_base_measure
  use factorized_block_kinematics, only: &
       generate_nbody => generate_factorized_nbody, &
       generate_nbody_rest => generate_factorized_nbody_rest, &
       boost_nbody_from_rest => boost_factorized_block_from_rest, &
       boost_from_rest => boost_factorized_momentum_from_rest, &
       minkowski_square => factorized_minkowski_square
  use phase_space_kinematics, only: phase_space_lambda
  use decay_chain_metadata, only: has_decay_chains, decay_node_count, &
       decay_leaf_count, decay_random_dimension, born_context, &
       context_for_fks, context_core_count, &
       core_target_kind, core_target_id, core_leg_pdg, &
       direct_leg_target, decay_node_target, leaf_visible_leg, &
       node_pdg, node_child_count, node_child_kind, node_child_id, &
       decay_leaf_child, decay_node_child, leaf_pdg, &
       visible_color_pair
  use decay_chain_parameters, only: decay_dummy_width_ratio, &
       decay_physical_width
  implicit none
  private

  logical, save :: initialized = .false.
  double precision, allocatable, save :: node_masses(:)
  double precision, allocatable, save :: leaf_masses(:)
  double precision, allocatable, save :: core_born_storage(:, :)
  double precision, allocatable, save :: node_rest_storage(:, :, :)
  logical, allocatable, save :: node_rest_valid(:)

  public :: initialize_decay_chain_kinematics
  public :: minimum_core_final_mass, core_mass
  public :: generate_core_born_and_decays, expand_real_decay_momenta
  public :: store_core_event_momenta
  public :: get_core_born_momenta, get_core_mass_buffer
  public :: active_core_count, fks_leg_mass
  public :: map_core_color_pair
  public :: contract_visible_momenta
  ! These Lorentz-covariant building blocks are also used by the dedicated
  ! NLO-decay phase-space path.  They do not depend on decay-chain metadata.
  public :: generate_nbody, generate_nbody_rest
  public :: boost_from_rest, boost_nbody_from_rest, minkowski_square
  public :: set_decay_cut_mask

  interface
    double precision function get_mass_from_id(id)
      integer, intent(in) :: id
    end function get_mass_from_id
  end interface

contains

  subroutine initialize_decay_chain_kinematics()
    integer :: node, leaf

    if (initialized) return
    call validate_process_dimensions()
    if (.not. has_decay_chains()) then
      return
    end if

    allocate(node_masses(decay_node_count()))
    allocate(leaf_masses(decay_leaf_count()))
    allocate(core_born_storage(0:3, nexternal - 1))
    allocate(node_rest_storage(0:3, nexternal, decay_node_count()))
    allocate(node_rest_valid(decay_node_count()))
    core_born_storage = 0d0
    node_rest_storage = 0d0
    node_rest_valid = .false.

    do node = 1, decay_node_count()
      node_masses(node) = abs(get_mass_from_id(node_pdg(node)))
      if (node_masses(node) <= 0d0) then
        call fail_kinematics('a forced decay parent has zero model mass')
      end if
    end do
    do leaf = 1, decay_leaf_count()
      leaf_masses(leaf) = abs(get_mass_from_id(leaf_pdg(leaf)))
    end do
    initialized = .true.
  end subroutine initialize_decay_chain_kinematics


  integer function production_random_dimension()
    integer :: final_count
    call require_enabled()
    final_count = context_core_count(born_context()) - nincoming
    if (final_count == 1) then
      production_random_dimension = 0
    else
      production_random_dimension = 3*final_count - 4
    end if
  end function production_random_dimension


  integer function decay_variable_start()
    call require_enabled()
    decay_variable_start = production_random_dimension() + 1
  end function decay_variable_start


  double precision function minimum_core_final_mass()
    integer :: context, leg
    call require_enabled()
    context = born_context()
    minimum_core_final_mass = 0d0
    do leg = nincoming + 1, context_core_count(context)
      minimum_core_final_mass = minimum_core_final_mass + &
           core_mass(context, leg)
    end do
  end function minimum_core_final_mass


  double precision function core_mass(context, leg)
    integer, intent(in) :: context, leg
    integer :: target
    call require_enabled()
    if (core_target_kind(context, leg) == decay_node_target) then
      target = core_target_id(context, leg)
      core_mass = node_masses(target)
    else
      core_mass = abs(get_mass_from_id(core_leg_pdg(context, leg)))
    end if
  end function core_mass


  subroutine generate_core_born_and_decays(x, shat, sqrtshat, xjac, &
                                            xpswgt, visible_momenta, pass)
    double precision, intent(in) :: x(99), shat, sqrtshat
    double precision, intent(inout) :: xjac, xpswgt
    double precision, intent(out) :: visible_momenta(0:3, nexternal - 1)
    logical, intent(out) :: pass

    integer :: context, core_count, final_count, leg, index
    double precision :: parent(0:3)
    double precision :: final_masses(nexternal)
    double precision :: final_momenta(0:3, nexternal)
    type(factorized_measure_state) :: production_measure
    logical :: measure_available

    call require_enabled()
    node_rest_storage = 0d0
    node_rest_valid = .false.
    context = born_context()
    core_count = context_core_count(context)
    final_count = core_count - nincoming
    parent = 0d0
    parent(0) = sqrtshat
    do leg = 1, final_count
      final_masses(leg) = core_mass(context, nincoming + leg)
    end do

    production_measure%jacobian = xjac
    production_measure%phase_space_weight = xpswgt
    call generate_nbody(parent, final_count, final_masses, x, 1, &
         final_momenta, production_measure%jacobian, &
         production_measure%phase_space_weight, pass)
    if (.not. pass) return
    if (final_count == 1) production_measure%phase_space_weight = &
         production_measure%phase_space_weight/(2d0*sqrtshat)
    call store_factorized_base_measure(0, production_measure)

    core_born_storage = 0d0
    call fill_incoming_momenta(shat, sqrtshat, context, &
                               core_born_storage, pass)
    if (.not. pass) return
    do leg = 1, final_count
      core_born_storage(:, nincoming + leg) = final_momenta(:, leg)
    end do
    call store_factorized_block_momenta(soft_counterevent, 0, core_count, &
                                        core_born_storage)

    visible_momenta = 0d0
    index = decay_variable_start()
    call expand_context(context, soft_counterevent, core_born_storage, x, &
                        index, visible_momenta, .true., pass)
    if (.not. pass) return
    if (index /= decay_variable_start() + decay_random_dimension()) then
      call fail_kinematics('decay random-variable accounting is inconsistent')
    end if
    call compose_factorized_base_measure( &
         xjac, xpswgt, measure_available)
    if (.not. measure_available) then
      call fail_kinematics('the production measure is unavailable')
    end if
  end subroutine generate_core_born_and_decays


  subroutine expand_real_decay_momenta(configuration, event_slot, x, &
                                       core_momenta, &
                                       visible_momenta, pass)
    integer, intent(in) :: configuration, event_slot
    double precision, intent(in) :: x(99)
    double precision, intent(in) :: core_momenta(0:3, nexternal)
    double precision, intent(out) :: visible_momenta(0:3, nexternal)
    logical, intent(out) :: pass
    integer :: context, index

    call require_enabled()
    context = context_for_fks(configuration)
    visible_momenta = 0d0
    index = decay_variable_start()
    call expand_context(context, event_slot, core_momenta, x, index, &
                        visible_momenta, .false., pass)
    if (.not. pass) return
    if (index /= decay_variable_start() + decay_random_dimension()) then
      call fail_kinematics('real decay random-variable accounting is inconsistent')
    end if
  end subroutine expand_real_decay_momenta


  subroutine expand_context(context, event_slot, core_momenta, x, index, &
                            visible_momenta, keep_weight, pass)
    integer, intent(in) :: context, event_slot
    double precision, intent(in) :: core_momenta(0:, :), x(99)
    integer, intent(inout) :: index
    double precision, intent(out) :: visible_momenta(0:, :)
    logical, intent(in) :: keep_weight
    logical, intent(out) :: pass

    integer :: leg, target
    call set_decay_cut_mask(context)
    pass = .true.
    do leg = 1, context_core_count(context)
      target = core_target_id(context, leg)
      if (core_target_kind(context, leg) == direct_leg_target) then
        visible_momenta(:, target) = core_momenta(:, leg)
      else
        call expand_node(context, event_slot, target, core_momenta(:, leg), &
                         x, index, visible_momenta, keep_weight, pass)
        if (.not. pass) return
      end if
    end do
  end subroutine expand_context


  subroutine set_decay_cut_mask(context)
    integer, intent(in) :: context
    integer :: leaf, visible

    call require_enabled()
    event_from_decay = .false.
    do leaf = 1, decay_leaf_count()
      visible = leaf_visible_leg(context, leaf)
      if (visible > nincoming) event_from_decay(visible) = .true.
    end do
  end subroutine set_decay_cut_mask


  recursive subroutine expand_node(context, event_slot, node, &
                                   parent_momentum, x, &
                                   index, visible_momenta, &
                                   keep_weight, pass)
    integer, intent(in) :: context, event_slot, node
    double precision, intent(in) :: parent_momentum(0:3), x(99)
    integer, intent(inout) :: index
    double precision, intent(inout) :: visible_momenta(0:, :)
    logical, intent(in) :: keep_weight
    logical, intent(out) :: pass

    integer :: child_count, child, identifier, child_kind
    double precision :: child_masses(nexternal)
    double precision :: child_momenta(0:3, nexternal)
    double precision :: rest_momenta(0:3, nexternal)
    double precision :: block_momenta(0:3, nexternal)
    double precision :: local_jacobian, local_weight
    type(factorized_measure_state) :: decay_measure

    child_count = node_child_count(node)
    do child = 1, child_count
      identifier = node_child_id(node, child)
      child_kind = node_child_kind(node, child)
      if (child_kind == decay_node_child) then
        child_masses(child) = node_masses(identifier)
      else
        child_masses(child) = leaf_masses(identifier)
      end if
    end do

    local_jacobian = 1d0
    local_weight = 1d0
    if (keep_weight) then
      call generate_nbody_rest(node_masses(node), child_count, &
           child_masses, x, index, rest_momenta, local_jacobian, &
           local_weight, pass)
      if (.not. pass) return
      node_rest_storage(:, :, node) = 0d0
      node_rest_storage(:, 1:child_count, node) = &
           rest_momenta(:, 1:child_count)
      node_rest_valid(node) = .true.
    else
      if (.not. node_rest_valid(node)) then
        call fail_kinematics( &
             'a decay rest-frame block is unavailable for an event')
      end if
      rest_momenta = node_rest_storage(:, :, node)
    end if
    call boost_nbody_from_rest(rest_momenta, child_count, parent_momentum, &
                               node_masses(node), child_momenta)
    index = index + 3*child_count - 4
    if (keep_weight) then
      decay_measure%jacobian = local_jacobian
      decay_measure%phase_space_weight = &
           local_weight*decay_node_nwa_weight(node)
      call store_factorized_base_measure(node, decay_measure)
    end if

    block_momenta = 0d0
    block_momenta(:, 1) = parent_momentum
    block_momenta(:, 2:child_count + 1) = &
         child_momenta(:, 1:child_count)
    ! Slot zero is the Born matrix-element layout.  Its projected real
    ! counterevent contains an extra zero-emission leg and must not replace
    ! the Born-local block cached above.
    if (keep_weight .or. event_slot /= soft_counterevent) then
      call store_factorized_block_momenta(event_slot, node, &
                                          child_count + 1, block_momenta)
    end if

    do child = 1, child_count
      identifier = node_child_id(node, child)
      child_kind = node_child_kind(node, child)
      if (child_kind == decay_node_child) then
        call expand_node(context, event_slot, identifier, &
                         child_momenta(:, child), x, index, &
                         visible_momenta, keep_weight, pass)
        if (.not. pass) return
      else
        visible_momenta(:, leaf_visible_leg(context, identifier)) = &
             child_momenta(:, child)
      end if
    end do
  end subroutine expand_node


  subroutine fill_incoming_momenta(shat, sqrtshat, context, momenta, pass)
    double precision, intent(in) :: shat, sqrtshat
    integer, intent(in) :: context
    double precision, intent(inout) :: momenta(0:, :)
    logical, intent(out) :: pass
    double precision :: mass1, mass2, lambda_value, momentum_length

    pass = .false.
    if (nincoming == 1) then
      momenta(:, 1) = 0d0
      momenta(0, 1) = sqrtshat
      pass = .true.
      return
    end if
    mass1 = core_mass(context, 1)
    mass2 = core_mass(context, 2)
    lambda_value = phase_space_lambda(shat, mass1**2, mass2**2)
    if (lambda_value < 0d0) return
    momentum_length = sqrt(lambda_value)/(2d0*sqrtshat)
    momenta(:, 1) = 0d0
    momenta(:, 2) = 0d0
    momenta(0, 1) = (shat + mass1**2 - mass2**2)/(2d0*sqrtshat)
    momenta(0, 2) = (shat + mass2**2 - mass1**2)/(2d0*sqrtshat)
    momenta(3, 1) = momentum_length
    momenta(3, 2) = -momentum_length
    pass = .true.
  end subroutine fill_incoming_momenta


  subroutine store_core_event_momenta(event_slot, momenta, particle_count)
    integer, intent(in) :: event_slot, particle_count
    double precision, intent(in) :: momenta(0:, :)
    call require_enabled()
    if (event_slot < 0 .or. event_slot > 3) then
      call fail_kinematics('event slot is out of range')
    end if
    if (particle_count < 1 .or. particle_count > nexternal) then
      call fail_kinematics('the production block size is out of range')
    end if
    if (size(momenta, 1) < 4 .or. size(momenta, 2) < particle_count) then
      call fail_kinematics('the production block momentum shape is invalid')
    end if
    call store_factorized_kernel_momenta( &
         event_slot, 0, particle_count, momenta)
    if (event_slot /= soft_counterevent) then
      call store_factorized_block_momenta( &
           event_slot, 0, particle_count, momenta)
    end if
  end subroutine store_core_event_momenta


  subroutine get_core_born_momenta(momenta)
    double precision, intent(out) :: momenta(0:3, nexternal - 1)
    call require_enabled()
    momenta = core_born_storage
  end subroutine get_core_born_momenta


  subroutine get_core_mass_buffer(context, masses)
    integer, intent(in) :: context
    double precision, intent(out) :: masses(nexternal)
    integer :: leg
    call require_enabled()
    masses = 0d0
    do leg = 1, context_core_count(context)
      masses(leg) = core_mass(context, leg)
    end do
  end subroutine get_core_mass_buffer


  integer function active_core_count(configuration)
    integer, intent(in) :: configuration
    call require_enabled()
    active_core_count = context_core_count(context_for_fks(configuration))
  end function active_core_count


  double precision function fks_leg_mass(configuration, leg)
    integer, intent(in) :: configuration, leg
    call require_enabled()
    fks_leg_mass = core_mass(context_for_fks(configuration), leg)
  end function fks_leg_mass


  subroutine map_core_color_pair(core_first, core_second, visible_first, &
                                 visible_second)
    integer, intent(in) :: core_first, core_second
    integer, intent(out) :: visible_first, visible_second
    call require_enabled()
    call visible_color_pair(core_first, core_second, visible_first, &
                            visible_second)
  end subroutine map_core_color_pair


  subroutine contract_visible_momenta(context, visible_momenta, &
                                      core_momenta)
    integer, intent(in) :: context
    double precision, intent(in) :: visible_momenta(0:3, nexternal)
    double precision, intent(out) :: core_momenta(0:3, nexternal)
    integer :: leg, target

    call require_enabled()
    core_momenta = 0d0
    do leg = 1, context_core_count(context)
      target = core_target_id(context, leg)
      if (core_target_kind(context, leg) == direct_leg_target) then
        core_momenta(:, leg) = visible_momenta(:, target)
      else
        call sum_decay_node_momenta(context, target, visible_momenta, &
                                    core_momenta(:, leg))
      end if
    end do
  end subroutine contract_visible_momenta


  recursive subroutine sum_decay_node_momenta(context, node, &
                                               visible_momenta, momentum)
    integer, intent(in) :: context, node
    double precision, intent(in) :: visible_momenta(0:3, nexternal)
    double precision, intent(out) :: momentum(0:3)
    double precision :: child_momentum(0:3)
    integer :: child, identifier

    momentum = 0d0
    do child = 1, node_child_count(node)
      identifier = node_child_id(node, child)
      if (node_child_kind(node, child) == decay_node_child) then
        call sum_decay_node_momenta(context, identifier, visible_momenta, &
                                    child_momentum)
      else
        child_momentum = &
             visible_momenta(:, leaf_visible_leg(context, identifier))
      end if
      momentum = momentum + child_momentum
    end do
  end subroutine sum_decay_node_momenta


  double precision function decay_node_nwa_weight(node)
    integer, intent(in) :: node
    double precision :: denominator_scale, physical_width
    include 'decay_matrix_factorization.inc'
    call require_enabled()
    if (node < 1 .or. node > decay_node_count()) then
      call fail_kinematics('a decay measure requested an invalid node')
    end if
    physical_width = decay_physical_width(node_pdg(node))
    if (factorized_decay_matrix_elements) then
      decay_node_nwa_weight = &
           1d0/(2d0*node_masses(node)*physical_width)
    else
      denominator_scale = decay_dummy_width_ratio()*node_masses(node)**2
      decay_node_nwa_weight = denominator_scale**2/ &
           (2d0*node_masses(node)*physical_width)
    end if
  end function decay_node_nwa_weight


  subroutine require_enabled()
    if (.not. initialized) call initialize_decay_chain_kinematics()
    if (.not. has_decay_chains()) then
      call fail_kinematics('no decay-chain metadata are present')
    end if
  end subroutine require_enabled


  subroutine fail_kinematics(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in decay_chain_kinematics: '//trim(message)
    stop 1
  end subroutine fail_kinematics

end module decay_chain_kinematics
