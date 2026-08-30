module decay_chain_kinematics
  use process_dimensions, only: nexternal, nincoming, max_branch, &
       validate_process_dimensions
  use fnlo_process_common, only: event_from_decay => from_decay, &
       soft_counterevent, this_config, config_map, config_forest, &
       config_sprop, config_mass, config_width
  use factorized_phase_space, only: factorized_measure_state, &
       store_factorized_block_momenta, store_factorized_kernel_momenta, &
       store_factorized_embedded_momenta, &
       fetch_factorized_embedded_momenta, &
       store_factorized_base_measure, compose_factorized_base_measure
  use factorized_block_kinematics, only: &
       generate_nbody => generate_factorized_nbody, &
       generate_nbody_rest => generate_factorized_nbody_rest, &
       generate_decay_tree_rest => generate_factorized_decay_tree_rest, &
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
  public :: decay_variable_start
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
    call store_factorized_embedded_momenta( &
         soft_counterevent, 0, core_count, core_born_storage)

    index = decay_variable_start()
    call sample_decay_context(context, x, index, pass)
    if (.not. pass) return
    if (index /= decay_variable_start() + decay_random_dimension()) then
      call fail_kinematics('decay random-variable accounting is inconsistent')
    end if
    call embed_decay_context( &
         context, soft_counterevent, core_born_storage, .true., pass)
    if (.not. pass) return
    call materialize_decay_event( &
         context, soft_counterevent, visible_momenta, pass)
    if (.not. pass) return
    call compose_factorized_base_measure( &
         xjac, xpswgt, measure_available)
    if (.not. measure_available) then
      call fail_kinematics('the production measure is unavailable')
    end if
  end subroutine generate_core_born_and_decays


  subroutine expand_real_decay_momenta(configuration, event_slot, &
                                       core_momenta, &
                                       visible_momenta, pass)
    integer, intent(in) :: configuration, event_slot
    double precision, intent(in) :: core_momenta(0:3, nexternal)
    double precision, intent(out) :: visible_momenta(0:3, nexternal)
    logical, intent(out) :: pass
    integer :: context

    call require_enabled()
    context = context_for_fks(configuration)
    call embed_decay_context( &
         context, event_slot, core_momenta, &
         event_slot /= soft_counterevent, pass)
    if (.not. pass) return
    call materialize_decay_event(context, event_slot, visible_momenta, pass)
  end subroutine expand_real_decay_momenta


  subroutine sample_decay_context(context, x, index, pass)
    integer, intent(in) :: context
    double precision, intent(in) :: x(99)
    integer, intent(inout) :: index
    logical, intent(out) :: pass

    integer :: leg, target
    pass = .true.
    do leg = 1, context_core_count(context)
      target = core_target_id(context, leg)
      if (core_target_kind(context, leg) /= decay_node_target) cycle
      call sample_decay_node(target, x, index, pass)
      if (.not. pass) return
    end do
  end subroutine sample_decay_context


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


  recursive subroutine sample_decay_node(node, x, index, pass)
    integer, intent(in) :: node
    double precision, intent(in) :: x(99)
    integer, intent(inout) :: index
    logical, intent(out) :: pass

    integer :: child_count, child, identifier, child_kind
    integer :: local_tree(2, -nexternal:-1)
    integer :: local_propagator_ids(-nexternal:-1)
    double precision :: child_masses(nexternal)
    double precision :: local_propagator_masses(-nexternal:-1)
    double precision :: local_propagator_widths(-nexternal:-1)
    double precision :: rest_momenta(0:3, nexternal)
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
    call build_decay_configuration_tree( &
         node, child_count, local_tree, local_propagator_masses, &
         local_propagator_widths, local_propagator_ids)
    call generate_decay_tree_rest( &
         node_masses(node), child_count, child_masses, local_tree, &
         local_propagator_masses, local_propagator_widths, &
         local_propagator_ids, x, index, rest_momenta, local_jacobian, &
         local_weight, pass)
    if (.not. pass) return
    node_rest_storage(:, :, node) = 0d0
    node_rest_storage(:, 1:child_count, node) = &
         rest_momenta(:, 1:child_count)
    node_rest_valid(node) = .true.
    index = index + 3*child_count - 4
    decay_measure%jacobian = local_jacobian
    decay_measure%phase_space_weight = &
         local_weight*decay_node_nwa_weight(node)
    call store_factorized_base_measure(node, decay_measure)

    do child = 1, child_count
      if (node_child_kind(node, child) /= decay_node_child) cycle
      call sample_decay_node(node_child_id(node, child), x, index, pass)
      if (.not. pass) return
    end do
  end subroutine sample_decay_node


  subroutine build_decay_configuration_tree( &
       node, child_count, local_tree, local_masses, local_widths, &
       local_propagator_ids)
    integer, intent(in) :: node, child_count
    integer, intent(out) :: local_tree(2, -nexternal:-1)
    double precision, intent(out) :: local_masses(-nexternal:-1)
    double precision, intent(out) :: local_widths(-nexternal:-1)
    integer, intent(out) :: local_propagator_ids(-nexternal:-1)

    logical :: atomic_masks(nexternal - 1, nexternal)
    logical :: node_mask(nexternal - 1), branch_mask(nexternal - 1)
    logical :: used_children(nexternal), valid
    integer :: child, branch, root_branch, root_matches
    integer :: internal_count, local_root

    if (this_config < 1 .or. this_config > config_map(0, 0)) then
      call fail_kinematics( &
           'the selected Born configuration is outside the generated range')
    end if
    if (child_count < 2 .or. child_count > nexternal) then
      call fail_kinematics('a decay node has an invalid child count')
    end if

    local_tree = 0
    local_masses = 0d0
    local_widths = 0d0
    local_propagator_ids = 0
    atomic_masks = .false.
    used_children = .false.
    call mark_decay_node_visible_legs(node, node_mask)
    do child = 1, child_count
      call mark_decay_child_visible_legs(node, child, &
                                          atomic_masks(:, child))
      if (.not. any(atomic_masks(:, child))) then
        call fail_kinematics('a decay child has no visible descendants')
      end if
    end do
    if (.not. all(any(atomic_masks(:, 1:child_count), dim=2) .eqv. &
         node_mask)) then
      call fail_kinematics('decay-child visible masks are inconsistent')
    end if

    root_branch = 0
    root_matches = 0
    do branch = -max_branch, -1
      if (all(config_forest(:, branch, this_config, 0) == 0)) cycle
      call configuration_descendant_mask( &
           branch, this_config, branch_mask, valid, 0)
      if (.not. valid) cycle
      if (all(branch_mask .eqv. node_mask)) then
        root_branch = branch
        root_matches = root_matches + 1
      end if
    end do
    if (root_matches /= 1) then
      call fail_kinematics( &
           'the unique LO decay topology is absent from the selected Born '// &
           'configuration')
    end if

    internal_count = 0
    call copy_decay_configuration_subtree( &
         root_branch, .false., child_count, atomic_masks, used_children, &
         internal_count, local_root, local_tree, local_masses, &
         local_widths, local_propagator_ids)
    if (internal_count /= child_count - 1 .or. &
        local_root /= -(child_count - 1) .or. &
        .not. all(used_children(1:child_count))) then
      call fail_kinematics( &
           'the selected Born configuration does not contain one binary '// &
           'topology for the decay node')
    end if
  end subroutine build_decay_configuration_tree


  recursive subroutine copy_decay_configuration_subtree( &
       source, allow_atomic, child_count, atomic_masks, used_children, &
       internal_count, local_label, local_tree, local_masses, &
       local_widths, local_propagator_ids)
    integer, intent(in) :: source, child_count
    logical, intent(in) :: allow_atomic
    logical, intent(in) :: atomic_masks(nexternal - 1, nexternal)
    logical, intent(inout) :: used_children(nexternal)
    integer, intent(inout) :: internal_count
    integer, intent(out) :: local_label
    integer, intent(inout) :: local_tree(2, -nexternal:-1)
    double precision, intent(inout) :: local_masses(-nexternal:-1)
    double precision, intent(inout) :: local_widths(-nexternal:-1)
    integer, intent(inout) :: local_propagator_ids(-nexternal:-1)

    logical :: source_mask(nexternal - 1), valid
    integer :: atomic_child, first_local, second_local

    call configuration_descendant_mask( &
         source, this_config, source_mask, valid, 0)
    if (.not. valid) then
      call fail_kinematics('a decay configuration subtree is malformed')
    end if

    if (allow_atomic) then
      atomic_child = matching_atomic_child( &
           source_mask, child_count, atomic_masks)
      if (atomic_child > 0) then
        if (used_children(atomic_child)) then
          call fail_kinematics( &
               'a decay topology uses one child more than once')
        end if
        used_children(atomic_child) = .true.
        local_label = atomic_child
        return
      end if
    end if
    if (source >= 0 .or. source < -max_branch .or. &
        source < -nexternal) then
      call fail_kinematics( &
           'a decay topology cannot be represented by generated branches')
    end if
    if (any(config_forest(:, source, this_config, 0) == 0)) then
      call fail_kinematics( &
           'factorized decay phase space requires a binary LO topology')
    end if

    call copy_decay_configuration_subtree( &
         config_forest(1, source, this_config, 0), .true., child_count, &
         atomic_masks, used_children, internal_count, first_local, &
         local_tree, local_masses, local_widths, local_propagator_ids)
    call copy_decay_configuration_subtree( &
         config_forest(2, source, this_config, 0), .true., child_count, &
         atomic_masks, used_children, internal_count, second_local, &
         local_tree, local_masses, local_widths, local_propagator_ids)
    internal_count = internal_count + 1
    if (internal_count > child_count - 1) then
      call fail_kinematics('a decay topology contains too many branches')
    end if
    local_label = -internal_count
    local_tree(1, local_label) = first_local
    local_tree(2, local_label) = second_local
    local_masses(local_label) = &
         config_mass(source, this_config, 0)
    local_widths(local_label) = &
         config_width(source, this_config, 0)
    local_propagator_ids(local_label) = &
         config_sprop(source, this_config, 0)
  end subroutine copy_decay_configuration_subtree


  integer function matching_atomic_child(mask, child_count, atomic_masks)
    logical, intent(in) :: mask(nexternal - 1)
    integer, intent(in) :: child_count
    logical, intent(in) :: atomic_masks(nexternal - 1, nexternal)
    integer :: child, matches

    matching_atomic_child = 0
    matches = 0
    do child = 1, child_count
      if (.not. all(mask .eqv. atomic_masks(:, child))) cycle
      matching_atomic_child = child
      matches = matches + 1
    end do
    if (matches > 1) then
      call fail_kinematics( &
           'two immediate decay children have identical visible content')
    end if
  end function matching_atomic_child


  recursive subroutine configuration_descendant_mask( &
       source, configuration, mask, valid, depth)
    integer, intent(in) :: source, configuration, depth
    logical, intent(out) :: mask(nexternal - 1), valid
    logical :: child_mask(nexternal - 1), child_valid
    integer :: child

    mask = .false.
    valid = .false.
    if (depth > max_branch) return
    if (source > 0) then
      if (source > nexternal - 1) return
      mask(source) = .true.
      valid = .true.
      return
    end if
    if (source >= 0 .or. source < -max_branch) return
    if (all(config_forest(:, source, configuration, 0) == 0)) return
    valid = .true.
    do child = 1, 2
      call configuration_descendant_mask( &
           config_forest(child, source, configuration, 0), &
           configuration, child_mask, child_valid, depth + 1)
      if (.not. child_valid .or. any(mask .and. child_mask)) then
        valid = .false.
        mask = .false.
        return
      end if
      mask = mask .or. child_mask
    end do
  end subroutine configuration_descendant_mask


  recursive subroutine mark_decay_child_visible_legs(node, child, mask)
    integer, intent(in) :: node, child
    logical, intent(out) :: mask(nexternal - 1)
    integer :: identifier, visible

    mask = .false.
    identifier = node_child_id(node, child)
    if (node_child_kind(node, child) == decay_node_child) then
      call mark_decay_node_visible_legs(identifier, mask)
    else
      visible = leaf_visible_leg(born_context(), identifier)
      if (visible < 1 .or. visible > nexternal - 1) then
        call fail_kinematics('a decay leaf has an invalid Born target')
      end if
      mask(visible) = .true.
    end if
  end subroutine mark_decay_child_visible_legs


  recursive subroutine mark_decay_node_visible_legs(node, mask)
    integer, intent(in) :: node
    logical, intent(out) :: mask(nexternal - 1)
    logical :: child_mask(nexternal - 1)
    integer :: child

    mask = .false.
    do child = 1, node_child_count(node)
      call mark_decay_child_visible_legs(node, child, child_mask)
      if (any(mask .and. child_mask)) then
        call fail_kinematics('decay nodes share a visible Born target')
      end if
      mask = mask .or. child_mask
    end do
  end subroutine mark_decay_node_visible_legs


  subroutine embed_decay_context(context, event_slot, core_momenta, &
                                 store_matrix_blocks, pass)
    integer, intent(in) :: context, event_slot
    double precision, intent(in) :: core_momenta(0:, :)
    logical, intent(in) :: store_matrix_blocks
    logical, intent(out) :: pass
    integer :: leg, target

    call store_factorized_embedded_momenta( &
         event_slot, 0, context_core_count(context), core_momenta)
    pass = .true.
    do leg = 1, context_core_count(context)
      if (core_target_kind(context, leg) /= decay_node_target) cycle
      target = core_target_id(context, leg)
      call embed_decay_node(event_slot, target, core_momenta(:, leg), &
                            store_matrix_blocks, pass)
      if (.not. pass) return
    end do
  end subroutine embed_decay_context


  recursive subroutine embed_decay_node(event_slot, node, parent_momentum, &
                                        store_matrix_blocks, pass)
    integer, intent(in) :: event_slot, node
    double precision, intent(in) :: parent_momentum(0:3)
    logical, intent(in) :: store_matrix_blocks
    logical, intent(out) :: pass
    integer :: child_count, child, identifier
    double precision :: child_momenta(0:3, nexternal)
    double precision :: block_momenta(0:3, nexternal)

    child_count = node_child_count(node)
    if (.not. node_rest_valid(node)) then
      call fail_kinematics( &
           'a decay rest-frame block is unavailable for embedding')
    end if
    call boost_nbody_from_rest( &
         node_rest_storage(:, :, node), child_count, parent_momentum, &
         node_masses(node), child_momenta)

    block_momenta = 0d0
    block_momenta(:, 1) = parent_momentum
    block_momenta(:, 2:child_count + 1) = &
         child_momenta(:, 1:child_count)
    call store_factorized_embedded_momenta( &
         event_slot, node, child_count + 1, block_momenta)
    if (store_matrix_blocks) then
      call store_factorized_block_momenta(event_slot, node, &
                                          child_count + 1, block_momenta)
    end if

    pass = .true.
    do child = 1, child_count
      identifier = node_child_id(node, child)
      if (node_child_kind(node, child) /= decay_node_child) cycle
      call embed_decay_node(event_slot, identifier, &
                            child_momenta(:, child), &
                            store_matrix_blocks, pass)
      if (.not. pass) return
    end do
  end subroutine embed_decay_node


  subroutine materialize_decay_event(context, event_slot, visible_momenta, &
                                     pass)
    integer, intent(in) :: context, event_slot
    double precision, intent(out) :: visible_momenta(0:, :)
    logical, intent(out) :: pass
    double precision :: core_momenta(0:3, nexternal)
    integer :: leg, target, core_count
    logical :: available

    call set_decay_cut_mask(context)
    visible_momenta = 0d0
    core_count = context_core_count(context)
    call fetch_factorized_embedded_momenta( &
         event_slot, 0, core_count, core_momenta(:, 1:core_count), available)
    if (.not. available) then
      call fail_kinematics('the embedded production block is unavailable')
    end if
    pass = .true.
    do leg = 1, core_count
      target = core_target_id(context, leg)
      if (core_target_kind(context, leg) == direct_leg_target) then
        visible_momenta(:, target) = core_momenta(:, leg)
      else
        call materialize_decay_node( &
             context, event_slot, target, visible_momenta, pass)
        if (.not. pass) return
      end if
    end do
  end subroutine materialize_decay_event


  recursive subroutine materialize_decay_node( &
       context, event_slot, node, visible_momenta, pass)
    integer, intent(in) :: context, event_slot, node
    double precision, intent(inout) :: visible_momenta(0:, :)
    logical, intent(out) :: pass
    double precision :: block_momenta(0:3, nexternal)
    integer :: child_count, child, identifier
    logical :: available

    child_count = node_child_count(node)
    call fetch_factorized_embedded_momenta( &
         event_slot, node, child_count + 1, &
         block_momenta(:, 1:child_count + 1), available)
    if (.not. available) then
      call fail_kinematics('an embedded decay block is unavailable')
    end if
    pass = .true.
    do child = 1, child_count
      identifier = node_child_id(node, child)
      if (node_child_kind(node, child) == decay_node_child) then
        call materialize_decay_node( &
             context, event_slot, identifier, visible_momenta, pass)
        if (.not. pass) return
      else
        visible_momenta(:, leaf_visible_leg(context, identifier)) = &
             block_momenta(:, child + 1)
      end if
    end do
  end subroutine materialize_decay_node


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
    call store_factorized_embedded_momenta( &
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
