module decay_chain_metadata
  use process_dimensions, only: nexternal, nincoming, fks_configs, &
                                validate_process_dimensions
  use fks_metadata, only: validate_fks_metadata, fks_i_d, fks_j_d
  use nlo_decay_metadata, only: has_nlo_decay, nlo_decay_node_count
  use nlo_contribution_bundle, only: has_nlo_contribution_bundle, &
       active_contribution_is_production, contribution_fks_first, &
       contribution_fks_last, global_fks_configuration, &
       local_fks_configuration
  implicit none
  private

  integer, parameter, public :: direct_leg_target = 1
  integer, parameter, public :: decay_node_target = 2
  integer, parameter, public :: decay_leaf_child = 1
  integer, parameter, public :: decay_node_child = 2

  logical, save :: initialized = .false.
  logical, save :: enabled = .false.
  integer, save :: metadata_format = 0
  integer, save :: number_of_nodes = 0
  integer, save :: number_of_leaves = 0
  integer, save :: number_of_contexts = 0
  integer, save :: number_of_fks_maps = 0
  integer, save :: number_of_color_records = 0
  integer, save :: number_of_generated_color_links = 0

  integer, allocatable, save :: node_parent_values(:)
  integer, allocatable, save :: node_pdg_values(:)
  integer, allocatable, save :: node_qcd_orders(:)
  integer, allocatable, save :: node_carrier_values(:)
  integer, allocatable, save :: node_child_counts(:)
  integer, allocatable, save :: node_child_kinds(:, :)
  integer, allocatable, save :: node_child_ids(:, :)
  integer, allocatable, save :: leaf_parent_values(:)
  integer, allocatable, save :: leaf_pdg_values(:)

  character(len=16), allocatable, save :: context_kind_values(:)
  integer, allocatable, save :: context_source_values(:)
  integer, allocatable, save :: context_core_counts(:)
  integer, allocatable, save :: context_visible_counts(:)
  integer, allocatable, save :: core_target_kinds(:, :)
  integer, allocatable, save :: core_target_ids(:, :)
  integer, allocatable, save :: core_pdg_values(:, :)
  logical, allocatable, save :: core_final_values(:, :)
  integer, allocatable, save :: leaf_visible_values(:, :)

  integer, allocatable, save :: fks_context_values(:)
  integer, allocatable, save :: fks_i_values(:)
  integer, allocatable, save :: fks_j_values(:)
  integer, allocatable, save :: fks_ij_values(:)
  integer, allocatable, save :: color_core_first(:)
  integer, allocatable, save :: color_core_second(:)
  integer, allocatable, save :: color_visible_first(:)
  integer, allocatable, save :: color_visible_second(:)
  integer, allocatable, save :: color_generated_index(:)

  public :: initialize_decay_chain_metadata
  public :: has_decay_chains, decay_node_count, decay_leaf_count
  public :: decay_random_dimension, real_phase_space_dimension
  public :: born_context, context_for_fks, context_core_count
  public :: core_target_kind, core_target_id
  public :: core_leg_pdg, leaf_visible_leg
  public :: node_pdg, node_qcd_order, node_child_count
  public :: node_child_kind, node_child_id, leaf_pdg
  public :: visible_color_pair

contains

  subroutine initialize_decay_chain_metadata()
    implicit none
    logical :: exists
    integer :: unit_number, ios, node, leaf, context, core_leg
    integer :: configuration, child, color_capacity
    integer :: parent, pdg, qcd_order, carrier, child_count, source
    integer :: core_count, visible_count, target, visible
    integer :: core_first, core_second, visible_first, visible_second
    integer :: generated_index
    character(len=512) :: line
    character(len=32) :: keyword, kind, state
    character(len=16) :: child_words(nexternal)
    integer :: child_ids(nexternal)
    logical :: end_seen

    if (initialized) return
    call validate_process_dimensions()
    inquire(file='decay_chain_info.dat', exist=exists)
    if (.not. exists) then
      initialized = .true.
      enabled = .false.
      return
    end if

    open(newunit=unit_number, file='decay_chain_info.dat', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) call fail_metadata('cannot open decay_chain_info.dat')

    do
      read(unit_number, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) call fail_metadata('cannot read decay_chain_info.dat')
      if (len_trim(line) == 0) cycle
      read(line, *, iostat=ios) keyword
      if (ios /= 0) call fail_metadata('malformed metadata record')
      select case (trim(keyword))
      case ('FORMAT')
        read(line, *, iostat=ios) keyword, metadata_format
      case ('COUNTS')
        read(line, *, iostat=ios) keyword, number_of_nodes, &
             number_of_leaves, number_of_contexts, number_of_fks_maps, &
             number_of_generated_color_links
      end select
      if (ios /= 0) call fail_metadata('malformed metadata header')
    end do

    if (metadata_format /= 4) then
      call fail_metadata('FORMAT 4 is required; regenerate the process')
    end if
    if (number_of_nodes < 1 .or. number_of_leaves < 1 .or. &
        number_of_contexts < 1 .or. number_of_fks_maps < 1 .or. &
        number_of_generated_color_links < 0) then
      call fail_metadata('invalid metadata counts')
    end if
    if (has_nlo_contribution_bundle()) then
      if (number_of_fks_maps /= contribution_fks_last(1) - &
          contribution_fks_first(1) + 1) then
        call fail_metadata('production contribution FKS count disagrees')
      end if
    else if (number_of_fks_maps /= fks_configs) then
      call fail_metadata('FKS map count disagrees with generated data')
    end if
    allocate(node_parent_values(number_of_nodes))
    allocate(node_pdg_values(number_of_nodes))
    allocate(node_qcd_orders(number_of_nodes))
    allocate(node_carrier_values(number_of_nodes))
    allocate(node_child_counts(number_of_nodes))
    allocate(node_child_kinds(nexternal, number_of_nodes))
    allocate(node_child_ids(nexternal, number_of_nodes))
    allocate(leaf_parent_values(number_of_leaves))
    allocate(leaf_pdg_values(number_of_leaves))
    allocate(context_kind_values(number_of_contexts))
    allocate(context_source_values(number_of_contexts))
    allocate(context_core_counts(number_of_contexts))
    allocate(context_visible_counts(number_of_contexts))
    allocate(core_target_kinds(nexternal, number_of_contexts))
    allocate(core_target_ids(nexternal, number_of_contexts))
    allocate(core_pdg_values(nexternal, number_of_contexts))
    allocate(core_final_values(nexternal, number_of_contexts))
    allocate(leaf_visible_values(number_of_leaves, number_of_contexts))
    allocate(fks_context_values(number_of_fks_maps))
    allocate(fks_i_values(number_of_fks_maps))
    allocate(fks_j_values(number_of_fks_maps))
    allocate(fks_ij_values(number_of_fks_maps))
    color_capacity = max(1, nexternal*nexternal)
    allocate(color_core_first(color_capacity))
    allocate(color_core_second(color_capacity))
    allocate(color_visible_first(color_capacity))
    allocate(color_visible_second(color_capacity))
    allocate(color_generated_index(color_capacity))

    node_parent_values = -1
    node_pdg_values = 0
    node_qcd_orders = -1
    node_carrier_values = 0
    node_child_counts = 0
    node_child_kinds = 0
    node_child_ids = 0
    leaf_parent_values = 0
    leaf_pdg_values = 0
    context_kind_values = ''
    context_source_values = 0
    context_core_counts = 0
    context_visible_counts = 0
    core_target_kinds = 0
    core_target_ids = 0
    core_pdg_values = 0
    core_final_values = .false.
    leaf_visible_values = 0
    fks_context_values = 0
    fks_i_values = 0
    fks_j_values = 0
    fks_ij_values = 0
    number_of_color_records = 0
    color_core_first = 0
    color_core_second = 0
    color_visible_first = 0
    color_visible_second = 0
    color_generated_index = 0

    rewind(unit_number)
    end_seen = .false.
    do
      read(unit_number, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) call fail_metadata('cannot read decay metadata body')
      if (len_trim(line) == 0) cycle
      if (end_seen) call fail_metadata('record found after END')
      read(line, *, iostat=ios) keyword
      if (ios /= 0) call fail_metadata('malformed metadata keyword')
      select case (trim(keyword))
      case ('NODE')
        child_words = ''
        child_ids = 0
        read(line, *, iostat=ios) keyword, node, parent, pdg, qcd_order, &
             carrier, child_count
        if (ios /= 0) call fail_metadata('malformed NODE record')
        call check_node(node)
        if (node_parent_values(node) /= -1) then
          call fail_metadata('duplicate NODE record')
        end if
        if (child_count < 2 .or. child_count > nexternal) then
          call fail_metadata('invalid decay-node child count')
        end if
        read(line, *, iostat=ios) keyword, node, parent, pdg, qcd_order, &
             carrier, child_count, (child_words(child), child_ids(child), &
             child=1, child_count)
        node_parent_values(node) = parent
        node_pdg_values(node) = pdg
        node_qcd_orders(node) = qcd_order
        node_carrier_values(node) = carrier
        node_child_counts(node) = child_count
        do child = 1, child_count
          select case (trim(child_words(child)))
          case ('LEAF')
            node_child_kinds(child, node) = decay_leaf_child
          case ('NODE')
            node_child_kinds(child, node) = decay_node_child
          case default
            call fail_metadata('invalid decay-node child kind')
          end select
          node_child_ids(child, node) = child_ids(child)
        end do
      case ('DECAY_LEAF')
        read(line, *, iostat=ios) keyword, leaf, parent, pdg
        call check_leaf(leaf)
        if (leaf_parent_values(leaf) /= 0) then
          call fail_metadata('duplicate DECAY_LEAF record')
        end if
        leaf_parent_values(leaf) = parent
        leaf_pdg_values(leaf) = pdg
      case ('CONTEXT')
        read(line, *, iostat=ios) keyword, context, kind, source, &
             core_count, visible_count
        call check_context(context)
        if (len_trim(context_kind_values(context)) /= 0) then
          call fail_metadata('duplicate CONTEXT record')
        end if
        context_kind_values(context) = kind
        context_source_values(context) = source
        context_core_counts(context) = core_count
        context_visible_counts(context) = visible_count
      case ('CORE_LEG')
        read(line, *, iostat=ios) keyword, context, core_leg, pdg, state
        call check_context(context)
        call check_core_leg_index(core_leg)
        if (core_pdg_values(core_leg, context) /= 0) then
          call fail_metadata('duplicate CORE_LEG record')
        end if
        core_pdg_values(core_leg, context) = pdg
        select case (trim(state))
        case ('I')
          core_final_values(core_leg, context) = .false.
        case ('F')
          core_final_values(core_leg, context) = .true.
        case default
          call fail_metadata('invalid CORE_LEG state')
        end select
      case ('CORE_MAP')
        read(line, *, iostat=ios) keyword, context, core_leg, kind, target
        call check_context(context)
        call check_core_leg_index(core_leg)
        if (core_target_kinds(core_leg, context) /= 0) then
          call fail_metadata('duplicate CORE_MAP record')
        end if
        select case (trim(kind))
        case ('LEG')
          core_target_kinds(core_leg, context) = direct_leg_target
        case ('NODE')
          core_target_kinds(core_leg, context) = decay_node_target
        case default
          call fail_metadata('invalid CORE_MAP target kind')
        end select
        core_target_ids(core_leg, context) = target
      case ('LEAF_MAP')
        read(line, *, iostat=ios) keyword, context, leaf, visible
        call check_context(context)
        call check_leaf(leaf)
        if (leaf_visible_values(leaf, context) /= 0) then
          call fail_metadata('duplicate LEAF_MAP record')
        end if
        leaf_visible_values(leaf, context) = visible
      case ('FKS_MAP')
        read(line, *, iostat=ios) keyword, configuration, context, &
             core_leg, target, parent
        call check_configuration(configuration)
        if (fks_context_values(configuration) /= 0) then
          call fail_metadata('duplicate FKS_MAP record')
        end if
        fks_context_values(configuration) = context
        fks_i_values(configuration) = core_leg
        fks_j_values(configuration) = target
        fks_ij_values(configuration) = parent
      case ('COLOR_LINK')
        read(line, *, iostat=ios) keyword, core_first, core_second, &
             visible_first, visible_second, generated_index
        number_of_color_records = number_of_color_records + 1
        if (number_of_color_records > size(color_core_first)) then
          call fail_metadata('too many COLOR_LINK records')
        end if
        color_core_first(number_of_color_records) = core_first
        color_core_second(number_of_color_records) = core_second
        color_visible_first(number_of_color_records) = visible_first
        color_visible_second(number_of_color_records) = visible_second
        color_generated_index(number_of_color_records) = generated_index
      case ('FORMAT', 'FORCED_SPECIES', 'COUNTS')
        continue
      case ('END')
        end_seen = .true.
      case default
        call fail_metadata('unknown metadata keyword '//trim(keyword))
      end select
      if (ios /= 0) call fail_metadata('malformed metadata record')
    end do
    close(unit_number)

    if (.not. end_seen) call fail_metadata('END record is absent')

    call validate_metadata_values()
    enabled = .true.
    initialized = .true.
  end subroutine initialize_decay_chain_metadata


  logical function has_decay_chains()
    call require_initialized()
    has_decay_chains = enabled
    if (has_nlo_contribution_bundle()) then
      has_decay_chains = has_decay_chains .and. &
           active_contribution_is_production()
    end if
  end function has_decay_chains


  integer function decay_node_count()
    call require_enabled()
    decay_node_count = number_of_nodes
  end function decay_node_count


  integer function decay_leaf_count()
    call require_enabled()
    decay_leaf_count = number_of_leaves
  end function decay_leaf_count


  integer function decay_random_dimension()
    integer :: node
    call require_enabled()
    decay_random_dimension = 0
    do node = 1, number_of_nodes
      decay_random_dimension = decay_random_dimension + &
           3*node_child_counts(node) - 4
    end do
  end function decay_random_dimension


  integer function real_phase_space_dimension()
    call require_initialized()
    real_phase_space_dimension = 3*(nexternal - nincoming) - 4
    if (has_decay_chains()) real_phase_space_dimension = &
         real_phase_space_dimension - number_of_nodes
    ! Every forced node in an NLO-decay forest removes one invariant-mass
    ! integration exactly as for an ordinary decay-chain node.
    if (has_nlo_decay()) real_phase_space_dimension = &
         real_phase_space_dimension - nlo_decay_node_count()
  end function real_phase_space_dimension


  integer function born_context()
    integer :: context
    call require_enabled()
    born_context = 0
    do context = 1, number_of_contexts
      if (trim(context_kind_values(context)) == 'BORN') then
        born_context = context
        exit
      end if
    end do
    if (born_context == 0) call fail_metadata('BORN context is absent')
  end function born_context


  integer function context_for_fks(configuration)
    integer, intent(in) :: configuration
    integer :: local_configuration
    call require_enabled()
    local_configuration = local_fks_configuration(configuration)
    call check_configuration(local_configuration)
    context_for_fks = fks_context_values(local_configuration)
  end function context_for_fks


  integer function context_core_count(context)
    integer, intent(in) :: context
    call require_enabled()
    call check_context(context)
    context_core_count = context_core_counts(context)
  end function context_core_count


  integer function core_target_kind(context, core_leg)
    integer, intent(in) :: context, core_leg
    call validate_core_lookup(context, core_leg)
    core_target_kind = core_target_kinds(core_leg, context)
  end function core_target_kind


  integer function core_target_id(context, core_leg)
    integer, intent(in) :: context, core_leg
    call validate_core_lookup(context, core_leg)
    core_target_id = core_target_ids(core_leg, context)
  end function core_target_id


  integer function core_leg_pdg(context, core_leg)
    integer, intent(in) :: context, core_leg
    call validate_core_lookup(context, core_leg)
    core_leg_pdg = core_pdg_values(core_leg, context)
  end function core_leg_pdg


  integer function leaf_visible_leg(context, leaf)
    integer, intent(in) :: context, leaf
    call require_enabled()
    call check_context(context)
    call check_leaf(leaf)
    leaf_visible_leg = leaf_visible_values(leaf, context)
  end function leaf_visible_leg


  integer function node_pdg(node)
    integer, intent(in) :: node
    call require_enabled()
    call check_node(node)
    node_pdg = node_pdg_values(node)
  end function node_pdg


  integer function node_qcd_order(node)
    integer, intent(in) :: node
    call require_enabled()
    call check_node(node)
    node_qcd_order = node_qcd_orders(node)
  end function node_qcd_order


  integer function node_child_count(node)
    integer, intent(in) :: node
    call require_enabled()
    call check_node(node)
    node_child_count = node_child_counts(node)
  end function node_child_count


  integer function node_child_kind(node, child)
    integer, intent(in) :: node, child
    call validate_child_lookup(node, child)
    node_child_kind = node_child_kinds(child, node)
  end function node_child_kind


  integer function node_child_id(node, child)
    integer, intent(in) :: node, child
    call validate_child_lookup(node, child)
    node_child_id = node_child_ids(child, node)
  end function node_child_id


  integer function leaf_pdg(leaf)
    integer, intent(in) :: leaf
    call require_enabled()
    call check_leaf(leaf)
    leaf_pdg = leaf_pdg_values(leaf)
  end function leaf_pdg


  subroutine visible_color_pair(core_first, core_second, visible_first, &
                                visible_second, generated_index)
    integer, intent(in) :: core_first, core_second
    integer, intent(out) :: visible_first, visible_second
    integer, intent(out), optional :: generated_index
    integer :: record, first, second
    call require_enabled()
    first = min(core_first, core_second)
    second = max(core_first, core_second)
    do record = 1, number_of_color_records
      if (color_core_first(record) == first .and. &
          color_core_second(record) == second) then
        visible_first = color_visible_first(record)
        visible_second = color_visible_second(record)
        if (present(generated_index)) &
             generated_index = color_generated_index(record)
        return
      end if
    end do
    call fail_metadata('requested production color link is absent')
  end subroutine visible_color_pair


  subroutine validate_metadata_values()
    integer :: context, core_leg, configuration, node, child, identifier
    integer :: leaf, record, previous, born_context_id, born_count
    integer :: target_kind, target_id, visible, generated_index
    integer :: generated_first(max(1, number_of_generated_color_links))
    integer :: generated_second(max(1, number_of_generated_color_links))
    logical :: node_referenced(number_of_nodes)
    logical :: leaf_referenced(number_of_leaves)
    logical :: root_used(number_of_nodes)
    logical :: visible_used(nexternal)
    logical :: generated_seen(max(1, number_of_generated_color_links))

    call validate_fks_metadata()

    node_referenced = .false.
    leaf_referenced = .false.
    do node = 1, number_of_nodes
      if (node_parent_values(node) < 0 .or. node_pdg_values(node) == 0 .or. &
          node_qcd_orders(node) < 0) then
        call fail_metadata('incomplete NODE metadata')
      end if
      if (node_parent_values(node) >= node) then
        call fail_metadata('decay nodes are not in parent-before-child order')
      end if
      if (node_carrier_values(node) < 0 .or. &
          node_carrier_values(node) > number_of_leaves) then
        call fail_metadata('decay-node carrier leaf is out of range')
      end if
      do child = 1, node_child_counts(node)
        identifier = node_child_ids(child, node)
        select case (node_child_kinds(child, node))
        case (decay_node_child)
          call check_node(identifier)
          if (identifier <= node) then
            call fail_metadata('nested decay node is not after its parent')
          end if
          if (node_parent_values(identifier) /= node) then
            call fail_metadata('nested decay node has an inconsistent parent')
          end if
          if (node_referenced(identifier)) then
            call fail_metadata('nested decay node is referenced more than once')
          end if
          node_referenced(identifier) = .true.
        case (decay_leaf_child)
          call check_leaf(identifier)
          if (leaf_parent_values(identifier) /= node) then
            call fail_metadata('decay leaf has an inconsistent parent')
          end if
          if (leaf_referenced(identifier)) then
            call fail_metadata('decay leaf is referenced more than once')
          end if
          leaf_referenced(identifier) = .true.
        case default
          call fail_metadata('invalid decay-node child kind')
        end select
      end do
      if (node_carrier_values(node) /= 0) then
        if (.not. carrier_is_descendant(node, node_carrier_values(node))) then
          call fail_metadata('carrier leaf is not a descendant of its node')
        end if
      end if
    end do
    do node = 1, number_of_nodes
      if ((node_parent_values(node) == 0 .and. node_referenced(node)) .or. &
          (node_parent_values(node) /= 0 .and. &
           .not. node_referenced(node))) then
        call fail_metadata('decay-node parentage is incomplete')
      end if
    end do
    do leaf = 1, number_of_leaves
      if (leaf_parent_values(leaf) < 1 .or. &
          leaf_parent_values(leaf) > number_of_nodes .or. &
          leaf_pdg_values(leaf) == 0 .or. .not. leaf_referenced(leaf)) then
        call fail_metadata('incomplete DECAY_LEAF metadata')
      end if
    end do

    born_context_id = 0
    born_count = 0
    do context = 1, number_of_contexts
      if (context_core_counts(context) < nincoming + 1 .or. &
          context_core_counts(context) > nexternal .or. &
          context_visible_counts(context) < context_core_counts(context) .or. &
          context_visible_counts(context) > nexternal) then
        call fail_metadata('invalid context dimensions')
      end if
      if (context_source_values(context) < 1) then
        call fail_metadata('invalid context source index')
      end if
      select case (trim(context_kind_values(context)))
      case ('BORN')
        born_count = born_count + 1
        born_context_id = context
        if (context_visible_counts(context) /= nexternal - 1) then
          call fail_metadata('BORN visible count differs from NEXTERNAL-1')
        end if
      case ('REAL')
        if (context_visible_counts(context) /= nexternal) then
          call fail_metadata('REAL visible count differs from NEXTERNAL')
        end if
      case ('COUNTERTERM')
        continue
      case default
        call fail_metadata('invalid context kind')
      end select

      visible_used = .false.
      root_used = .false.
      do core_leg = 1, context_core_counts(context)
        if (core_pdg_values(core_leg, context) == 0 .or. &
            core_target_kinds(core_leg, context) == 0 .or. &
            core_target_ids(core_leg, context) == 0) then
          call fail_metadata('incomplete core-leg metadata')
        end if
        if ((core_leg <= nincoming) .eqv. &
            core_final_values(core_leg, context)) then
          call fail_metadata('CORE_LEG state disagrees with NINCOMING')
        end if
        target_kind = core_target_kinds(core_leg, context)
        target_id = core_target_ids(core_leg, context)
        select case (target_kind)
        case (direct_leg_target)
          if (target_id < 1 .or. &
              target_id > context_visible_counts(context)) then
            call fail_metadata('direct CORE_MAP target is out of range')
          end if
          if (visible_used(target_id)) then
            call fail_metadata('visible leg is mapped more than once')
          end if
          visible_used(target_id) = .true.
        case (decay_node_target)
          call check_node(target_id)
          if (.not. core_final_values(core_leg, context) .or. &
              node_parent_values(target_id) /= 0 .or. &
              core_pdg_values(core_leg, context) /= &
              node_pdg_values(target_id)) then
            call fail_metadata('invalid decay-node CORE_MAP target')
          end if
          if (root_used(target_id)) then
            call fail_metadata('root decay node is mapped more than once')
          end if
          root_used(target_id) = .true.
        case default
          call fail_metadata('invalid CORE_MAP target kind')
        end select
      end do
      if (context_core_counts(context) < nexternal) then
        if (any(core_pdg_values(context_core_counts(context) + 1:, &
                                context) /= 0) .or. &
            any(core_target_kinds(context_core_counts(context) + 1:, &
                                  context) /= 0)) then
          call fail_metadata('core record lies beyond its context extent')
        end if
      end if
      do node = 1, number_of_nodes
        if (node_parent_values(node) == 0 .and. .not. root_used(node)) then
          call fail_metadata('root decay node is absent from a context')
        end if
      end do
      do leaf = 1, number_of_leaves
        visible = leaf_visible_values(leaf, context)
        if (visible < 1 .or. visible > context_visible_counts(context)) then
          call fail_metadata('LEAF_MAP target is out of range')
        end if
        if (visible_used(visible)) then
          call fail_metadata('visible leg is mapped more than once')
        end if
        visible_used(visible) = .true.
      end do
      if (.not. all(visible_used(1:context_visible_counts(context)))) then
        call fail_metadata('context does not cover every visible leg')
      end if
    end do
    if (born_count /= 1) then
      call fail_metadata('exactly one BORN context is required')
    end if

    do configuration = 1, number_of_fks_maps
      context = fks_context_values(configuration)
      call check_context(context)
      if (trim(context_kind_values(context)) /= 'REAL') then
        call fail_metadata('an FKS map does not reference a REAL context')
      end if
      if (fks_i_values(configuration) < 1 .or. &
          fks_i_values(configuration) > context_core_counts(context) .or. &
          fks_j_values(configuration) < 1 .or. &
          fks_j_values(configuration) > context_core_counts(context) .or. &
          fks_i_values(configuration) == fks_j_values(configuration)) then
        call fail_metadata('an FKS map contains an invalid core index')
      end if
      if (fks_ij_values(configuration) < 1 .or. &
          fks_ij_values(configuration) > &
          context_core_counts(born_context_id)) then
        call fail_metadata('an FKS map contains an invalid Born index')
      end if
      if (fks_i_values(configuration) /= fks_i_d( &
              global_fks_configuration(1, configuration)) .or. &
          fks_j_values(configuration) /= fks_j_d( &
              global_fks_configuration(1, configuration))) then
        call fail_metadata('decay and generated FKS maps disagree')
      end if
    end do

    generated_seen = .false.
    generated_first = 0
    generated_second = 0
    do record = 1, number_of_color_records
      if (color_core_first(record) < 1 .or. &
          color_core_second(record) > &
          context_core_counts(born_context_id) .or. &
          color_core_first(record) > color_core_second(record)) then
        call fail_metadata('production color-link pair is invalid')
      end if
      if (color_visible_first(record) < 1 .or. &
          color_visible_second(record) > &
          context_visible_counts(born_context_id) .or. &
          color_visible_first(record) > color_visible_second(record)) then
        call fail_metadata('visible color-link pair is invalid')
      end if
      generated_index = color_generated_index(record)
      if (generated_index < 1 .or. &
          generated_index > number_of_generated_color_links) then
        call fail_metadata('generated color-link index is invalid')
      end if
      do previous = 1, record - 1
        if (color_core_first(previous) == color_core_first(record) .and. &
            color_core_second(previous) == color_core_second(record)) then
          call fail_metadata('duplicate production color-link pair')
        end if
      end do
      if (generated_seen(generated_index)) then
        if (generated_first(generated_index) /= &
            color_visible_first(record) .or. &
            generated_second(generated_index) /= &
            color_visible_second(record)) then
          call fail_metadata('generated color-link mapping is inconsistent')
        end if
      else
        generated_seen(generated_index) = .true.
        generated_first(generated_index) = color_visible_first(record)
        generated_second(generated_index) = color_visible_second(record)
      end if
    end do
    if (number_of_generated_color_links > 0) then
      if (.not. all(generated_seen(1:number_of_generated_color_links))) then
        call fail_metadata('generated color-link indices are incomplete')
      end if
      do generated_index = 1, number_of_generated_color_links
        do previous = 1, generated_index - 1
          if (generated_first(previous) == generated_first(generated_index) &
              .and. generated_second(previous) == &
              generated_second(generated_index)) then
            call fail_metadata('visible color-link pair has two indices')
          end if
        end do
      end do
    else if (number_of_color_records /= 0) then
      call fail_metadata('color-link records exist without generated links')
    end if
  end subroutine validate_metadata_values


  logical function carrier_is_descendant(node, leaf)
    integer, intent(in) :: node, leaf
    integer :: parent

    carrier_is_descendant = .false.
    parent = leaf_parent_values(leaf)
    do while (parent > 0)
      if (parent > number_of_nodes) return
      if (parent == node) then
        carrier_is_descendant = .true.
        return
      end if
      parent = node_parent_values(parent)
    end do
  end function carrier_is_descendant


  subroutine validate_core_lookup(context, core_leg)
    integer, intent(in) :: context, core_leg
    call require_enabled()
    call check_context(context)
    if (core_leg < 1 .or. core_leg > context_core_counts(context)) then
      call fail_metadata('core leg is outside its context')
    end if
  end subroutine validate_core_lookup


  subroutine validate_child_lookup(node, child)
    integer, intent(in) :: node, child
    call require_enabled()
    call check_node(node)
    if (child < 1 .or. child > node_child_counts(node)) then
      call fail_metadata('decay child is outside its node')
    end if
  end subroutine validate_child_lookup


  subroutine require_initialized()
    if (.not. initialized) call initialize_decay_chain_metadata()
  end subroutine require_initialized


  subroutine require_enabled()
    call require_initialized()
    if (.not. has_decay_chains()) then
      call fail_metadata('no active decay-chain metadata are present')
    end if
  end subroutine require_enabled


  subroutine check_node(node)
    integer, intent(in) :: node
    if (node < 1 .or. node > number_of_nodes) then
      call fail_metadata('decay-node index is out of range')
    end if
  end subroutine check_node


  subroutine check_leaf(leaf)
    integer, intent(in) :: leaf
    if (leaf < 1 .or. leaf > number_of_leaves) then
      call fail_metadata('decay-leaf index is out of range')
    end if
  end subroutine check_leaf


  subroutine check_context(context)
    integer, intent(in) :: context
    if (context < 1 .or. context > number_of_contexts) then
      call fail_metadata('decay context is out of range')
    end if
  end subroutine check_context


  subroutine check_core_leg_index(core_leg)
    integer, intent(in) :: core_leg
    if (core_leg < 1 .or. core_leg > nexternal) then
      call fail_metadata('CORE_LEG index is out of range')
    end if
  end subroutine check_core_leg_index


  subroutine check_configuration(configuration)
    integer, intent(in) :: configuration
    if (configuration < 1 .or. configuration > number_of_fks_maps) then
      call fail_metadata('FKS configuration is out of range')
    end if
  end subroutine check_configuration


  subroutine fail_metadata(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in decay_chain_metadata: '//trim(message)
    stop 1
  end subroutine fail_metadata

end module decay_chain_metadata
