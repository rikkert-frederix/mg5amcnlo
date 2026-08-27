module nlo_decay_metadata
  use process_dimensions, only: nexternal, nincoming, fks_configs, &
                                validate_process_dimensions
  use fks_metadata, only: fks_i_d, fks_j_d, validate_fks_metadata
  implicit none
  private

  integer, parameter, public :: nlo_decay_leg_target = 1
  integer, parameter, public :: nlo_decay_node_target = 2
  integer, parameter, public :: nlo_decay_leaf_child = 1
  integer, parameter, public :: nlo_decay_node_child = 2

  logical, save :: initialized = .false.
  logical, save :: enabled = .false.
  integer, save :: metadata_format = 0
  integer, save :: corrected_parent_pdg_value = 0
  integer, save :: corrected_parent_occurrence_value = 0
  integer, save :: number_of_contexts = 0
  integer, save :: number_of_fks_maps = 0
  integer, save :: number_of_fks_partners = 0
  integer, save :: number_of_generated_color_links = 0
  integer, save :: number_of_color_link_records = 0
  integer, save :: number_of_production_legs = 0
  integer, save :: number_of_nodes = 0
  integer, save :: number_of_leaves = 0
  integer, save :: corrected_node_value = 0
  integer, save :: born_context_id = 0
  integer, save :: production_born_qcd_order_value = -1
  integer, save :: decay_born_qcd_order_value = -1

  integer, allocatable, save :: production_pdg_values(:)
  logical, allocatable, save :: production_final_values(:)
  character(len=16), allocatable, save :: context_kind_values(:)
  integer, allocatable, save :: context_source_values(:)
  integer, allocatable, save :: context_local_counts(:)
  integer, allocatable, save :: context_visible_counts(:)
  integer, allocatable, save :: production_target_kinds(:, :)
  integer, allocatable, save :: production_target_ids(:, :)
  integer, allocatable, save :: local_pdg_values(:, :)
  logical, allocatable, save :: local_final_values(:, :)
  integer, allocatable, save :: local_target_kinds(:, :)
  integer, allocatable, save :: local_target_ids(:, :)
  integer, allocatable, save :: fks_context_values(:)
  integer, allocatable, save :: fks_i_values(:)
  integer, allocatable, save :: fks_j_values(:)
  integer, allocatable, save :: fks_ij_values(:)
  integer, allocatable, save :: fks_target_kinds(:, :)
  integer, allocatable, save :: fks_target_ids(:, :)
  integer, allocatable, save :: fks_target_local_legs(:, :)
  integer, allocatable, save :: fks_partner_target_kinds(:, :)
  integer, allocatable, save :: fks_partner_target_ids(:, :)
  integer, allocatable, save :: fks_partner_counts(:)
  integer, allocatable, save :: fks_partner_local_values(:, :)
  integer, allocatable, save :: real_to_born_values(:, :)
  integer, allocatable, save :: color_local_first_values(:)
  integer, allocatable, save :: color_local_second_values(:)
  integer, allocatable, save :: color_visible_first_values(:)
  integer, allocatable, save :: color_visible_second_values(:)
  integer, allocatable, save :: color_generated_index_values(:)
  integer, allocatable, save :: node_parent_values(:)
  integer, allocatable, save :: node_pdg_values(:)
  integer, allocatable, save :: node_qcd_orders(:)
  integer, allocatable, save :: node_carrier_values(:)
  integer, allocatable, save :: node_child_counts(:)
  integer, allocatable, save :: node_child_kinds(:, :)
  integer, allocatable, save :: node_child_ids(:, :)
  integer, allocatable, save :: leaf_parent_values(:)
  integer, allocatable, save :: leaf_pdg_values(:)
  integer, allocatable, save :: leaf_visible_values(:, :)

  public :: initialize_nlo_decay_metadata, has_nlo_decay
  public :: corrected_parent_pdg, corrected_parent_occurrence
  public :: nlo_decay_production_born_qcd_order
  public :: nlo_decay_born_qcd_order
  public :: nlo_decay_born_context, nlo_decay_context_for_fks
  public :: nlo_decay_production_count, nlo_decay_production_pdg
  public :: nlo_decay_production_is_final
  public :: nlo_decay_production_target_kind
  public :: nlo_decay_production_target_id
  public :: nlo_decay_local_count, nlo_decay_local_pdg
  public :: nlo_decay_local_is_final, nlo_decay_local_target_kind
  public :: nlo_decay_local_target_id, nlo_decay_visible_count
  public :: nlo_decay_fks_i, nlo_decay_fks_j, nlo_decay_fks_ij
  public :: nlo_decay_partner_count, nlo_decay_partner_local
  public :: nlo_decay_map_color_link
  public :: nlo_decay_real_to_born
  public :: nlo_decay_corrected_node, nlo_decay_node_count
  public :: nlo_decay_leaf_count, nlo_decay_node_pdg
  public :: nlo_decay_node_qcd_order, nlo_decay_node_child_count
  public :: nlo_decay_node_child_kind, nlo_decay_node_child_id
  public :: nlo_decay_leaf_pdg, nlo_decay_leaf_visible

contains

  subroutine initialize_nlo_decay_metadata()
    logical :: exists, end_seen
    integer :: unit_number, ios, production_records, color_record
    integer :: context, configuration, leg, pdg, source, local_count
    integer :: visible_count, local_i, local_j, local_ij, target
    integer :: target_position, node, leaf, parent, qcd_order, carrier
    integer :: child, child_count, visible
    character(len=512) :: line
    character(len=32) :: keyword, kind, state, name
    character(len=16) :: child_words(nexternal)
    integer :: child_ids(nexternal)

    if (initialized) return
    call validate_process_dimensions()
    call validate_fks_metadata()
    inquire(file='nlo_decay_info.dat', exist=exists)
    if (.not. exists) then
      initialized = .true.
      enabled = .false.
      return
    end if

    open(newunit=unit_number, file='nlo_decay_info.dat', status='old', &
         action='read', iostat=ios)
    if (ios /= 0) call fail_metadata('cannot open nlo_decay_info.dat')

    production_records = 0
    number_of_color_link_records = 0
    do
      read(unit_number, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) call fail_metadata('cannot read nlo_decay_info.dat')
      if (len_trim(line) == 0) cycle
      read(line, *, iostat=ios) keyword
      if (ios /= 0) call fail_metadata('malformed metadata record')
      select case (trim(keyword))
      case ('FORMAT')
        read(line, *, iostat=ios) keyword, metadata_format
      case ('PARENT')
        read(line, *, iostat=ios) keyword, corrected_parent_pdg_value, &
             corrected_parent_occurrence_value
      case ('COUNTS')
        read(line, *, iostat=ios) keyword, number_of_contexts, &
             number_of_fks_maps, number_of_generated_color_links, &
             number_of_fks_partners
      case ('QCD_ORDERS')
        read(line, *, iostat=ios) keyword, &
             production_born_qcd_order_value, &
             decay_born_qcd_order_value
      case ('TOPOLOGY')
        read(line, *, iostat=ios) keyword, number_of_nodes, &
             number_of_leaves, corrected_node_value
      case ('PRODUCTION_LEG')
        production_records = production_records + 1
      case ('COLOR_LINK')
        number_of_color_link_records = number_of_color_link_records + 1
      end select
      if (ios /= 0) call fail_metadata('malformed metadata header')
    end do

    if (metadata_format /= 5) then
      call fail_metadata('FORMAT 5 is required; regenerate the process')
    end if
    if (corrected_parent_pdg_value == 0 .or. &
        corrected_parent_occurrence_value < 1) then
      call fail_metadata('invalid corrected-decay parent')
    end if
    if (number_of_contexts < 2 .or. &
        number_of_fks_maps /= fks_configs .or. &
        number_of_fks_partners < number_of_fks_maps .or. &
        number_of_generated_color_links < 1 .or. &
        number_of_color_link_records < number_of_generated_color_links) then
      call fail_metadata('invalid metadata counts')
    end if
    if (production_records < nincoming + 1 .or. &
        production_records > nexternal) then
      call fail_metadata('invalid production-leg count')
    end if
    if (number_of_nodes < 1 .or. number_of_leaves < 1 .or. &
        corrected_node_value < 1 .or. &
        corrected_node_value > number_of_nodes) then
      call fail_metadata('invalid decay topology counts')
    end if
    number_of_production_legs = production_records

    allocate(production_pdg_values(number_of_production_legs))
    allocate(production_final_values(number_of_production_legs))
    allocate(context_kind_values(number_of_contexts))
    allocate(context_source_values(number_of_contexts))
    allocate(context_local_counts(number_of_contexts))
    allocate(context_visible_counts(number_of_contexts))
    allocate(production_target_kinds(number_of_production_legs, &
                                     number_of_contexts))
    allocate(production_target_ids(number_of_production_legs, &
                                   number_of_contexts))
    allocate(local_pdg_values(nexternal, number_of_contexts))
    allocate(local_final_values(nexternal, number_of_contexts))
    allocate(local_target_kinds(nexternal, number_of_contexts))
    allocate(local_target_ids(nexternal, number_of_contexts))
    allocate(fks_context_values(fks_configs))
    allocate(fks_i_values(fks_configs))
    allocate(fks_j_values(fks_configs))
    allocate(fks_ij_values(fks_configs))
    allocate(fks_target_kinds(3, fks_configs))
    allocate(fks_target_ids(3, fks_configs))
    allocate(fks_target_local_legs(3, fks_configs))
    allocate(fks_partner_target_kinds(nexternal, fks_configs))
    allocate(fks_partner_target_ids(nexternal, fks_configs))
    allocate(fks_partner_counts(fks_configs))
    allocate(fks_partner_local_values(nexternal, fks_configs))
    allocate(real_to_born_values(nexternal, fks_configs))
    allocate(color_local_first_values(number_of_color_link_records))
    allocate(color_local_second_values(number_of_color_link_records))
    allocate(color_visible_first_values(number_of_color_link_records))
    allocate(color_visible_second_values(number_of_color_link_records))
    allocate(color_generated_index_values(number_of_color_link_records))
    allocate(node_parent_values(number_of_nodes))
    allocate(node_pdg_values(number_of_nodes))
    allocate(node_qcd_orders(number_of_nodes))
    allocate(node_carrier_values(number_of_nodes))
    allocate(node_child_counts(number_of_nodes))
    allocate(node_child_kinds(nexternal, number_of_nodes))
    allocate(node_child_ids(nexternal, number_of_nodes))
    allocate(leaf_parent_values(number_of_leaves))
    allocate(leaf_pdg_values(number_of_leaves))
    allocate(leaf_visible_values(number_of_leaves, number_of_contexts))

    production_pdg_values = 0
    production_final_values = .false.
    context_kind_values = ''
    context_source_values = 0
    context_local_counts = 0
    context_visible_counts = 0
    production_target_kinds = 0
    production_target_ids = 0
    local_pdg_values = 0
    local_final_values = .false.
    local_target_kinds = 0
    local_target_ids = 0
    fks_context_values = 0
    fks_i_values = 0
    fks_j_values = 0
    fks_ij_values = 0
    fks_target_kinds = 0
    fks_target_ids = 0
    fks_target_local_legs = 0
    fks_partner_target_kinds = 0
    fks_partner_target_ids = 0
    fks_partner_counts = 0
    fks_partner_local_values = 0
    real_to_born_values = 0
    color_local_first_values = 0
    color_local_second_values = 0
    color_visible_first_values = 0
    color_visible_second_values = 0
    color_generated_index_values = 0
    node_parent_values = -1
    node_pdg_values = 0
    node_qcd_orders = -1
    node_carrier_values = 0
    node_child_counts = 0
    node_child_kinds = 0
    node_child_ids = 0
    leaf_parent_values = 0
    leaf_pdg_values = 0
    leaf_visible_values = 0

    rewind(unit_number)
    end_seen = .false.
    color_record = 0
    do
      read(unit_number, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) call fail_metadata('cannot read metadata body')
      if (len_trim(line) == 0) cycle
      if (end_seen) call fail_metadata('record found after END')
      read(line, *, iostat=ios) keyword
      if (ios /= 0) call fail_metadata('malformed metadata keyword')
      select case (trim(keyword))
      case ('FORMAT', 'STATUS', 'CORRECTION', 'PARENT', 'HAS_VIRTUAL', &
            'VIRTUAL_COMPOSITION', 'VIRTUAL_CURRENT_COUNT', &
            'QCD_ORDERS', 'FORCED_SPECIES', 'TOPOLOGY', 'COUNTS')
        continue
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
            node_child_kinds(child, node) = nlo_decay_leaf_child
          case ('NODE')
            node_child_kinds(child, node) = nlo_decay_node_child
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
      case ('PRODUCTION_LEG')
        read(line, *, iostat=ios) keyword, leg, pdg, state
        call check_production_leg(leg)
        if (production_pdg_values(leg) /= 0) then
          call fail_metadata('duplicate PRODUCTION_LEG record')
        end if
        production_pdg_values(leg) = pdg
        call parse_state(state, production_final_values(leg))
      case ('CONTEXT')
        read(line, *, iostat=ios) keyword, context, kind, source, &
             local_count, visible_count
        call check_context(context)
        if (context_local_counts(context) /= 0) then
          call fail_metadata('duplicate CONTEXT record')
        end if
        context_kind_values(context) = kind(1:16)
        context_source_values(context) = source
        context_local_counts(context) = local_count
        context_visible_counts(context) = visible_count
      case ('PRODUCTION_MAP')
        read(line, *, iostat=ios) keyword, context, leg, kind, target
        call check_context(context)
        call check_production_leg(leg)
        if (production_target_kinds(leg, context) /= 0) then
          call fail_metadata('duplicate PRODUCTION_MAP record')
        end if
        production_target_kinds(leg, context) = target_kind(kind)
        production_target_ids(leg, context) = target
      case ('LOCAL_LEG')
        read(line, *, iostat=ios) keyword, context, leg, pdg, state
        call check_context(context)
        call check_local_storage_index(leg)
        if (local_pdg_values(leg, context) /= 0) then
          call fail_metadata('duplicate LOCAL_LEG record')
        end if
        local_pdg_values(leg, context) = pdg
        call parse_state(state, local_final_values(leg, context))
      case ('LOCAL_MAP')
        read(line, *, iostat=ios) keyword, context, leg, kind, target
        call check_context(context)
        call check_local_storage_index(leg)
        if (local_target_kinds(leg, context) /= 0) then
          call fail_metadata('duplicate LOCAL_MAP record')
        end if
        local_target_kinds(leg, context) = target_kind(kind)
        local_target_ids(leg, context) = target
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
             local_i, local_j, local_ij
        call check_configuration(configuration)
        if (fks_context_values(configuration) /= 0) then
          call fail_metadata('duplicate FKS_MAP record')
        end if
        fks_context_values(configuration) = context
        fks_i_values(configuration) = local_i
        fks_j_values(configuration) = local_j
        fks_ij_values(configuration) = local_ij
      case ('FKS_TARGET')
        read(line, *, iostat=ios) keyword, configuration, name, leg, &
             kind, target
        call check_configuration(configuration)
        select case (trim(name))
        case ('I')
          target_position = 1
        case ('J')
          target_position = 2
        case ('IJ')
          target_position = 3
        case default
          call fail_metadata('unknown FKS target name')
        end select
        if (fks_target_kinds(target_position, configuration) /= 0) then
          call fail_metadata('duplicate FKS_TARGET record')
        end if
        fks_target_kinds(target_position, configuration) = target_kind(kind)
        fks_target_ids(target_position, configuration) = target
        fks_target_local_legs(target_position, configuration) = leg
      case ('FKS_PARTNER')
        read(line, *, iostat=ios) keyword, configuration, leg, kind, target
        call check_configuration(configuration)
        call check_local_storage_index(leg)
        if (fks_partner_target_kinds(leg, configuration) /= 0) then
          call fail_metadata('duplicate FKS_PARTNER record')
        end if
        fks_partner_target_kinds(leg, configuration) = target_kind(kind)
        fks_partner_target_ids(leg, configuration) = target
        fks_partner_counts(configuration) = &
             fks_partner_counts(configuration) + 1
        if (fks_partner_counts(configuration) > nexternal) then
          call fail_metadata('too many FKS partners in one configuration')
        end if
        fks_partner_local_values(fks_partner_counts(configuration), &
                                 configuration) = leg
      case ('REAL_BORN_MAP')
        read(line, *, iostat=ios) keyword, configuration, leg, target
        call check_configuration(configuration)
        call check_local_storage_index(leg)
        if (real_to_born_values(leg, configuration) /= 0) then
          call fail_metadata('duplicate REAL_BORN_MAP record')
        end if
        real_to_born_values(leg, configuration) = target
      case ('COLOR_LINK')
        color_record = color_record + 1
        if (color_record > number_of_color_link_records) then
          call fail_metadata('too many COLOR_LINK records')
        end if
        read(line, *, iostat=ios) keyword, &
             color_local_first_values(color_record), &
             color_local_second_values(color_record), &
             color_visible_first_values(color_record), &
             color_visible_second_values(color_record), &
             color_generated_index_values(color_record)
      case ('END')
        end_seen = .true.
      case default
        call fail_metadata('unknown metadata keyword '//trim(keyword))
      end select
      if (ios /= 0) call fail_metadata('malformed metadata record')
    end do
    close(unit_number)
    if (.not. end_seen) call fail_metadata('END record is absent')
    if (color_record /= number_of_color_link_records) then
      call fail_metadata('COLOR_LINK record count changed while reading')
    end if

    call validate_metadata_values()
    enabled = .true.
    initialized = .true.
  end subroutine initialize_nlo_decay_metadata


  subroutine validate_metadata_values()
    integer :: leg, context, configuration, target, born_count
    integer :: node_maps, expected_visible, node, child, identifier, leaf
    logical :: visible_used(nexternal), born_local_used(nexternal)
    logical :: root_used(number_of_nodes)
    logical :: node_referenced(number_of_nodes)
    logical :: leaf_referenced(number_of_leaves)

    node_referenced = .false.
    leaf_referenced = .false.
    do node = 1, number_of_nodes
      if (node_parent_values(node) < 0 .or. node_pdg_values(node) == 0 .or. &
          node_qcd_orders(node) < 0 .or. node_child_counts(node) < 2) then
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
        case (nlo_decay_node_child)
          call check_node(identifier)
          if (identifier <= node .or. &
              node_parent_values(identifier) /= node .or. &
              node_referenced(identifier)) then
            call fail_metadata('invalid nested decay-node reference')
          end if
          node_referenced(identifier) = .true.
        case (nlo_decay_leaf_child)
          call check_leaf(identifier)
          if (leaf_parent_values(identifier) /= node .or. &
              leaf_referenced(identifier)) then
            call fail_metadata('invalid decay-leaf reference')
          end if
          leaf_referenced(identifier) = .true.
        case default
          call fail_metadata('invalid decay-node child kind')
        end select
      end do
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
    if (node_pdg_values(corrected_node_value) /= &
        corrected_parent_pdg_value) then
      call fail_metadata('corrected node and PARENT records disagree')
    end if

    if (production_born_qcd_order_value < 0 .or. &
        decay_born_qcd_order_value < 0 .or. &
        mod(production_born_qcd_order_value, 2) /= 0 .or. &
        mod(decay_born_qcd_order_value, 2) /= 0) then
      call fail_metadata('invalid Born QCD orders')
    end if
    if (any(production_pdg_values == 0)) then
      call fail_metadata('production-leg records are incomplete')
    end if
    if (count(fks_partner_target_kinds /= 0) /= number_of_fks_partners) then
      call fail_metadata('FKS-partner records disagree with COUNTS')
    end if
    if (sum(fks_partner_counts) /= number_of_fks_partners) then
      call fail_metadata('ordered FKS-partner records disagree with COUNTS')
    end if
    do leg = 1, number_of_production_legs
      if ((leg <= nincoming) .eqv. production_final_values(leg)) then
        call fail_metadata('production-leg state disagrees with NINCOMING')
      end if
    end do

    born_count = 0
    do context = 1, number_of_contexts
      if (context_source_values(context) < 1 .or. &
          context_local_counts(context) < 2 .or. &
          context_local_counts(context) > nexternal) then
        call fail_metadata('invalid context extent')
      end if
      select case (trim(context_kind_values(context)))
      case ('BORN')
        born_count = born_count + 1
        born_context_id = context
        expected_visible = nexternal - 1
      case ('REAL')
        expected_visible = nexternal
      case default
        call fail_metadata('invalid context kind')
      end select
      if (context_visible_counts(context) /= expected_visible) then
        call fail_metadata('context visible count has the wrong extent')
      end if
      if (any(production_target_kinds(:, context) == 0) .or. &
          any(production_target_ids(:, context) == 0)) then
        call fail_metadata('production maps are incomplete')
      end if
      if (any(local_pdg_values(1:context_local_counts(context), context) == 0) &
          .or. any(local_target_kinds(1:context_local_counts(context), &
                                     context) == 0) .or. &
          any(local_target_ids(1:context_local_counts(context), context) == 0)) then
        call fail_metadata('local decay records are incomplete')
      end if
      if (context_local_counts(context) < nexternal) then
        if (any(local_pdg_values(context_local_counts(context) + 1:, &
                                 context) /= 0)) then
          call fail_metadata('local decay record lies beyond its context')
        end if
      end if

      visible_used = .false.
      root_used = .false.
      do leg = 1, number_of_production_legs
        target = production_target_ids(leg, context)
        select case (production_target_kinds(leg, context))
        case (nlo_decay_leg_target)
          call mark_visible_target(target, context_visible_counts(context), &
                                   visible_used)
        case (nlo_decay_node_target)
          call check_node(target)
          if (.not. production_final_values(leg) .or. &
              node_parent_values(target) /= 0 .or. &
              production_pdg_values(leg) /= node_pdg_values(target) .or. &
              root_used(target)) then
            call fail_metadata('invalid production decay-node target')
          end if
          root_used(target) = .true.
        case default
          call fail_metadata('invalid production target kind')
        end select
      end do
      do node = 1, number_of_nodes
        if (node_parent_values(node) == 0 .and. .not. root_used(node)) then
          call fail_metadata('a root decay node is absent from the context')
        end if
      end do

      node_maps = 0
      do leg = 1, context_local_counts(context)
        target = local_target_ids(leg, context)
        select case (local_target_kinds(leg, context))
        case (nlo_decay_leg_target)
          if (.not. local_final_values(leg, context)) then
            call fail_metadata('incoming decay leg maps to a visible leg')
          end if
          call mark_visible_target(target, context_visible_counts(context), &
                                   visible_used)
        case (nlo_decay_node_target)
          call check_node(target)
          if (local_final_values(leg, context)) then
            if (node_parent_values(target) /= corrected_node_value .or. &
                local_pdg_values(leg, context) /= node_pdg_values(target)) then
              call fail_metadata('invalid local nested-node target')
            end if
          else if (target /= corrected_node_value .or. &
                   local_pdg_values(leg, context) /= &
                   corrected_parent_pdg_value) then
            call fail_metadata('invalid local corrected-parent target')
          else
            node_maps = node_maps + 1
          end if
        case default
          call fail_metadata('invalid local target kind')
        end select
      end do
      if (node_maps /= 1) then
        call fail_metadata('corrected parent is not mapped exactly once')
      end if
      do leaf = 1, number_of_leaves
        target = leaf_visible_values(leaf, context)
        if (leaf_parent_values(leaf) == corrected_node_value) then
          if (target /= 0) then
            call fail_metadata('a direct corrected leaf has a LEAF_MAP')
          end if
          cycle
        end if
        if (target == 0) then
          call fail_metadata('a static decay leaf has no visible target')
        end if
        call mark_visible_target(target, context_visible_counts(context), &
                                 visible_used)
      end do
      if (.not. all(visible_used(1:context_visible_counts(context)))) then
        call fail_metadata('context does not cover the visible event')
      end if
    end do
    if (born_count /= 1) then
      call fail_metadata('exactly one Born context is required')
    end if
    call validate_color_link_values()

    do configuration = 1, fks_configs
      context = fks_context_values(configuration)
      call check_context(context)
      if (trim(context_kind_values(context)) /= 'REAL') then
        call fail_metadata('FKS map does not reference a real context')
      end if
      call check_local_leg(context, fks_i_values(configuration))
      call check_local_leg(context, fks_j_values(configuration))
      call check_local_leg(born_context_id, fks_ij_values(configuration))
      if (context_local_counts(context) /= &
          context_local_counts(born_context_id) + 1) then
        call fail_metadata('real decay context must add exactly one leg')
      end if
      if (fks_i_values(configuration) == fks_j_values(configuration)) then
        call fail_metadata('FKS i and j are identical')
      end if
      if (.not. local_final_values(fks_i_values(configuration), context) .or. &
          .not. local_final_values(fks_j_values(configuration), context) .or. &
          .not. local_final_values(fks_ij_values(configuration), &
                                   born_context_id)) then
        call fail_metadata('decay-local FKS legs must be final state')
      end if
      if (any(fks_target_kinds(:, configuration) /= &
              nlo_decay_leg_target)) then
        call fail_metadata('FKS i, j and ij must be visible targets')
      end if
      if (any(fks_target_local_legs(:, configuration) /= &
              (/fks_i_values(configuration), fks_j_values(configuration), &
                fks_ij_values(configuration)/))) then
        call fail_metadata('FKS target labels disagree with the local map')
      end if
      if (fks_target_ids(1, configuration) /= &
          local_target_ids(fks_i_values(configuration), context) .or. &
          fks_target_ids(2, configuration) /= &
          local_target_ids(fks_j_values(configuration), context) .or. &
          fks_target_ids(3, configuration) /= &
          local_target_ids(fks_ij_values(configuration), born_context_id)) then
        call fail_metadata('FKS targets disagree with the context maps')
      end if
      if (fks_target_ids(1, configuration) /= fks_i_d(configuration) .or. &
          fks_target_ids(2, configuration) /= fks_j_d(configuration)) then
        call fail_metadata('local and generated visible FKS maps disagree')
      end if
      do leg = 1, context_local_counts(context)
        if (fks_partner_target_kinds(leg, configuration) == 0) cycle
        if (fks_partner_target_kinds(leg, configuration) /= &
            local_target_kinds(leg, context) .or. &
            fks_partner_target_ids(leg, configuration) /= &
            local_target_ids(leg, context)) then
          call fail_metadata('FKS partner disagrees with the local map')
        end if
      end do
      if (context_local_counts(context) < nexternal) then
        if (any(fks_partner_target_kinds( &
                context_local_counts(context) + 1:, configuration) /= 0)) then
          call fail_metadata('FKS partner lies beyond its real context')
        end if
      end if
      if (fks_partner_target_kinds( &
          fks_j_values(configuration), configuration) == 0) then
        call fail_metadata('FKS sister is absent from the partner list')
      end if
      if (fks_partner_counts(configuration) < 1) then
        call fail_metadata('FKS configuration has no soft partners')
      end if
      do target = 1, fks_partner_counts(configuration)
        leg = fks_partner_local_values(target, configuration)
        call check_local_leg(context, leg)
        if (leg == fks_i_values(configuration)) then
          call fail_metadata('radiated leg appears in its soft-partner list')
        end if
        if (fks_partner_target_kinds(leg, configuration) == 0) then
          call fail_metadata('ordered FKS partner is absent from its map')
        end if
      end do
      born_local_used = .false.
      do leg = 1, context_local_counts(context)
        if (leg == fks_i_values(configuration)) then
          if (real_to_born_values(leg, configuration) /= 0) then
            call fail_metadata('radiated leg has a Born-map target')
          end if
        else
          target = real_to_born_values(leg, configuration)
          call check_local_leg(born_context_id, target)
          if (born_local_used(target)) then
            call fail_metadata('real-to-Born map is not one-to-one')
          end if
          born_local_used(target) = .true.
        end if
      end do
      if (.not. all(born_local_used( &
          1:context_local_counts(born_context_id)))) then
        call fail_metadata('real-to-Born map does not cover the Born decay')
      end if
      if (real_to_born_values(fks_j_values(configuration), configuration) /= &
          fks_ij_values(configuration)) then
        call fail_metadata('decay-local j does not map to Born ij')
      end if
    end do
  end subroutine validate_metadata_values


  subroutine validate_color_link_values()
    integer :: record, other, local_first, local_second
    integer :: visible_first, visible_second, generated_index
    integer :: expected_first, expected_second, endpoint_first
    integer :: endpoint_second, node_leg, node_carrier, direct_target
    logical :: node_used
    logical, allocatable :: generated_used(:)

    allocate(generated_used(number_of_generated_color_links))
    generated_used = .false.
    node_leg = 0
    do local_first = 1, context_local_counts(born_context_id)
      if (local_target_kinds(local_first, born_context_id) == &
          nlo_decay_node_target .and. &
          .not. local_final_values(local_first, born_context_id)) then
        node_leg = local_first
      end if
    end do
    if (node_leg == 0) then
      call fail_metadata('Born decay context has no parent node')
    end if

    ! A massive coloured parent has a self link.  Its visible image fixes the
    ! unique colour carrier used when the decay is embedded in production.
    node_carrier = 0
    node_used = .false.
    do record = 1, number_of_color_link_records
      local_first = color_local_first_values(record)
      local_second = color_local_second_values(record)
      visible_first = color_visible_first_values(record)
      visible_second = color_visible_second_values(record)
      if (local_first /= node_leg .and. local_second /= node_leg) cycle
      node_used = .true.
      if (local_first == node_leg .and. local_second == node_leg) then
        if (visible_first /= visible_second) then
          call fail_metadata('parent self link maps to two visible legs')
        end if
        call set_node_carrier(visible_first, node_carrier)
      else
        if (local_first == node_leg) then
          direct_target = local_target_ids(local_second, born_context_id)
        else
          direct_target = local_target_ids(local_first, born_context_id)
        end if
        if (visible_first == direct_target) then
          call set_node_carrier(visible_second, node_carrier)
        else if (visible_second == direct_target) then
          call set_node_carrier(visible_first, node_carrier)
        else
          call fail_metadata('parent color link omits its direct-leg target')
        end if
      end if
    end do
    if (node_used .and. node_carrier == 0) then
      call fail_metadata('cannot identify the parent color carrier')
    end if

    do record = 1, number_of_color_link_records
      local_first = color_local_first_values(record)
      local_second = color_local_second_values(record)
      visible_first = color_visible_first_values(record)
      visible_second = color_visible_second_values(record)
      generated_index = color_generated_index_values(record)
      call check_local_leg(born_context_id, local_first)
      call check_local_leg(born_context_id, local_second)
      if (local_first > local_second) then
        call fail_metadata('COLOR_LINK local pair is not ordered')
      end if
      if (visible_first < 1 .or. visible_second > nexternal - 1 .or. &
          visible_first > visible_second) then
        call fail_metadata('COLOR_LINK visible pair is invalid')
      end if
      if (generated_index < 1 .or. &
          generated_index > number_of_generated_color_links) then
        call fail_metadata('COLOR_LINK generated index is out of range')
      end if
      generated_used(generated_index) = .true.

      if (local_first == node_leg) then
        endpoint_first = node_carrier
      else
        if (local_target_kinds(local_first, born_context_id) /= &
            nlo_decay_leg_target) then
          call fail_metadata('COLOR_LINK endpoint has an invalid target')
        end if
        endpoint_first = local_target_ids(local_first, born_context_id)
      end if
      if (local_second == node_leg) then
        endpoint_second = node_carrier
      else
        if (local_target_kinds(local_second, born_context_id) /= &
            nlo_decay_leg_target) then
          call fail_metadata('COLOR_LINK endpoint has an invalid target')
        end if
        endpoint_second = local_target_ids(local_second, born_context_id)
      end if
      expected_first = min(endpoint_first, endpoint_second)
      expected_second = max(endpoint_first, endpoint_second)
      if (visible_first /= expected_first .or. &
          visible_second /= expected_second) then
        call fail_metadata('COLOR_LINK disagrees with the decay color map')
      end if

      do other = 1, record - 1
        if (color_local_first_values(other) == local_first .and. &
            color_local_second_values(other) == local_second) then
          call fail_metadata('duplicate local COLOR_LINK pair')
        end if
        if (color_generated_index_values(other) == generated_index) then
          if (color_visible_first_values(other) /= visible_first .or. &
              color_visible_second_values(other) /= visible_second) then
            call fail_metadata('generated COLOR_LINK index is ambiguous')
          end if
        else if (color_visible_first_values(other) == visible_first .and. &
                 color_visible_second_values(other) == visible_second) then
          call fail_metadata('visible COLOR_LINK pair has two indices')
        end if
      end do
    end do
    if (.not. all(generated_used)) then
      call fail_metadata('generated COLOR_LINK indices are incomplete')
    end if
    deallocate(generated_used)
  end subroutine validate_color_link_values


  subroutine set_node_carrier(candidate, node_carrier)
    integer, intent(in) :: candidate
    integer, intent(inout) :: node_carrier
    if (candidate < 1 .or. candidate > nexternal - 1) then
      call fail_metadata('parent color carrier is out of range')
    end if
    if (node_carrier /= 0 .and. node_carrier /= candidate) then
      call fail_metadata('parent has inconsistent visible color carriers')
    end if
    node_carrier = candidate
  end subroutine set_node_carrier


  subroutine mark_visible_target(target, visible_count, used)
    integer, intent(in) :: target, visible_count
    logical, intent(inout) :: used(:)
    if (target < 1 .or. target > visible_count) then
      call fail_metadata('visible target is out of range')
    end if
    if (used(target)) call fail_metadata('visible target is mapped twice')
    used(target) = .true.
  end subroutine mark_visible_target


  subroutine parse_state(word, is_final)
    character(len=*), intent(in) :: word
    logical, intent(out) :: is_final
    select case (trim(word))
    case ('I')
      is_final = .false.
    case ('F')
      is_final = .true.
    case default
      call fail_metadata('particle state must be I or F')
    end select
  end subroutine parse_state


  integer function target_kind(word)
    character(len=*), intent(in) :: word
    select case (trim(word))
    case ('LEG')
      target_kind = nlo_decay_leg_target
    case ('NODE')
      target_kind = nlo_decay_node_target
    case default
      call fail_metadata('target kind must be LEG or NODE')
    end select
  end function target_kind


  logical function has_nlo_decay()
    if (.not. initialized) call initialize_nlo_decay_metadata()
    has_nlo_decay = enabled
  end function has_nlo_decay


  integer function corrected_parent_pdg()
    call require_enabled()
    corrected_parent_pdg = corrected_parent_pdg_value
  end function corrected_parent_pdg


  integer function corrected_parent_occurrence()
    call require_enabled()
    corrected_parent_occurrence = corrected_parent_occurrence_value
  end function corrected_parent_occurrence


  integer function nlo_decay_production_born_qcd_order()
    call require_enabled()
    nlo_decay_production_born_qcd_order = &
         production_born_qcd_order_value
  end function nlo_decay_production_born_qcd_order


  integer function nlo_decay_born_qcd_order()
    call require_enabled()
    nlo_decay_born_qcd_order = decay_born_qcd_order_value
  end function nlo_decay_born_qcd_order


  integer function nlo_decay_born_context()
    call require_enabled()
    nlo_decay_born_context = born_context_id
  end function nlo_decay_born_context


  integer function nlo_decay_context_for_fks(configuration)
    integer, intent(in) :: configuration
    call require_enabled()
    call check_configuration(configuration)
    nlo_decay_context_for_fks = fks_context_values(configuration)
  end function nlo_decay_context_for_fks


  integer function nlo_decay_production_count()
    call require_enabled()
    nlo_decay_production_count = number_of_production_legs
  end function nlo_decay_production_count


  integer function nlo_decay_production_pdg(leg)
    integer, intent(in) :: leg
    call require_enabled()
    call check_production_leg(leg)
    nlo_decay_production_pdg = production_pdg_values(leg)
  end function nlo_decay_production_pdg


  logical function nlo_decay_production_is_final(leg)
    integer, intent(in) :: leg
    call require_enabled()
    call check_production_leg(leg)
    nlo_decay_production_is_final = production_final_values(leg)
  end function nlo_decay_production_is_final


  integer function nlo_decay_production_target_kind(context, leg)
    integer, intent(in) :: context, leg
    call require_enabled()
    call check_context(context)
    call check_production_leg(leg)
    nlo_decay_production_target_kind = production_target_kinds(leg, context)
  end function nlo_decay_production_target_kind


  integer function nlo_decay_production_target_id(context, leg)
    integer, intent(in) :: context, leg
    call require_enabled()
    call check_context(context)
    call check_production_leg(leg)
    nlo_decay_production_target_id = production_target_ids(leg, context)
  end function nlo_decay_production_target_id


  integer function nlo_decay_local_count(context)
    integer, intent(in) :: context
    call require_enabled()
    call check_context(context)
    nlo_decay_local_count = context_local_counts(context)
  end function nlo_decay_local_count


  integer function nlo_decay_local_pdg(context, leg)
    integer, intent(in) :: context, leg
    call require_enabled()
    call check_local_leg(context, leg)
    nlo_decay_local_pdg = local_pdg_values(leg, context)
  end function nlo_decay_local_pdg


  logical function nlo_decay_local_is_final(context, leg)
    integer, intent(in) :: context, leg
    call require_enabled()
    call check_local_leg(context, leg)
    nlo_decay_local_is_final = local_final_values(leg, context)
  end function nlo_decay_local_is_final


  integer function nlo_decay_local_target_kind(context, leg)
    integer, intent(in) :: context, leg
    call require_enabled()
    call check_local_leg(context, leg)
    nlo_decay_local_target_kind = local_target_kinds(leg, context)
  end function nlo_decay_local_target_kind


  integer function nlo_decay_local_target_id(context, leg)
    integer, intent(in) :: context, leg
    call require_enabled()
    call check_local_leg(context, leg)
    nlo_decay_local_target_id = local_target_ids(leg, context)
  end function nlo_decay_local_target_id


  integer function nlo_decay_visible_count(context)
    integer, intent(in) :: context
    call require_enabled()
    call check_context(context)
    nlo_decay_visible_count = context_visible_counts(context)
  end function nlo_decay_visible_count


  integer function nlo_decay_fks_i(configuration)
    integer, intent(in) :: configuration
    call require_enabled()
    call check_configuration(configuration)
    nlo_decay_fks_i = fks_i_values(configuration)
  end function nlo_decay_fks_i


  integer function nlo_decay_fks_j(configuration)
    integer, intent(in) :: configuration
    call require_enabled()
    call check_configuration(configuration)
    nlo_decay_fks_j = fks_j_values(configuration)
  end function nlo_decay_fks_j


  integer function nlo_decay_fks_ij(configuration)
    integer, intent(in) :: configuration
    call require_enabled()
    call check_configuration(configuration)
    nlo_decay_fks_ij = fks_ij_values(configuration)
  end function nlo_decay_fks_ij


  integer function nlo_decay_partner_count(configuration)
    integer, intent(in) :: configuration
    call require_enabled()
    call check_configuration(configuration)
    nlo_decay_partner_count = fks_partner_counts(configuration)
  end function nlo_decay_partner_count


  integer function nlo_decay_partner_local(configuration, position)
    integer, intent(in) :: configuration, position
    call require_enabled()
    call check_configuration(configuration)
    if (position < 1 .or. position > fks_partner_counts(configuration)) then
      call fail_metadata('FKS-partner position is out of range')
    end if
    nlo_decay_partner_local = &
         fks_partner_local_values(position, configuration)
  end function nlo_decay_partner_local


  subroutine nlo_decay_map_color_link(configuration, real_first, real_second, &
                                      visible_first, visible_second, multiplier)
    integer, intent(in) :: configuration, real_first, real_second
    integer, intent(out) :: visible_first, visible_second
    double precision, intent(out) :: multiplier
    integer :: born_first, born_second, record, temporary

    call require_enabled()
    call check_configuration(configuration)
    call check_local_leg(fks_context_values(configuration), real_first)
    call check_local_leg(fks_context_values(configuration), real_second)
    born_first = real_to_born_values(real_first, configuration)
    born_second = real_to_born_values(real_second, configuration)
    if (born_first == 0 .or. born_second == 0) then
      call fail_metadata('radiated leg has no decay color-link image')
    end if
    if (born_first > born_second) then
      temporary = born_first
      born_first = born_second
      born_second = temporary
    end if

    record = find_color_link_record(born_first, born_second)
    if (record == 0) then
      call fail_metadata('decay-local color link is absent from metadata')
    end if
    visible_first = color_visible_first_values(record)
    visible_second = color_visible_second_values(record)

    ! Generated self links contain the conventional factor 1/2.  Crossing
    ! an incoming decay charge to its unique final-state carrier changes its
    ! sign.  Together these reproduce the standalone decay color links when
    ! distinct local endpoints collapse onto one visible carrier.
    multiplier = 1d0
    if (.not. local_final_values(born_first, born_context_id)) then
      multiplier = -multiplier
    end if
    if (.not. local_final_values(born_second, born_context_id)) then
      multiplier = -multiplier
    end if
    if (born_first /= born_second .and. visible_first == visible_second) then
      multiplier = 2d0*multiplier
    end if
  end subroutine nlo_decay_map_color_link


  integer function find_color_link_record(local_first, local_second)
    integer, intent(in) :: local_first, local_second
    integer :: record
    find_color_link_record = 0
    do record = 1, number_of_color_link_records
      if (color_local_first_values(record) == local_first .and. &
          color_local_second_values(record) == local_second) then
        find_color_link_record = record
        return
      end if
    end do
  end function find_color_link_record


  integer function nlo_decay_real_to_born(configuration, real_leg)
    integer, intent(in) :: configuration, real_leg
    call require_enabled()
    call check_configuration(configuration)
    call check_local_leg(fks_context_values(configuration), real_leg)
    nlo_decay_real_to_born = real_to_born_values(real_leg, configuration)
  end function nlo_decay_real_to_born


  integer function nlo_decay_corrected_node()
    call require_enabled()
    nlo_decay_corrected_node = corrected_node_value
  end function nlo_decay_corrected_node


  integer function nlo_decay_node_count()
    call require_enabled()
    nlo_decay_node_count = number_of_nodes
  end function nlo_decay_node_count


  integer function nlo_decay_leaf_count()
    call require_enabled()
    nlo_decay_leaf_count = number_of_leaves
  end function nlo_decay_leaf_count


  integer function nlo_decay_node_pdg(node)
    integer, intent(in) :: node
    call require_enabled()
    call check_node(node)
    nlo_decay_node_pdg = node_pdg_values(node)
  end function nlo_decay_node_pdg


  integer function nlo_decay_node_qcd_order(node)
    integer, intent(in) :: node
    call require_enabled()
    call check_node(node)
    nlo_decay_node_qcd_order = node_qcd_orders(node)
  end function nlo_decay_node_qcd_order


  integer function nlo_decay_node_child_count(node)
    integer, intent(in) :: node
    call require_enabled()
    call check_node(node)
    nlo_decay_node_child_count = node_child_counts(node)
  end function nlo_decay_node_child_count


  integer function nlo_decay_node_child_kind(node, child)
    integer, intent(in) :: node, child
    call require_enabled()
    call check_node(node)
    if (child < 1 .or. child > node_child_counts(node)) then
      call fail_metadata('decay-node child is out of range')
    end if
    nlo_decay_node_child_kind = node_child_kinds(child, node)
  end function nlo_decay_node_child_kind


  integer function nlo_decay_node_child_id(node, child)
    integer, intent(in) :: node, child
    call require_enabled()
    call check_node(node)
    if (child < 1 .or. child > node_child_counts(node)) then
      call fail_metadata('decay-node child is out of range')
    end if
    nlo_decay_node_child_id = node_child_ids(child, node)
  end function nlo_decay_node_child_id


  integer function nlo_decay_leaf_pdg(leaf)
    integer, intent(in) :: leaf
    call require_enabled()
    call check_leaf(leaf)
    nlo_decay_leaf_pdg = leaf_pdg_values(leaf)
  end function nlo_decay_leaf_pdg


  integer function nlo_decay_leaf_visible(context, leaf)
    integer, intent(in) :: context, leaf
    call require_enabled()
    call check_context(context)
    call check_leaf(leaf)
    nlo_decay_leaf_visible = leaf_visible_values(leaf, context)
  end function nlo_decay_leaf_visible


  subroutine check_production_leg(leg)
    integer, intent(in) :: leg
    if (leg < 1 .or. leg > number_of_production_legs) then
      call fail_metadata('production leg is out of range')
    end if
  end subroutine check_production_leg


  subroutine check_node(node)
    integer, intent(in) :: node
    if (node < 1 .or. node > number_of_nodes) then
      call fail_metadata('decay node is out of range')
    end if
  end subroutine check_node


  subroutine check_leaf(leaf)
    integer, intent(in) :: leaf
    if (leaf < 1 .or. leaf > number_of_leaves) then
      call fail_metadata('decay leaf is out of range')
    end if
  end subroutine check_leaf


  subroutine check_context(context)
    integer, intent(in) :: context
    if (context < 1 .or. context > number_of_contexts) then
      call fail_metadata('context is out of range')
    end if
  end subroutine check_context


  subroutine check_local_storage_index(leg)
    integer, intent(in) :: leg
    if (leg < 1 .or. leg > nexternal) then
      call fail_metadata('local decay leg is out of storage range')
    end if
  end subroutine check_local_storage_index


  subroutine check_local_leg(context, leg)
    integer, intent(in) :: context, leg
    call check_context(context)
    if (leg < 1 .or. leg > context_local_counts(context)) then
      call fail_metadata('local decay leg is out of range')
    end if
  end subroutine check_local_leg


  subroutine check_configuration(configuration)
    integer, intent(in) :: configuration
    if (configuration < 1 .or. configuration > fks_configs) then
      call fail_metadata('FKS configuration is out of range')
    end if
  end subroutine check_configuration


  subroutine require_enabled()
    if (.not. initialized) call initialize_nlo_decay_metadata()
    if (.not. enabled) then
      call fail_metadata('nlo_decay_info.dat is absent')
    end if
  end subroutine require_enabled


  subroutine fail_metadata(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in nlo_decay_metadata: '//trim(message)
    stop 1
  end subroutine fail_metadata

end module nlo_decay_metadata
