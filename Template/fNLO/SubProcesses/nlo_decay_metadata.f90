module nlo_decay_metadata
  use process_dimensions, only: nexternal, nincoming, fks_configs, &
                                validate_process_dimensions
  use fks_metadata, only: fks_i_d, fks_j_d, validate_fks_metadata
  implicit none
  private

  integer, parameter, public :: nlo_decay_leg_target = 1
  integer, parameter, public :: nlo_decay_node_target = 2

  logical, save :: initialized = .false.
  logical, save :: enabled = .false.
  integer, save :: metadata_format = 0
  integer, save :: corrected_parent_pdg_value = 0
  integer, save :: corrected_parent_occurrence_value = 0
  integer, save :: number_of_contexts = 0
  integer, save :: number_of_fks_maps = 0
  integer, save :: number_of_production_legs = 0
  integer, save :: born_context_id = 0

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
  integer, allocatable, save :: real_to_born_values(:, :)

  public :: initialize_nlo_decay_metadata, has_nlo_decay
  public :: corrected_parent_pdg, corrected_parent_occurrence
  public :: nlo_decay_born_context, nlo_decay_context_for_fks
  public :: nlo_decay_production_count, nlo_decay_production_pdg
  public :: nlo_decay_production_is_final
  public :: nlo_decay_production_target_kind
  public :: nlo_decay_production_target_id
  public :: nlo_decay_local_count, nlo_decay_local_pdg
  public :: nlo_decay_local_is_final, nlo_decay_local_target_kind
  public :: nlo_decay_local_target_id, nlo_decay_visible_count
  public :: nlo_decay_fks_i, nlo_decay_fks_j, nlo_decay_fks_ij
  public :: nlo_decay_real_to_born

contains

  subroutine initialize_nlo_decay_metadata()
    logical :: exists, end_seen
    integer :: unit_number, ios, production_records
    integer :: context, configuration, leg, pdg, source, local_count
    integer :: visible_count, local_i, local_j, local_ij, target
    integer :: unused_first, unused_second
    character(len=512) :: line
    character(len=32) :: keyword, kind, state, name

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
             number_of_fks_maps, unused_first, unused_second
      case ('PRODUCTION_LEG')
        production_records = production_records + 1
      end select
      if (ios /= 0) call fail_metadata('malformed metadata header')
    end do

    if (metadata_format /= 3) then
      call fail_metadata('FORMAT 3 is required; regenerate the process')
    end if
    if (corrected_parent_pdg_value == 0 .or. &
        corrected_parent_occurrence_value < 1) then
      call fail_metadata('invalid corrected-decay parent')
    end if
    if (number_of_contexts < 2 .or. &
        number_of_fks_maps /= fks_configs) then
      call fail_metadata('invalid metadata counts')
    end if
    if (production_records < nincoming + 1 .or. &
        production_records > nexternal) then
      call fail_metadata('invalid production-leg count')
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
    allocate(real_to_born_values(nexternal, fks_configs))

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
    real_to_born_values = 0

    rewind(unit_number)
    end_seen = .false.
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
            'VIRTUAL_COMPOSITION', 'VIRTUAL_CURRENT_COUNT', 'COUNTS', &
            'FKS_PARTNER', 'COLOR_LINK')
        continue
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
          unused_first = 1
        case ('J')
          unused_first = 2
        case ('IJ')
          unused_first = 3
        case default
          call fail_metadata('unknown FKS target name')
        end select
        if (fks_target_kinds(unused_first, configuration) /= 0) then
          call fail_metadata('duplicate FKS_TARGET record')
        end if
        fks_target_kinds(unused_first, configuration) = target_kind(kind)
        fks_target_ids(unused_first, configuration) = target
        fks_target_local_legs(unused_first, configuration) = leg
      case ('REAL_BORN_MAP')
        read(line, *, iostat=ios) keyword, configuration, leg, target
        call check_configuration(configuration)
        call check_local_storage_index(leg)
        if (real_to_born_values(leg, configuration) /= 0) then
          call fail_metadata('duplicate REAL_BORN_MAP record')
        end if
        real_to_born_values(leg, configuration) = target
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
  end subroutine initialize_nlo_decay_metadata


  subroutine validate_metadata_values()
    integer :: leg, context, configuration, target, born_count
    integer :: node_maps, expected_visible
    logical :: visible_used(nexternal), born_local_used(nexternal)

    if (any(production_pdg_values == 0)) then
      call fail_metadata('production-leg records are incomplete')
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
      node_maps = 0
      do leg = 1, number_of_production_legs
        target = production_target_ids(leg, context)
        select case (production_target_kinds(leg, context))
        case (nlo_decay_leg_target)
          call mark_visible_target(target, context_visible_counts(context), &
                                   visible_used)
        case (nlo_decay_node_target)
          if (.not. production_final_values(leg) .or. target /= 1 .or. &
              production_pdg_values(leg) /= corrected_parent_pdg_value) then
            call fail_metadata('invalid corrected production node target')
          end if
          node_maps = node_maps + 1
        case default
          call fail_metadata('invalid production target kind')
        end select
      end do
      if (node_maps /= 1) then
        call fail_metadata('corrected parent is not mapped exactly once')
      end if

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
          if (local_final_values(leg, context) .or. target /= 1 .or. &
              local_pdg_values(leg, context) /= corrected_parent_pdg_value) then
            call fail_metadata('invalid local corrected-parent target')
          end if
          node_maps = node_maps + 1
        case default
          call fail_metadata('invalid local target kind')
        end select
      end do
      if (node_maps /= 1 .or. &
          .not. all(visible_used(1:context_visible_counts(context)))) then
        call fail_metadata('context does not cover the visible event')
      end if
    end do
    if (born_count /= 1) then
      call fail_metadata('exactly one Born context is required')
    end if

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


  integer function nlo_decay_real_to_born(configuration, real_leg)
    integer, intent(in) :: configuration, real_leg
    call require_enabled()
    call check_configuration(configuration)
    call check_local_leg(fks_context_values(configuration), real_leg)
    nlo_decay_real_to_born = real_to_born_values(real_leg, configuration)
  end function nlo_decay_real_to_born


  subroutine check_production_leg(leg)
    integer, intent(in) :: leg
    if (leg < 1 .or. leg > number_of_production_legs) then
      call fail_metadata('production leg is out of range')
    end if
  end subroutine check_production_leg


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
