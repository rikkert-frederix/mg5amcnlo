module multiplicative_product
  use iso_fortran_env, only: int64, real64
  use ieee_arithmetic, only: ieee_is_finite
  use process_dimensions, only: nsplitorders, validate_process_dimensions
  implicit none
  private

  integer, parameter, public :: product_state_born = 1
  integer, parameter, public :: product_state_finite = 2
  integer, parameter, public :: product_state_real = 3

  integer, parameter, public :: product_stage_production = 1
  integer, parameter, public :: product_stage_nlo_decay = 2

  integer, parameter, public :: product_slot_none = 0
  integer, parameter, public :: product_slot_real = 1
  integer, parameter, public :: product_slot_soft = 2
  integer, parameter, public :: product_slot_collinear = 3
  integer, parameter, public :: product_slot_soft_collinear = 4

  integer, parameter :: label_length = 64
  integer, parameter :: word_length = 32

  ! A stage event is a self-contained dispatch record.  Generated mapping,
  ! carrier and kernel dispatchers must not consult the legacy global
  ! nfksprocess while several stages are active at the same point.  The
  ! carrier dispatcher owns all simultaneous colour/spin insertions; kernel
  ! callbacks return only their scalar kinematic factors.
  type, public :: product_stage_event
    integer :: stage_id = 0
    integer :: stage_kind = 0
    integer :: parent_pdg = 0
    integer :: parent_occurrence = 0
    integer :: corrected_node = 0
    integer :: choice_id = 0
    integer :: state = 0
    integer :: source_index = 0
    integer :: configuration_index = 0
    integer :: fks_i = 0
    integer :: fks_j = 0
    integer :: fks_ij = 0
    integer :: slot = product_slot_none
    integer :: coordinate_offset = 0
    character(len=label_length) :: label = ''
    real(real64) :: xi = 0._real64
    real(real64) :: y = 0._real64
    real(real64) :: phi = 0._real64
  end type product_stage_event

  type, public :: product_event_descriptor
    integer(int64) :: sector_id = 0_int64
    integer(int64) :: counterevent_id = 0_int64
    integer :: perturbative_order = 0
    integer :: real_order = 0
    integer :: finite_order = 0
    integer :: inclusion_sign = 1
    real(real64), allocatable :: born_coordinates(:)
    type(product_stage_event), allocatable :: stages(:)
  end type product_event_descriptor

  abstract interface
    subroutine product_radiation_mapper(stage, momenta, masses, jacobian, &
                                        pass)
      import :: product_stage_event, real64
      type(product_stage_event), intent(in) :: stage
      real(real64), intent(inout) :: momenta(0:, :)
      real(real64), intent(inout) :: masses(:)
      real(real64), intent(out) :: jacobian
      logical, intent(out) :: pass
    end subroutine product_radiation_mapper

    subroutine product_carrier_evaluator(event, momenta, masses, value, pass)
      import :: product_event_descriptor, real64
      type(product_event_descriptor), intent(in) :: event
      real(real64), intent(in) :: momenta(0:, :)
      real(real64), intent(in) :: masses(:)
      real(real64), intent(out) :: value
      logical, intent(out) :: pass
    end subroutine product_carrier_evaluator

    subroutine product_local_kernel(stage, momenta, masses, value, pass)
      import :: product_stage_event, real64
      type(product_stage_event), intent(in) :: stage
      real(real64), intent(in) :: momenta(0:, :)
      real(real64), intent(in) :: masses(:)
      real(real64), intent(out) :: value
      logical, intent(out) :: pass
    end subroutine product_local_kernel
  end interface

  logical, save :: initialized = .false.
  logical, save :: enabled = .false.
  integer, save :: number_of_stages = 0
  integer(int64), save :: number_of_sectors = 0_int64
  integer(int64), save :: number_of_first_order_sectors = 0_int64
  integer, save :: maximum_radiations = 0
  integer, save :: maximum_choices = 0
  integer, save :: maximum_virtual_orders = 0

  character(len=label_length), allocatable, save :: stage_labels(:)
  integer, allocatable, save :: stage_kind_values(:)
  integer, allocatable, save :: stage_parent_values(:)
  integer, allocatable, save :: stage_occurrence_values(:)
  integer, allocatable, save :: stage_node_values(:)
  logical, allocatable, save :: stage_finite_values(:)
  integer, allocatable, save :: stage_real_source_counts(:)
  integer, allocatable, save :: stage_choice_counts(:)
  integer, allocatable, save :: stage_virtual_order_counts(:)

  integer, allocatable, save :: choice_state_values(:, :)
  integer, allocatable, save :: choice_source_values(:, :)
  integer, allocatable, save :: choice_configuration_values(:, :)
  integer, allocatable, save :: choice_i_values(:, :)
  integer, allocatable, save :: choice_j_values(:, :)
  integer, allocatable, save :: choice_ij_values(:, :)
  logical, allocatable, save :: choice_soft_values(:, :)
  logical, allocatable, save :: choice_collinear_values(:, :)
  integer, allocatable, save :: virtual_order_values(:, :, :)

  public :: initialize_multiplicative_product
  public :: has_multiplicative_product
  public :: multiplicative_product_stage_count
  public :: multiplicative_product_sector_count
  public :: multiplicative_product_first_order_sector_count
  public :: multiplicative_product_max_radiations
  public :: multiplicative_product_stage
  public :: multiplicative_product_choice
  public :: multiplicative_product_virtual_order
  public :: decode_multiplicative_sector
  public :: multiplicative_counterevent_count
  public :: decode_multiplicative_counterevent
  public :: multiplicative_phase_space_dimension
  public :: build_multiplicative_event
  public :: compose_product_radiation_maps
  public :: evaluate_product_local_kernels
  public :: evaluate_product_counterevent
  public :: product_radiation_mapper
  public :: product_carrier_evaluator
  public :: product_local_kernel

contains

  subroutine initialize_multiplicative_product()
    logical :: exists, end_seen
    logical :: format_seen, prescription_seen, enumeration_seen
    logical :: counterevents_seen, stages_seen, sectors_seen
    logical :: first_order_seen, max_radiations_seen
    logical, allocatable :: virtual_seen(:, :)
    integer :: unit_number, ios, metadata_format
    integer :: stage_records, choice_records, virtual_records
    integer :: stage, choice, virtual_index
    integer :: parent, occurrence, corrected_node, has_finite
    integer :: real_sources, choices, virtual_orders
    integer :: source, configuration, local_i, local_j, local_ij
    integer :: soft_limit, collinear_limit
    integer :: maximum_choices_found, maximum_virtual_found
    integer, allocatable :: order_buffer(:)
    integer(int64) :: sectors, first_order_sectors
    character(len=512) :: line
    character(len=label_length) :: label
    character(len=word_length) :: keyword, word, kind, state

    if (initialized) return
    call validate_process_dimensions()
    inquire(file='multiplicative_product_info.dat', exist=exists)
    if (.not. exists) then
      initialized = .true.
      enabled = .false.
      return
    end if

    open(newunit=unit_number, file='multiplicative_product_info.dat', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) call fail_product('cannot open product metadata')

    metadata_format = 0
    number_of_stages = 0
    sectors = 0_int64
    first_order_sectors = 0_int64
    maximum_radiations = -1
    stage_records = 0
    choice_records = 0
    virtual_records = 0
    maximum_choices_found = 0
    maximum_virtual_found = 0
    end_seen = .false.
    format_seen = .false.
    prescription_seen = .false.
    enumeration_seen = .false.
    counterevents_seen = .false.
    stages_seen = .false.
    sectors_seen = .false.
    first_order_seen = .false.
    max_radiations_seen = .false.

    do
      read(unit_number, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) call fail_product('cannot read product metadata')
      if (len_trim(line) == 0) cycle
      read(line, *, iostat=ios) keyword
      if (ios /= 0) call fail_product('malformed product record')
      if (end_seen) call fail_product('record found after END')
      select case (trim(keyword))
      case ('FORMAT')
        if (format_seen) call fail_product('duplicate FORMAT record')
        read(line, *, iostat=ios) keyword, metadata_format
        format_seen = .true.
      case ('PRESCRIPTION')
        if (prescription_seen) &
             call fail_product('duplicate PRESCRIPTION record')
        read(line, *, iostat=ios) keyword, word
        if (trim(word) /= 'STAGEWISE_NLO_PRODUCT') &
             call fail_product('unknown product prescription')
        prescription_seen = .true.
      case ('ENUMERATION')
        if (enumeration_seen) &
             call fail_product('duplicate ENUMERATION record')
        read(line, *, iostat=ios) keyword, word
        if (trim(word) /= 'CARTESIAN_LAZY') &
             call fail_product('unknown sector enumeration')
        enumeration_seen = .true.
      case ('COUNTEREVENTS')
        if (counterevents_seen) &
             call fail_product('duplicate COUNTEREVENTS record')
        read(line, *, iostat=ios) keyword, word
        if (trim(word) /= 'TENSOR_PRODUCT') &
             call fail_product('unknown counterevent prescription')
        counterevents_seen = .true.
      case ('STAGES')
        if (stages_seen) call fail_product('duplicate STAGES record')
        read(line, *, iostat=ios) keyword, number_of_stages
        stages_seen = .true.
      case ('SECTORS')
        if (sectors_seen) call fail_product('duplicate SECTORS record')
        read(line, *, iostat=ios) keyword, sectors
        sectors_seen = .true.
      case ('FIRST_ORDER_SECTORS')
        if (first_order_seen) &
             call fail_product('duplicate FIRST_ORDER_SECTORS record')
        read(line, *, iostat=ios) keyword, first_order_sectors
        first_order_seen = .true.
      case ('MAX_RADIATIONS')
        if (max_radiations_seen) &
             call fail_product('duplicate MAX_RADIATIONS record')
        read(line, *, iostat=ios) keyword, maximum_radiations
        max_radiations_seen = .true.
      case ('STAGE')
        read(line, *, iostat=ios) keyword, stage, label, kind, parent, &
             occurrence, corrected_node, has_finite, real_sources, &
             choices, virtual_orders
        stage_records = stage_records + 1
        maximum_choices_found = max(maximum_choices_found, choices)
        maximum_virtual_found = max(maximum_virtual_found, virtual_orders)
      case ('CHOICE')
        choice_records = choice_records + 1
      case ('VIRTUAL_ORDER')
        virtual_records = virtual_records + 1
      case ('END')
        end_seen = .true.
      case default
        call fail_product('unknown product metadata keyword')
      end select
      if (ios /= 0) call fail_product('malformed product metadata record')
    end do

    if (.not. end_seen) call fail_product('END record is absent')
    if (.not. (format_seen .and. prescription_seen .and. &
         enumeration_seen .and. counterevents_seen .and. stages_seen .and. &
         sectors_seen .and. first_order_seen .and. max_radiations_seen)) then
      call fail_product('product metadata header is incomplete')
    end if
    if (metadata_format /= 1) call fail_product('FORMAT 1 is required')
    if (number_of_stages < 1) call fail_product('invalid stage count')
    if (stage_records /= number_of_stages) &
         call fail_product('STAGE records are incomplete')
    if (sectors < 1_int64 .or. first_order_sectors < 0_int64) &
         call fail_product('invalid sector count')
    if (maximum_radiations < 0 .or. &
        maximum_radiations > number_of_stages) then
      call fail_product('invalid maximum radiation order')
    end if
    if (maximum_choices_found < 1) &
         call fail_product('a stage has no choices')
    if (nsplitorders < 1) &
         call fail_product('process split-order metadata are unavailable')

    maximum_choices = maximum_choices_found
    maximum_virtual_orders = maximum_virtual_found
    number_of_sectors = sectors
    number_of_first_order_sectors = first_order_sectors
    allocate(stage_labels(number_of_stages))
    allocate(stage_kind_values(number_of_stages))
    allocate(stage_parent_values(number_of_stages))
    allocate(stage_occurrence_values(number_of_stages))
    allocate(stage_node_values(number_of_stages))
    allocate(stage_finite_values(number_of_stages))
    allocate(stage_real_source_counts(number_of_stages))
    allocate(stage_choice_counts(number_of_stages))
    allocate(stage_virtual_order_counts(number_of_stages))
    allocate(choice_state_values(maximum_choices, number_of_stages))
    allocate(choice_source_values(maximum_choices, number_of_stages))
    allocate(choice_configuration_values(&
         maximum_choices, number_of_stages))
    allocate(choice_i_values(maximum_choices, number_of_stages))
    allocate(choice_j_values(maximum_choices, number_of_stages))
    allocate(choice_ij_values(maximum_choices, number_of_stages))
    allocate(choice_soft_values(maximum_choices, number_of_stages))
    allocate(choice_collinear_values(maximum_choices, number_of_stages))
    allocate(virtual_order_values(nsplitorders, &
         max(1, maximum_virtual_orders), number_of_stages))
    allocate(virtual_seen(max(1, maximum_virtual_orders), number_of_stages))
    allocate(order_buffer(nsplitorders))

    stage_labels = ''
    stage_kind_values = 0
    stage_parent_values = 0
    stage_occurrence_values = 0
    stage_node_values = 0
    stage_finite_values = .false.
    stage_real_source_counts = -1
    stage_choice_counts = 0
    stage_virtual_order_counts = -1
    choice_state_values = 0
    choice_source_values = 0
    choice_configuration_values = 0
    choice_i_values = 0
    choice_j_values = 0
    choice_ij_values = 0
    choice_soft_values = .false.
    choice_collinear_values = .false.
    virtual_order_values = 0
    virtual_seen = .false.

    rewind(unit_number)
    end_seen = .false.
    do
      read(unit_number, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) call fail_product('cannot reread product metadata')
      if (len_trim(line) == 0) cycle
      read(line, *, iostat=ios) keyword
      if (ios /= 0) call fail_product('malformed product body')
      if (end_seen) call fail_product('record found after END')
      select case (trim(keyword))
      case ('FORMAT', 'PRESCRIPTION', 'ENUMERATION', 'COUNTEREVENTS', &
            'STAGES', 'SECTORS', 'FIRST_ORDER_SECTORS', 'MAX_RADIATIONS')
        continue
      case ('STAGE')
        read(line, *, iostat=ios) keyword, stage, label, kind, parent, &
             occurrence, corrected_node, has_finite, real_sources, &
             choices, virtual_orders
        if (ios /= 0) call fail_product('malformed STAGE record')
        call check_stage_index(stage)
        if (stage_choice_counts(stage) /= 0) &
             call fail_product('duplicate STAGE record')
        select case (trim(kind))
        case ('PRODUCTION')
          stage_kind_values(stage) = product_stage_production
        case ('NLO_DECAY')
          stage_kind_values(stage) = product_stage_nlo_decay
        case default
          call fail_product('unknown stage kind')
        end select
        if (len_trim(label) == 0 .or. choices < 1 .or. &
            choices > maximum_choices .or. real_sources < 0 .or. &
            virtual_orders < 0 .or. &
            virtual_orders > maximum_virtual_orders .or. &
            (has_finite /= 0 .and. has_finite /= 1)) then
          call fail_product('invalid STAGE record')
        end if
        stage_labels(stage) = trim(label)
        stage_parent_values(stage) = parent
        stage_occurrence_values(stage) = occurrence
        stage_node_values(stage) = corrected_node
        stage_finite_values(stage) = has_finite == 1
        stage_real_source_counts(stage) = real_sources
        stage_choice_counts(stage) = choices
        stage_virtual_order_counts(stage) = virtual_orders
      case ('CHOICE')
        read(line, *, iostat=ios) keyword, stage, choice, state, source, &
             configuration, local_i, local_j, local_ij, soft_limit, &
             collinear_limit
        if (ios /= 0) call fail_product('malformed CHOICE record')
        call check_stage_index(stage)
        if (stage_choice_counts(stage) == 0) &
             call fail_product('CHOICE precedes its STAGE record')
        if (choice < 1 .or. choice > stage_choice_counts(stage)) &
             call fail_product('choice index is out of range')
        if (choice_state_values(choice, stage) /= 0) &
             call fail_product('duplicate CHOICE record')
        select case (trim(state))
        case ('BORN')
          choice_state_values(choice, stage) = product_state_born
        case ('FINITE')
          choice_state_values(choice, stage) = product_state_finite
        case ('REAL')
          choice_state_values(choice, stage) = product_state_real
        case default
          call fail_product('unknown stage choice')
        end select
        if ((soft_limit /= 0 .and. soft_limit /= 1) .or. &
            (collinear_limit /= 0 .and. collinear_limit /= 1)) then
          call fail_product('invalid CHOICE limit flags')
        end if
        choice_source_values(choice, stage) = source
        choice_configuration_values(choice, stage) = configuration
        choice_i_values(choice, stage) = local_i
        choice_j_values(choice, stage) = local_j
        choice_ij_values(choice, stage) = local_ij
        choice_soft_values(choice, stage) = soft_limit == 1
        choice_collinear_values(choice, stage) = collinear_limit == 1
      case ('VIRTUAL_ORDER')
        order_buffer = 0
        read(line, *, iostat=ios) keyword, stage, virtual_index, order_buffer
        if (ios /= 0) call fail_product('malformed VIRTUAL_ORDER record')
        call check_stage_index(stage)
        if (stage_virtual_order_counts(stage) < 0) &
             call fail_product('VIRTUAL_ORDER precedes its STAGE record')
        if (virtual_index < 1 .or. &
            virtual_index > stage_virtual_order_counts(stage)) then
          call fail_product('virtual-order index is out of range')
        end if
        if (virtual_seen(virtual_index, stage)) &
             call fail_product('duplicate VIRTUAL_ORDER record')
        if (any(order_buffer < 0)) &
             call fail_product('negative virtual split order')
        virtual_order_values(:, virtual_index, stage) = order_buffer
        virtual_seen(virtual_index, stage) = .true.
      case ('END')
        end_seen = .true.
      case default
        call fail_product('unknown product body keyword')
      end select
    end do
    close(unit_number)

    if (.not. end_seen) call fail_product('END record is absent')
    if (sum(stage_choice_counts) /= choice_records) &
         call fail_product('CHOICE records are incomplete')
    if (sum(stage_virtual_order_counts) /= virtual_records) &
         call fail_product('VIRTUAL_ORDER records are incomplete')
    do stage = 1, number_of_stages
      if (stage_kind_values(stage) == 0 .or. &
          stage_real_source_counts(stage) < 0 .or. &
          stage_virtual_order_counts(stage) < 0) then
        call fail_product('STAGE records are incomplete')
      end if
      if (stage_virtual_order_counts(stage) > 0) then
        if (.not. all(virtual_seen(&
             1:stage_virtual_order_counts(stage), stage))) then
          call fail_product('VIRTUAL_ORDER records are incomplete')
        end if
      end if
    end do

    call validate_product_metadata()
    deallocate(order_buffer)
    deallocate(virtual_seen)
    initialized = .true.
    enabled = .true.
  end subroutine initialize_multiplicative_product


  subroutine validate_product_metadata()
    integer :: stage, other_stage, choice, first_real_choice
    integer :: previous_source, previous_configuration
    integer(int64) :: expected_sectors, radix
    integer(int64) :: expected_first_order
    integer :: expected_maximum_radiations

    if (stage_kind_values(1) /= product_stage_production .or. &
        count(stage_kind_values == product_stage_production) /= 1) then
      call fail_product('stage one must be the unique production stage')
    end if
    if (stage_parent_values(1) /= 0 .or. &
        stage_occurrence_values(1) /= 0 .or. stage_node_values(1) /= 0) then
      call fail_product('production stage has decay-node metadata')
    end if
    do stage = 2, number_of_stages
      if (stage_kind_values(stage) /= product_stage_nlo_decay .or. &
          stage_parent_values(stage) == 0 .or. &
          stage_occurrence_values(stage) < 1 .or. &
          stage_node_values(stage) < 1) then
        call fail_product('invalid corrected-decay stage metadata')
      end if
    end do
    do stage = 1, number_of_stages
      do other_stage = stage + 1, number_of_stages
        if (trim(stage_labels(stage)) == trim(stage_labels(other_stage))) &
             call fail_product('stage labels are not unique')
      end do

      if (any(choice_state_values(&
           1:stage_choice_counts(stage), stage) == 0)) then
        call fail_product('CHOICE records are incomplete')
      end if
      if (choice_state_values(1, stage) /= product_state_born) &
           call fail_product('the first stage choice must be BORN')
      first_real_choice = 2
      if (stage_finite_values(stage)) then
        if (stage_virtual_order_counts(stage) < 1 .or. &
            stage_choice_counts(stage) < 2 .or. &
            choice_state_values(2, stage) /= product_state_finite) then
          call fail_product('FINITE choice metadata are inconsistent')
        end if
        first_real_choice = 3
      else if (stage_virtual_order_counts(stage) /= 0) then
        call fail_product('virtual orders exist without a FINITE choice')
      end if
      if (count(choice_state_values(&
           1:stage_choice_counts(stage), stage) == product_state_born) /= 1 &
          .or. count(choice_state_values(&
           1:stage_choice_counts(stage), stage) == product_state_finite) /= &
           merge(1, 0, stage_finite_values(stage))) then
        call fail_product('BORN/FINITE stage choices are not unique')
      end if
      if (first_real_choice <= stage_choice_counts(stage)) then
        if (any(choice_state_values(&
             first_real_choice:stage_choice_counts(stage), stage) /= &
             product_state_real)) then
          call fail_product('REAL stage choices are not contiguous')
        end if
      end if
      if (count(choice_state_values(&
           1:stage_choice_counts(stage), stage) == product_state_real) < &
          stage_real_source_counts(stage)) then
        call fail_product('a real source has no FKS configuration')
      end if

      previous_source = 0
      previous_configuration = 0
      do choice = 1, stage_choice_counts(stage)
        select case (choice_state_values(choice, stage))
        case (product_state_born, product_state_finite)
          if (choice_source_values(choice, stage) /= 0 .or. &
              choice_configuration_values(choice, stage) /= 0 .or. &
              choice_i_values(choice, stage) /= 0 .or. &
              choice_j_values(choice, stage) /= 0 .or. &
              choice_ij_values(choice, stage) /= 0 .or. &
              choice_soft_values(choice, stage) .or. &
              choice_collinear_values(choice, stage)) then
            call fail_product('non-real CHOICE carries FKS metadata')
          end if
        case (product_state_real)
          if (choice_source_values(choice, stage) < 1 .or. &
              choice_source_values(choice, stage) > &
                   stage_real_source_counts(stage) .or. &
              choice_configuration_values(choice, stage) < 1 .or. &
              choice_i_values(choice, stage) < 1 .or. &
              choice_j_values(choice, stage) < 1 .or. &
              choice_ij_values(choice, stage) < 1 .or. &
              choice_i_values(choice, stage) == &
                   choice_j_values(choice, stage)) then
            call fail_product('invalid real CHOICE metadata')
          end if
          if (choice_source_values(choice, stage) == previous_source) then
            if (choice_configuration_values(choice, stage) /= &
                previous_configuration + 1) then
              call fail_product('real configurations are not contiguous')
            end if
          else
            if (choice_source_values(choice, stage) /= previous_source + 1 &
                .or. choice_configuration_values(choice, stage) /= 1) then
              call fail_product('real sources are not contiguous')
            end if
          end if
          previous_source = choice_source_values(choice, stage)
          previous_configuration = &
               choice_configuration_values(choice, stage)
        end select
      end do
      if (previous_source /= stage_real_source_counts(stage)) &
           call fail_product('real source records are incomplete')
    end do

    expected_sectors = 1_int64
    expected_first_order = 0_int64
    expected_maximum_radiations = 0
    do stage = 1, number_of_stages
      radix = int(stage_choice_counts(stage), int64)
      if (expected_sectors > huge(expected_sectors)/radix) &
           call fail_product('sector count overflows a 64-bit index')
      expected_sectors = expected_sectors*radix
      expected_first_order = expected_first_order + radix - 1_int64
      if (stage_real_source_counts(stage) > 0) &
           expected_maximum_radiations = expected_maximum_radiations + 1
    end do
    if (expected_sectors /= number_of_sectors) &
         call fail_product('SECTORS does not match the Cartesian product')
    if (expected_first_order /= number_of_first_order_sectors) then
      call fail_product('FIRST_ORDER_SECTORS is inconsistent')
    end if
    if (expected_maximum_radiations /= maximum_radiations) &
         call fail_product('MAX_RADIATIONS is inconsistent')
  end subroutine validate_product_metadata


  logical function has_multiplicative_product()
    call initialize_multiplicative_product()
    has_multiplicative_product = enabled
  end function has_multiplicative_product


  integer function multiplicative_product_stage_count()
    call initialize_multiplicative_product()
    multiplicative_product_stage_count = number_of_stages
  end function multiplicative_product_stage_count


  integer(int64) function multiplicative_product_sector_count()
    call initialize_multiplicative_product()
    multiplicative_product_sector_count = number_of_sectors
  end function multiplicative_product_sector_count


  integer(int64) function multiplicative_product_first_order_sector_count()
    call initialize_multiplicative_product()
    multiplicative_product_first_order_sector_count = &
         number_of_first_order_sectors
  end function multiplicative_product_first_order_sector_count


  integer function multiplicative_product_max_radiations()
    call initialize_multiplicative_product()
    multiplicative_product_max_radiations = maximum_radiations
  end function multiplicative_product_max_radiations


  subroutine multiplicative_product_stage(stage, label, kind, parent, &
       occurrence, corrected_node, has_finite, real_sources, choices, &
       virtual_orders)
    integer, intent(in) :: stage
    character(len=*), intent(out) :: label
    integer, intent(out) :: kind, parent, occurrence, corrected_node
    logical, intent(out) :: has_finite
    integer, intent(out) :: real_sources, choices, virtual_orders

    call require_enabled()
    call check_stage_index(stage)
    label = trim(stage_labels(stage))
    kind = stage_kind_values(stage)
    parent = stage_parent_values(stage)
    occurrence = stage_occurrence_values(stage)
    corrected_node = stage_node_values(stage)
    has_finite = stage_finite_values(stage)
    real_sources = stage_real_source_counts(stage)
    choices = stage_choice_counts(stage)
    virtual_orders = stage_virtual_order_counts(stage)
  end subroutine multiplicative_product_stage


  subroutine multiplicative_product_choice(stage, choice, state, source, &
       configuration, local_i, local_j, local_ij, soft_limit, &
       collinear_limit)
    integer, intent(in) :: stage, choice
    integer, intent(out) :: state, source, configuration
    integer, intent(out) :: local_i, local_j, local_ij
    logical, intent(out) :: soft_limit, collinear_limit

    call require_enabled()
    call check_choice_index(stage, choice)
    state = choice_state_values(choice, stage)
    source = choice_source_values(choice, stage)
    configuration = choice_configuration_values(choice, stage)
    local_i = choice_i_values(choice, stage)
    local_j = choice_j_values(choice, stage)
    local_ij = choice_ij_values(choice, stage)
    soft_limit = choice_soft_values(choice, stage)
    collinear_limit = choice_collinear_values(choice, stage)
  end subroutine multiplicative_product_choice


  subroutine multiplicative_product_virtual_order(stage, virtual_index, &
                                                   orders)
    integer, intent(in) :: stage, virtual_index
    integer, intent(out) :: orders(:)

    call require_enabled()
    call check_stage_index(stage)
    if (virtual_index < 1 .or. &
        virtual_index > stage_virtual_order_counts(stage)) then
      call fail_product('virtual-order index is out of range')
    end if
    if (size(orders) /= nsplitorders) &
         call fail_product('virtual-order output has the wrong size')
    orders = virtual_order_values(:, virtual_index, stage)
  end subroutine multiplicative_product_virtual_order


  subroutine decode_multiplicative_sector(sector_id, choices)
    integer(int64), intent(in) :: sector_id
    integer, allocatable, intent(out) :: choices(:)
    integer(int64) :: remainder, radix
    integer :: stage

    call require_enabled()
    if (sector_id < 1_int64 .or. sector_id > number_of_sectors) &
         call fail_product('sector index is out of range')
    allocate(choices(number_of_stages))
    remainder = sector_id - 1_int64
    do stage = number_of_stages, 1, -1
      radix = int(stage_choice_counts(stage), int64)
      choices(stage) = int(mod(remainder, radix)) + 1
      remainder = remainder/radix
    end do
    if (remainder /= 0_int64) &
         call fail_product('mixed-radix sector decoding failed')
  end subroutine decode_multiplicative_sector


  integer(int64) function multiplicative_counterevent_count(sector_id)
    integer(int64), intent(in) :: sector_id
    integer, allocatable :: choices(:)

    call decode_multiplicative_sector(sector_id, choices)
    multiplicative_counterevent_count = &
         counterevent_count_for_choices(choices)
    deallocate(choices)
  end function multiplicative_counterevent_count


  subroutine decode_multiplicative_counterevent(sector_id, &
       counterevent_id, slots, inclusion_sign)
    integer(int64), intent(in) :: sector_id, counterevent_id
    integer, allocatable, intent(out) :: slots(:)
    integer, intent(out) :: inclusion_sign
    integer, allocatable :: choices(:)
    integer(int64) :: count, remainder, radix
    integer :: stage, local_index

    call decode_multiplicative_sector(sector_id, choices)
    count = counterevent_count_for_choices(choices)
    if (counterevent_id < 1_int64 .or. counterevent_id > count) &
         call fail_product('counterevent index is out of range')
    allocate(slots(number_of_stages))
    slots = product_slot_none
    remainder = counterevent_id - 1_int64
    do stage = number_of_stages, 1, -1
      if (choice_state_values(choices(stage), stage) /= &
          product_state_real) cycle
      radix = int(local_counterevent_count(stage, choices(stage)), int64)
      local_index = int(mod(remainder, radix)) + 1
      remainder = remainder/radix
      slots(stage) = local_counterevent_slot(&
           stage, choices(stage), local_index)
    end do
    if (remainder /= 0_int64) &
         call fail_product('mixed-radix counterevent decoding failed')
    inclusion_sign = 1
    do stage = 1, number_of_stages
      if (slots(stage) == product_slot_soft .or. &
          slots(stage) == product_slot_collinear) then
        inclusion_sign = -inclusion_sign
      end if
    end do
    deallocate(choices)
  end subroutine decode_multiplicative_counterevent


  integer function multiplicative_phase_space_dimension(sector_id, &
                                                         born_dimension)
    integer(int64), intent(in) :: sector_id
    integer, intent(in) :: born_dimension
    integer, allocatable :: choices(:)
    integer :: stage, real_order

    if (born_dimension < 0) &
         call fail_product('negative Born phase-space dimension')
    call decode_multiplicative_sector(sector_id, choices)
    real_order = 0
    do stage = 1, number_of_stages
      if (choice_state_values(choices(stage), stage) == product_state_real) &
           real_order = real_order + 1
    end do
    multiplicative_phase_space_dimension = born_dimension + 3*real_order
    deallocate(choices)
  end function multiplicative_phase_space_dimension


  subroutine build_multiplicative_event(sector_id, counterevent_id, &
       born_dimension, coordinates, event)
    integer(int64), intent(in) :: sector_id, counterevent_id
    integer, intent(in) :: born_dimension
    real(real64), intent(in) :: coordinates(:)
    type(product_event_descriptor), intent(out) :: event
    integer, allocatable :: choices(:), slots(:)
    integer :: inclusion_sign, stage, offset, state

    if (born_dimension < 0) &
         call fail_product('negative Born phase-space dimension')
    if (size(coordinates) /= &
        multiplicative_phase_space_dimension(sector_id, born_dimension)) &
         call fail_product('product coordinate vector has the wrong size')

    call decode_multiplicative_sector(sector_id, choices)
    call decode_multiplicative_counterevent(&
         sector_id, counterevent_id, slots, inclusion_sign)
    event%sector_id = sector_id
    event%counterevent_id = counterevent_id
    event%inclusion_sign = inclusion_sign
    allocate(event%born_coordinates(born_dimension))
    allocate(event%stages(number_of_stages))
    if (born_dimension > 0) &
         event%born_coordinates = coordinates(1:born_dimension)

    offset = born_dimension + 1
    do stage = 1, number_of_stages
      state = choice_state_values(choices(stage), stage)
      event%stages(stage)%stage_id = stage
      event%stages(stage)%stage_kind = stage_kind_values(stage)
      event%stages(stage)%parent_pdg = stage_parent_values(stage)
      event%stages(stage)%parent_occurrence = &
           stage_occurrence_values(stage)
      event%stages(stage)%corrected_node = stage_node_values(stage)
      event%stages(stage)%choice_id = choices(stage)
      event%stages(stage)%state = state
      event%stages(stage)%source_index = &
           choice_source_values(choices(stage), stage)
      event%stages(stage)%configuration_index = &
           choice_configuration_values(choices(stage), stage)
      event%stages(stage)%fks_i = choice_i_values(choices(stage), stage)
      event%stages(stage)%fks_j = choice_j_values(choices(stage), stage)
      event%stages(stage)%fks_ij = choice_ij_values(choices(stage), stage)
      event%stages(stage)%slot = slots(stage)
      event%stages(stage)%label = trim(stage_labels(stage))
      if (state == product_state_born) then
        continue
      else
        event%perturbative_order = event%perturbative_order + 1
      end if
      if (state == product_state_finite) then
        event%finite_order = event%finite_order + 1
      else if (state == product_state_real) then
        event%real_order = event%real_order + 1
        event%stages(stage)%coordinate_offset = offset
        event%stages(stage)%xi = coordinates(offset)
        event%stages(stage)%y = coordinates(offset + 1)
        event%stages(stage)%phi = coordinates(offset + 2)
        if (slots(stage) == product_slot_soft .or. &
            slots(stage) == product_slot_soft_collinear) &
             event%stages(stage)%xi = 0._real64
        if (slots(stage) == product_slot_collinear .or. &
            slots(stage) == product_slot_soft_collinear) &
             event%stages(stage)%y = 1._real64
        offset = offset + 3
      end if
    end do
    if (offset /= size(coordinates) + 1) &
         call fail_product('radiation-coordinate accounting failed')
    deallocate(choices)
    deallocate(slots)
  end subroutine build_multiplicative_event


  subroutine compose_product_radiation_maps(event, momenta, masses, &
                                            jacobian, mapper, pass)
    type(product_event_descriptor), intent(in) :: event
    real(real64), intent(inout) :: momenta(0:, :)
    real(real64), intent(inout) :: masses(:)
    real(real64), intent(out) :: jacobian
    procedure(product_radiation_mapper) :: mapper
    logical, intent(out) :: pass
    real(real64), allocatable :: trial_momenta(:, :), trial_masses(:)
    real(real64) :: local_jacobian
    integer :: stage
    logical :: local_pass

    call validate_event_descriptor(event)
    call validate_momentum_shapes(momenta, masses)
    allocate(trial_momenta(0:3, size(momenta, 2)))
    allocate(trial_masses(size(masses)))
    trial_momenta = momenta
    trial_masses = masses
    jacobian = 1._real64
    pass = .false.

    ! Stage IDs are topology ordered: production first, followed by the
    ! corrected root decays.  Each dispatcher sees all earlier recoil maps.
    ! Recursive nested-node ordering is intentionally a later metadata format.
    do stage = 1, size(event%stages)
      if (event%stages(stage)%state /= product_state_real) cycle
      local_jacobian = 0._real64
      local_pass = .false.
      call mapper(event%stages(stage), trial_momenta, trial_masses, &
                  local_jacobian, local_pass)
      if (.not. local_pass) then
        jacobian = 0._real64
        return
      end if
      if (.not. ieee_is_finite(local_jacobian) .or. &
          local_jacobian < 0._real64 .or. &
          .not. all(ieee_is_finite(trial_momenta)) .or. &
          .not. all(ieee_is_finite(trial_masses))) then
        jacobian = 0._real64
        return
      end if
      jacobian = jacobian*local_jacobian
      if (.not. ieee_is_finite(jacobian)) then
        jacobian = 0._real64
        return
      end if
    end do

    momenta = trial_momenta
    masses = trial_masses
    pass = .true.
  end subroutine compose_product_radiation_maps


  subroutine evaluate_product_local_kernels(event, momenta, masses, kernel, &
       stage_kernels, kernel_product, pass)
    type(product_event_descriptor), intent(in) :: event
    real(real64), intent(in) :: momenta(0:, :), masses(:)
    procedure(product_local_kernel) :: kernel
    real(real64), intent(out) :: stage_kernels(:)
    real(real64), intent(out) :: kernel_product
    logical, intent(out) :: pass
    real(real64) :: value
    integer :: stage
    logical :: local_pass

    call validate_event_descriptor(event)
    call validate_momentum_shapes(momenta, masses)
    if (size(stage_kernels) /= number_of_stages) &
         call fail_product('stage-kernel output has the wrong size')
    stage_kernels = 1._real64
    kernel_product = 1._real64
    pass = .false.
    do stage = 1, size(event%stages)
      if (event%stages(stage)%state /= product_state_real) cycle
      value = 0._real64
      local_pass = .false.
      call kernel(event%stages(stage), momenta, masses, value, local_pass)
      if (.not. local_pass .or. .not. ieee_is_finite(value)) then
        kernel_product = 0._real64
        return
      end if
      stage_kernels(stage) = value
      kernel_product = kernel_product*value
      if (.not. ieee_is_finite(kernel_product)) then
        kernel_product = 0._real64
        return
      end if
    end do
    kernel_product = real(event%inclusion_sign, real64)*kernel_product
    pass = .true.
  end subroutine evaluate_product_local_kernels


  subroutine evaluate_product_counterevent(event, seed_momenta, seed_masses, &
       mapper, carrier, kernel, mapped_momenta, mapped_masses, &
       mapping_jacobian, carrier_value, stage_kernels, weight, pass)
    type(product_event_descriptor), intent(in) :: event
    real(real64), intent(in) :: seed_momenta(0:, :), seed_masses(:)
    procedure(product_radiation_mapper) :: mapper
    procedure(product_carrier_evaluator) :: carrier
    procedure(product_local_kernel) :: kernel
    real(real64), intent(out) :: mapped_momenta(0:, :), mapped_masses(:)
    real(real64), intent(out) :: mapping_jacobian, carrier_value
    real(real64), intent(out) :: stage_kernels(:), weight
    logical, intent(out) :: pass
    real(real64) :: kernel_product
    logical :: local_pass

    call validate_event_descriptor(event)
    call validate_momentum_shapes(seed_momenta, seed_masses)
    call validate_momentum_shapes(mapped_momenta, mapped_masses)
    if (size(mapped_momenta, 2) /= size(seed_momenta, 2)) &
         call fail_product('mapped momentum output has the wrong size')
    mapped_momenta = seed_momenta
    mapped_masses = seed_masses
    mapping_jacobian = 0._real64
    carrier_value = 0._real64
    stage_kernels = 0._real64
    weight = 0._real64
    pass = .false.

    call compose_product_radiation_maps(event, mapped_momenta, &
         mapped_masses, mapping_jacobian, mapper, local_pass)
    if (.not. local_pass) then
      mapped_momenta = seed_momenta
      mapped_masses = seed_masses
      return
    end if
    call carrier(event, mapped_momenta, mapped_masses, carrier_value, &
                 local_pass)
    if (.not. local_pass .or. .not. ieee_is_finite(carrier_value)) then
      mapped_momenta = seed_momenta
      mapped_masses = seed_masses
      carrier_value = 0._real64
      mapping_jacobian = 0._real64
      return
    end if
    call evaluate_product_local_kernels(event, mapped_momenta, &
         mapped_masses, kernel, stage_kernels, kernel_product, local_pass)
    if (.not. local_pass) then
      mapped_momenta = seed_momenta
      mapped_masses = seed_masses
      carrier_value = 0._real64
      mapping_jacobian = 0._real64
      return
    end if
    weight = mapping_jacobian*carrier_value*kernel_product
    if (.not. ieee_is_finite(weight)) then
      mapped_momenta = seed_momenta
      mapped_masses = seed_masses
      mapping_jacobian = 0._real64
      carrier_value = 0._real64
      stage_kernels = 0._real64
      weight = 0._real64
      return
    end if
    pass = .true.
  end subroutine evaluate_product_counterevent


  integer(int64) function counterevent_count_for_choices(choices)
    integer, intent(in) :: choices(:)
    integer(int64) :: radix
    integer :: stage

    if (size(choices) /= number_of_stages) &
         call fail_product('sector choice vector has the wrong size')
    counterevent_count_for_choices = 1_int64
    do stage = 1, number_of_stages
      call check_choice_index(stage, choices(stage))
      if (choice_state_values(choices(stage), stage) /= &
          product_state_real) cycle
      radix = int(local_counterevent_count(stage, choices(stage)), int64)
      if (counterevent_count_for_choices > &
          huge(counterevent_count_for_choices)/radix) then
        call fail_product('counterevent count overflows a 64-bit index')
      end if
      counterevent_count_for_choices = &
           counterevent_count_for_choices*radix
    end do
  end function counterevent_count_for_choices


  integer function local_counterevent_count(stage, choice)
    integer, intent(in) :: stage, choice
    call check_choice_index(stage, choice)
    if (choice_state_values(choice, stage) /= product_state_real) &
         call fail_product('counterevent basis requested for non-real choice')
    local_counterevent_count = 1
    if (choice_soft_values(choice, stage)) &
         local_counterevent_count = local_counterevent_count + 1
    if (choice_collinear_values(choice, stage)) &
         local_counterevent_count = local_counterevent_count + 1
    if (choice_soft_values(choice, stage) .and. &
        choice_collinear_values(choice, stage)) &
         local_counterevent_count = local_counterevent_count + 1
  end function local_counterevent_count


  integer function local_counterevent_slot(stage, choice, local_index)
    integer, intent(in) :: stage, choice, local_index
    integer :: next_index

    if (local_index < 1 .or. &
        local_index > local_counterevent_count(stage, choice)) &
         call fail_product('local counterevent index is out of range')
    if (local_index == 1) then
      local_counterevent_slot = product_slot_real
      return
    end if
    next_index = 2
    if (choice_soft_values(choice, stage)) then
      if (local_index == next_index) then
        local_counterevent_slot = product_slot_soft
        return
      end if
      next_index = next_index + 1
    end if
    if (choice_collinear_values(choice, stage)) then
      if (local_index == next_index) then
        local_counterevent_slot = product_slot_collinear
        return
      end if
      next_index = next_index + 1
    end if
    if (choice_soft_values(choice, stage) .and. &
        choice_collinear_values(choice, stage) .and. &
        local_index == next_index) then
      local_counterevent_slot = product_slot_soft_collinear
      return
    end if
    call fail_product('local counterevent decoding failed')
  end function local_counterevent_slot


  subroutine validate_event_descriptor(event)
    type(product_event_descriptor), intent(in) :: event
    integer, allocatable :: choices(:), slots(:)
    integer :: inclusion_sign, stage

    call require_enabled()
    if (.not. allocated(event%born_coordinates) .or. &
        .not. allocated(event%stages) .or. &
        size(event%stages) /= number_of_stages) then
      call fail_product('product event descriptor is incomplete')
    end if
    call decode_multiplicative_sector(event%sector_id, choices)
    call decode_multiplicative_counterevent(event%sector_id, &
         event%counterevent_id, slots, inclusion_sign)
    if (event%inclusion_sign /= inclusion_sign) &
         call fail_product('product event inclusion sign is inconsistent')
    do stage = 1, number_of_stages
      if (event%stages(stage)%stage_id /= stage .or. &
          event%stages(stage)%choice_id /= choices(stage) .or. &
          event%stages(stage)%state /= &
               choice_state_values(choices(stage), stage) .or. &
          event%stages(stage)%slot /= slots(stage)) then
        call fail_product('product event stage dispatch is inconsistent')
      end if
    end do
    deallocate(choices)
    deallocate(slots)
  end subroutine validate_event_descriptor


  subroutine validate_momentum_shapes(momenta, masses)
    real(real64), intent(in) :: momenta(0:, :), masses(:)
    if (size(momenta, 1) /= 4 .or. &
        size(momenta, 2) /= size(masses)) then
      call fail_product('product momentum and mass buffers disagree')
    end if
  end subroutine validate_momentum_shapes


  subroutine check_stage_index(stage)
    integer, intent(in) :: stage
    if (stage < 1 .or. stage > number_of_stages) &
         call fail_product('stage index is out of range')
  end subroutine check_stage_index


  subroutine check_choice_index(stage, choice)
    integer, intent(in) :: stage, choice
    call check_stage_index(stage)
    if (choice < 1 .or. choice > stage_choice_counts(stage)) &
         call fail_product('choice index is out of range')
  end subroutine check_choice_index


  subroutine require_enabled()
    call initialize_multiplicative_product()
    if (.not. enabled) &
         call fail_product('no multiplicative-product metadata are present')
  end subroutine require_enabled


  subroutine fail_product(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in multiplicative_product: '//trim(message)
    stop 1
  end subroutine fail_product

end module multiplicative_product
