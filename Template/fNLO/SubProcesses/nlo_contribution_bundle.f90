module nlo_contribution_bundle
  use process_dimensions, only: fks_configs, nsplitorders, &
       amp_split_size, amp_split_orders, validate_process_dimensions
  implicit none
  private

  integer, parameter :: production_contribution = 1
  integer, parameter :: nlo_decay_contribution = 2
  logical, save :: initialized = .false.
  logical, save :: enabled = .false.
  integer, save :: number_of_contributions = 1
  integer, save :: number_of_virtual_grids = 0
  integer, allocatable, save :: contribution_kind_values(:)
  integer, allocatable, save :: contribution_first_values(:)
  integer, allocatable, save :: contribution_last_values(:)
  integer, allocatable, save :: contribution_representative_values(:)
  integer, allocatable, save :: contribution_parent_values(:)
  integer, allocatable, save :: contribution_occurrence_values(:)
  integer, allocatable, save :: contribution_node_values(:)
  logical, allocatable, save :: contribution_virtual_values(:)
  integer, allocatable, save :: configuration_owner_values(:)
  integer, allocatable, save :: virtual_grid_values(:, :)

  integer :: nfksprocess
  common /c_nfksprocess/ nfksprocess

  public :: initialize_nlo_contribution_bundle
  public :: has_nlo_contribution_bundle, nlo_contribution_count
  public :: active_nlo_contribution, contribution_for_fks
  public :: local_fks_configuration, global_fks_configuration
  public :: contribution_fks_first, contribution_fks_last
  public :: contribution_representative_fks
  public :: contribution_is_nlo_decay, contribution_parent_pdg
  public :: contribution_parent_occurrence, contribution_corrected_node
  public :: active_contribution_fks_first, active_contribution_fks_last
  public :: active_contribution_is_production
  public :: active_contribution_is_nlo_decay
  public :: contribution_has_virtual, active_contribution_has_virtual
  public :: nlo_virtual_grid_count, active_virtual_grid_index
  public :: factorized_radiation_block_count
  public :: factorized_shared_dimension
  public :: factorized_integration_dimension
  public :: factorized_radiation_start
  public :: factorized_radiation_block
  public :: bundle_species_is_nlo
  public :: bundle_component_count, bundle_born_component
  public :: bundle_production_nlo_component, bundle_width_component
  public :: bundle_decay_component, bundle_nlo_component
  public :: bundle_component_label

contains

  subroutine initialize_nlo_contribution_bundle()
    logical :: exists, end_seen
    integer :: unit_number, ios, metadata_format, contribution_count
    integer :: contribution, first, last, representative, has_virtual
    integer :: parent, occurrence, corrected_node
    integer :: configuration, expected_first, virtual_grid_count
    integer :: virtual_grid, amp_position
    integer, allocatable :: virtual_orders(:)
    character(len=512) :: line
    character(len=32) :: keyword, kind

    if (initialized) return
    call validate_process_dimensions()
    inquire(file='nlo_contribution_info.dat', exist=exists)
    if (.not. exists) then
      initialized = .true.
      enabled = .false.
      return
    end if

    open(newunit=unit_number, file='nlo_contribution_info.dat', &
         status='old', action='read', iostat=ios)
    if (ios /= 0) call fail_bundle('cannot open contribution metadata')
    metadata_format = 0
    contribution_count = 0
    virtual_grid_count = -1
    do
      read(unit_number, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) call fail_bundle('cannot read contribution metadata')
      if (len_trim(line) == 0) cycle
      read(line, *, iostat=ios) keyword
      if (ios /= 0) call fail_bundle('malformed contribution metadata')
      select case (trim(keyword))
      case ('FORMAT')
        read(line, *, iostat=ios) keyword, metadata_format
      case ('COUNT')
        read(line, *, iostat=ios) keyword, contribution_count
      case ('VIRTUAL_GRIDS')
        read(line, *, iostat=ios) keyword, virtual_grid_count
      end select
      if (ios /= 0) call fail_bundle('malformed contribution header')
    end do
    if (metadata_format /= 2 .and. metadata_format /= 3) then
      call fail_bundle('FORMAT 2 or FORMAT 3 is required')
    end if
    if (contribution_count < 2 .or. contribution_count > fks_configs) then
      call fail_bundle('invalid contribution count')
    end if
    if (virtual_grid_count < 0 .or. &
        virtual_grid_count > contribution_count*amp_split_size) then
      call fail_bundle('invalid virtual-grid count')
    end if

    number_of_contributions = contribution_count
    number_of_virtual_grids = virtual_grid_count
    allocate(contribution_kind_values(number_of_contributions))
    allocate(contribution_first_values(number_of_contributions))
    allocate(contribution_last_values(number_of_contributions))
    allocate(contribution_representative_values(number_of_contributions))
    allocate(contribution_parent_values(number_of_contributions))
    allocate(contribution_occurrence_values(number_of_contributions))
    allocate(contribution_node_values(number_of_contributions))
    allocate(contribution_virtual_values(number_of_contributions))
    allocate(configuration_owner_values(fks_configs))
    allocate(virtual_grid_values(&
         amp_split_size, number_of_contributions))
    allocate(virtual_orders(nsplitorders))
    contribution_kind_values = 0
    contribution_first_values = 0
    contribution_last_values = 0
    contribution_representative_values = 0
    contribution_parent_values = 0
    contribution_occurrence_values = 0
    contribution_node_values = 0
    contribution_virtual_values = .false.
    configuration_owner_values = 0
    virtual_grid_values = 0

    rewind(unit_number)
    end_seen = .false.
    do
      read(unit_number, '(a)', iostat=ios) line
      if (ios < 0) exit
      if (ios /= 0) call fail_bundle('cannot read contribution body')
      if (len_trim(line) == 0) cycle
      if (end_seen) call fail_bundle('record found after END')
      read(line, *, iostat=ios) keyword
      if (ios /= 0) call fail_bundle('malformed contribution keyword')
      select case (trim(keyword))
      case ('FORMAT', 'COUNT', 'VIRTUAL_GRIDS')
        continue
      case ('CONTRIBUTION')
        occurrence = 0
        corrected_node = 0
        if (metadata_format == 3) then
          read(line, *, iostat=ios) keyword, contribution, kind, first, &
               last, representative, has_virtual, parent, occurrence, &
               corrected_node
        else
          read(line, *, iostat=ios) keyword, contribution, kind, first, &
               last, representative, has_virtual, parent
        end if
        if (ios /= 0) call fail_bundle('malformed CONTRIBUTION record')
        call check_contribution(contribution)
        if (contribution_kind_values(contribution) /= 0) then
          call fail_bundle('duplicate CONTRIBUTION record')
        end if
        select case (trim(kind))
        case ('PRODUCTION')
          contribution_kind_values(contribution) = production_contribution
        case ('NLO_DECAY')
          contribution_kind_values(contribution) = nlo_decay_contribution
        case default
          call fail_bundle('unknown contribution kind')
        end select
        if (first < 1 .or. last > fks_configs .or. first > last .or. &
            representative < first .or. representative > last) then
          call fail_bundle('invalid contribution FKS range')
        end if
        if (any(configuration_owner_values(first:last) /= 0)) then
          call fail_bundle('overlapping contribution FKS ranges')
        end if
        contribution_first_values(contribution) = first
        contribution_last_values(contribution) = last
        contribution_representative_values(contribution) = representative
        contribution_virtual_values(contribution) = has_virtual /= 0
        contribution_parent_values(contribution) = parent
        contribution_occurrence_values(contribution) = occurrence
        contribution_node_values(contribution) = corrected_node
        configuration_owner_values(first:last) = contribution
      case ('VIRTUAL_GRID')
        read(line, *, iostat=ios) keyword, contribution, virtual_grid, &
             virtual_orders
        if (ios /= 0) call fail_bundle('malformed VIRTUAL_GRID record')
        call check_contribution(contribution)
        if (virtual_grid < 1 .or. &
            virtual_grid > number_of_virtual_grids) then
          call fail_bundle('virtual-grid index is out of range')
        end if
        amp_position = 0
        do configuration = 1, amp_split_size
          if (all(amp_split_orders(configuration, :) == &
                  virtual_orders)) then
            amp_position = configuration
            exit
          end if
        end do
        if (amp_position == 0) then
          call fail_bundle('a virtual order has no amplitude-order slot')
        end if
        if (virtual_grid_values(amp_position, contribution) /= 0 .or. &
            any(virtual_grid_values == virtual_grid)) then
          call fail_bundle('duplicate virtual-grid mapping')
        end if
        virtual_grid_values(amp_position, contribution) = virtual_grid
      case ('END')
        end_seen = .true.
      case default
        call fail_bundle('unknown contribution metadata keyword')
      end select
    end do
    close(unit_number)

    if (.not. end_seen) call fail_bundle('END record is absent')
    if (any(contribution_kind_values == 0) .or. &
        any(configuration_owner_values == 0)) then
      call fail_bundle('contribution metadata do not cover all FKS regions')
    end if
    if (count(virtual_grid_values > 0) /= number_of_virtual_grids) then
      call fail_bundle('virtual-grid metadata are incomplete')
    end if
    do virtual_grid = 1, number_of_virtual_grids
      if (.not. any(virtual_grid_values == virtual_grid)) then
        call fail_bundle('virtual-grid indices are not contiguous')
      end if
    end do
    if (contribution_kind_values(1) /= production_contribution .or. &
        count(contribution_kind_values == production_contribution) /= 1) then
      call fail_bundle('the bundle must have exactly one production member')
    end if
    expected_first = 1
    do contribution = 1, number_of_contributions
      if (contribution_first_values(contribution) /= expected_first) then
        call fail_bundle('contribution FKS ranges are not contiguous')
      end if
      if (contribution_kind_values(contribution) == &
          nlo_decay_contribution .and. &
          contribution_parent_values(contribution) == 0) then
        call fail_bundle('an NLO-decay member has no corrected parent')
      end if
      if (metadata_format == 3 .and. &
          contribution_kind_values(contribution) == &
          nlo_decay_contribution .and. &
          (contribution_occurrence_values(contribution) < 1 .or. &
           contribution_node_values(contribution) < 1)) then
        call fail_bundle('an NLO-decay member has no occurrence or node')
      end if
      if (contribution_kind_values(contribution) == &
          production_contribution .and. &
          contribution_parent_values(contribution) /= 0) then
        call fail_bundle('the production member has a decay parent')
      end if
      if (contribution_kind_values(contribution) == &
          production_contribution .and. &
          (contribution_occurrence_values(contribution) /= 0 .or. &
           contribution_node_values(contribution) /= 0)) then
        call fail_bundle('the production member has decay identifiers')
      end if
      if (contribution_virtual_values(contribution) .neqv. &
          any(virtual_grid_values(:, contribution) > 0)) then
        call fail_bundle('a contribution has inconsistent virtual grids')
      end if
      expected_first = contribution_last_values(contribution) + 1
    end do
    enabled = .true.
    initialized = .true.
  end subroutine initialize_nlo_contribution_bundle


  logical function has_nlo_contribution_bundle()
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    has_nlo_contribution_bundle = enabled
  end function has_nlo_contribution_bundle


  integer function nlo_contribution_count()
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    nlo_contribution_count = number_of_contributions
  end function nlo_contribution_count


  integer function active_nlo_contribution()
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    if (.not. enabled) then
      active_nlo_contribution = 1
    else if (nfksprocess >= 1 .and. nfksprocess <= fks_configs) then
      active_nlo_contribution = configuration_owner_values(nfksprocess)
    else
      active_nlo_contribution = 1
    end if
  end function active_nlo_contribution


  integer function contribution_for_fks(configuration)
    integer, intent(in) :: configuration
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    call check_configuration(configuration)
    if (enabled) then
      contribution_for_fks = configuration_owner_values(configuration)
    else
      contribution_for_fks = 1
    end if
  end function contribution_for_fks


  integer function local_fks_configuration(configuration)
    integer, intent(in) :: configuration
    integer :: owner
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    call check_configuration(configuration)
    if (.not. enabled) then
      local_fks_configuration = configuration
      return
    end if
    owner = configuration_owner_values(configuration)
    local_fks_configuration = configuration - &
         contribution_first_values(owner) + 1
  end function local_fks_configuration


  integer function global_fks_configuration(contribution, local)
    integer, intent(in) :: contribution, local
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    if (.not. enabled) then
      global_fks_configuration = local
      return
    end if
    call check_contribution(contribution)
    global_fks_configuration = contribution_first_values(contribution) + &
         local - 1
    if (global_fks_configuration < &
        contribution_first_values(contribution) .or. &
        global_fks_configuration > contribution_last_values(contribution)) then
      call fail_bundle('local FKS configuration is out of range')
    end if
  end function global_fks_configuration


  integer function contribution_fks_first(contribution)
    integer, intent(in) :: contribution
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    if (.not. enabled) then
      contribution_fks_first = 1
      return
    end if
    call check_contribution(contribution)
    contribution_fks_first = contribution_first_values(contribution)
  end function contribution_fks_first


  integer function contribution_fks_last(contribution)
    integer, intent(in) :: contribution
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    if (.not. enabled) then
      contribution_fks_last = fks_configs
      return
    end if
    call check_contribution(contribution)
    contribution_fks_last = contribution_last_values(contribution)
  end function contribution_fks_last


  integer function contribution_representative_fks(contribution)
    integer, intent(in) :: contribution
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    if (.not. enabled) then
      contribution_representative_fks = 1
      return
    end if
    call check_contribution(contribution)
    contribution_representative_fks = &
         contribution_representative_values(contribution)
  end function contribution_representative_fks


  logical function contribution_is_nlo_decay(contribution)
    integer, intent(in) :: contribution
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    call check_contribution(contribution)
    contribution_is_nlo_decay = enabled .and. &
         contribution_kind_values(contribution) == nlo_decay_contribution
  end function contribution_is_nlo_decay


  integer function contribution_parent_pdg(contribution)
    integer, intent(in) :: contribution
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    call check_contribution(contribution)
    if (.not. enabled .or. &
        contribution_kind_values(contribution) /= nlo_decay_contribution) then
      contribution_parent_pdg = 0
    else
      contribution_parent_pdg = contribution_parent_values(contribution)
    end if
  end function contribution_parent_pdg


  integer function contribution_parent_occurrence(contribution)
    integer, intent(in) :: contribution
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    call check_contribution(contribution)
    if (.not. enabled .or. &
        contribution_kind_values(contribution) /= nlo_decay_contribution) then
      contribution_parent_occurrence = 0
    else
      contribution_parent_occurrence = &
           contribution_occurrence_values(contribution)
    end if
  end function contribution_parent_occurrence


  integer function contribution_corrected_node(contribution)
    integer, intent(in) :: contribution
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    call check_contribution(contribution)
    if (.not. enabled .or. &
        contribution_kind_values(contribution) /= nlo_decay_contribution) then
      contribution_corrected_node = 0
    else
      contribution_corrected_node = contribution_node_values(contribution)
    end if
  end function contribution_corrected_node


  integer function active_contribution_fks_first()
    active_contribution_fks_first = contribution_fks_first(&
         active_nlo_contribution())
  end function active_contribution_fks_first


  integer function active_contribution_fks_last()
    active_contribution_fks_last = contribution_fks_last(&
         active_nlo_contribution())
  end function active_contribution_fks_last


  logical function active_contribution_is_production()
    integer :: contribution
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    contribution = active_nlo_contribution()
    active_contribution_is_production = enabled .and. &
         contribution_kind_values(contribution) == production_contribution
  end function active_contribution_is_production


  logical function active_contribution_is_nlo_decay()
    integer :: contribution
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    contribution = active_nlo_contribution()
    active_contribution_is_nlo_decay = enabled .and. &
         contribution_kind_values(contribution) == nlo_decay_contribution
  end function active_contribution_is_nlo_decay


  logical function contribution_has_virtual(contribution)
    integer, intent(in) :: contribution
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    if (.not. enabled) then
      contribution_has_virtual = .true.
      return
    end if
    call check_contribution(contribution)
    contribution_has_virtual = contribution_virtual_values(contribution)
  end function contribution_has_virtual


  logical function active_contribution_has_virtual()
    active_contribution_has_virtual = contribution_has_virtual(&
         active_nlo_contribution())
  end function active_contribution_has_virtual


  integer function nlo_virtual_grid_count(local_grid_count)
    integer, intent(in) :: local_grid_count
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    if (local_grid_count < 0) then
      call fail_bundle('the local virtual-grid count is negative')
    end if
    if (enabled) then
      if (local_grid_count /= amp_split_size) then
        call fail_bundle('the local virtual-grid count changed')
      end if
      nlo_virtual_grid_count = number_of_virtual_grids
    else
      nlo_virtual_grid_count = local_grid_count
    end if
  end function nlo_virtual_grid_count


  integer function active_virtual_grid_index(local_grid, local_grid_count)
    integer, intent(in) :: local_grid, local_grid_count
    integer :: contribution

    if (.not. initialized) call initialize_nlo_contribution_bundle()
    if (local_grid < 1 .or. local_grid > local_grid_count) then
      call fail_bundle('the local virtual-grid index is out of range')
    end if
    if (.not. enabled) then
      active_virtual_grid_index = local_grid
      return
    end if
    if (local_grid_count /= amp_split_size) then
      call fail_bundle('the local virtual-grid count changed')
    end if
    contribution = active_nlo_contribution()
    if (.not. contribution_virtual_values(contribution)) then
      call fail_bundle('a real-only contribution requested a virtual grid')
    end if
    active_virtual_grid_index = &
         virtual_grid_values(local_grid, contribution)
  end function active_virtual_grid_index


  integer function factorized_radiation_block_count()
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    if (enabled) then
      factorized_radiation_block_count = number_of_contributions
    else
      factorized_radiation_block_count = 1
    end if
  end function factorized_radiation_block_count


  integer function factorized_shared_dimension(integration_dimension)
    integer, intent(in) :: integration_dimension
    integer :: radiation_dimensions

    radiation_dimensions = 3*factorized_radiation_block_count()
    factorized_shared_dimension = integration_dimension - &
         radiation_dimensions
    if (factorized_shared_dimension < 0) then
      call fail_bundle('the integration dimension has no radiative block')
    end if
  end function factorized_shared_dimension


  integer function factorized_integration_dimension(canonical_dimension)
    integer, intent(in) :: canonical_dimension

    if (canonical_dimension < 3) then
      call fail_bundle('the canonical phase space has no FKS variables')
    end if
    factorized_integration_dimension = canonical_dimension + &
         3*(factorized_radiation_block_count() - 1)
  end function factorized_integration_dimension


  integer function factorized_radiation_start(canonical_dimension, block)
    integer, intent(in) :: canonical_dimension, block
    integer :: block_count

    block_count = factorized_radiation_block_count()
    if (block < 1 .or. block > block_count) then
      call fail_bundle('the radiative phase-space block is out of range')
    end if
    if (canonical_dimension < 3) then
      call fail_bundle('the canonical phase space has no FKS variables')
    end if
    factorized_radiation_start = canonical_dimension - 2 + 3*(block - 1)
  end function factorized_radiation_start


  integer function factorized_radiation_block(integration_dimension, &
                                                dimension)
    integer, intent(in) :: integration_dimension, dimension
    integer :: shared_dimension, block_count

    if (dimension < 1 .or. dimension > integration_dimension) then
      call fail_bundle('an integration-grid dimension is out of range')
    end if
    shared_dimension = factorized_shared_dimension(integration_dimension)
    if (dimension <= shared_dimension) then
      factorized_radiation_block = 0
      return
    end if
    factorized_radiation_block = &
         (dimension - shared_dimension - 1)/3 + 1
    block_count = factorized_radiation_block_count()
    if (factorized_radiation_block > block_count) then
      call fail_bundle('an integration grid has no phase-space block')
    end if
  end function factorized_radiation_block


  logical function bundle_species_is_nlo(pdg)
    integer, intent(in) :: pdg
    integer :: contribution
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    bundle_species_is_nlo = .false.
    if (.not. enabled) return
    do contribution = 1, number_of_contributions
      if (contribution_kind_values(contribution) /= &
          nlo_decay_contribution) cycle
      if (abs(contribution_parent_values(contribution)) == abs(pdg)) then
        bundle_species_is_nlo = .true.
        return
      end if
    end do
  end function bundle_species_is_nlo


  integer function bundle_component_count()
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    if (enabled) then
      ! Born, production NLO, one component per decay correction, and the
      ! strict width-expansion counterterm.
      bundle_component_count = number_of_contributions + 2
    else
      bundle_component_count = 0
    end if
  end function bundle_component_count


  integer function bundle_born_component()
    bundle_born_component = 1
  end function bundle_born_component


  integer function bundle_production_nlo_component()
    bundle_production_nlo_component = 2
  end function bundle_production_nlo_component


  integer function bundle_decay_component(contribution)
    integer, intent(in) :: contribution
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    call check_contribution(contribution)
    if (.not. enabled .or. &
        contribution_kind_values(contribution) /= nlo_decay_contribution) then
      call fail_bundle('a decay component requires an NLO-decay member')
    end if
    bundle_decay_component = contribution + 1
  end function bundle_decay_component


  integer function bundle_nlo_component(contribution)
    integer, intent(in) :: contribution

    if (.not. initialized) call initialize_nlo_contribution_bundle()
    call check_contribution(contribution)
    if (.not. enabled) then
      call fail_bundle('an NLO component requires a contribution bundle')
    end if
    ! Component one is the common Born term.  The NLO numerator of each
    ! factorized production/decay block follows in contribution order.
    bundle_nlo_component = contribution + 1
  end function bundle_nlo_component


  integer function bundle_width_component()
    if (.not. initialized) call initialize_nlo_contribution_bundle()
    if (.not. enabled) then
      bundle_width_component = 0
    else
      bundle_width_component = number_of_contributions + 2
    end if
  end function bundle_width_component


  subroutine bundle_component_label(component, label)
    integer, intent(in) :: component
    character(len=*), intent(out) :: label
    integer :: contribution

    if (.not. initialized) call initialize_nlo_contribution_bundle()
    label = ''
    if (.not. enabled) call fail_bundle('component labels require a bundle')
    if (component == bundle_born_component()) then
      label = 'BORN'
    else if (component == bundle_production_nlo_component()) then
      label = 'PRODUCTION_NLO'
    else if (component == bundle_width_component()) then
      label = 'WIDTH_COUNTERTERM'
    else
      contribution = component - 1
      call check_contribution(contribution)
      if (contribution_kind_values(contribution) /= &
          nlo_decay_contribution) then
        call fail_bundle('unknown contribution component')
      end if
      write(label, '(a,i0,a,i0,a,i0,a,i0)') &
           'DECAY_NLO_', contribution - 1, '_PDG_', &
           contribution_parent_values(contribution), '_OCC_', &
           contribution_occurrence_values(contribution), '_NODE_', &
           contribution_node_values(contribution)
    end if
  end subroutine bundle_component_label


  subroutine check_configuration(configuration)
    integer, intent(in) :: configuration
    if (configuration < 1 .or. configuration > fks_configs) then
      call fail_bundle('FKS configuration is out of range')
    end if
  end subroutine check_configuration


  subroutine check_contribution(contribution)
    integer, intent(in) :: contribution
    if (contribution < 1 .or. &
        contribution > number_of_contributions) then
      call fail_bundle('contribution index is out of range')
    end if
  end subroutine check_contribution


  subroutine fail_bundle(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in nlo_contribution_bundle: '//trim(message)
    stop 1
  end subroutine fail_bundle

end module nlo_contribution_bundle
