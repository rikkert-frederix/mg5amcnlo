module multiplicative_nlo_decay
  use spin_density_matrix_results, only: spin_density_branch_result, &
       spin_density_bornlike_branch, spin_density_real_branch, &
       initialize_spin_density_branches, add_spin_density_branch, &
       spin_density_branch_leaf_count
  use factorized_phase_space, only: factorized_branch_snapshot, &
       reset_factorized_phase_space, capture_factorized_branch_snapshot, &
       restore_factorized_branch_snapshot
  implicit none
  private

  type, public :: multiplicative_nlo_workspace
    integer :: component_count = 0
    integer :: corrected_count = 0
    integer :: weight_count = 0
    integer(kind=8) :: next_mask = 0_8
    integer(kind=8) :: leaf_count = 0_8
    integer, allocatable :: component_ids(:)
    integer, allocatable :: component_open_sizes(:)
    integer, allocatable :: contribution_positions(:)
    integer, allocatable :: branch_by_component(:)
    integer, allocatable :: real_configuration_by_component(:)
    logical, allocatable :: has_snapshot(:, :)
    type(spin_density_branch_result), allocatable :: branches(:)
    type(factorized_branch_snapshot), allocatable :: snapshots(:, :)
  end type multiplicative_nlo_workspace

  type, public :: multiplicative_component_cache
    logical :: valid = .false.
    integer :: component_position = 0
    integer :: weight_count = 0
    double precision :: vegas_weight = 0d0
    type(spin_density_branch_result) :: branch
  end type multiplicative_component_cache

  public :: initialize_multiplicative_workspace
  public :: initialize_generated_multiplicative_workspace
  public :: prepare_generated_multiplicative_workspace
  public :: set_multiplicative_weight_count
  public :: reset_multiplicative_leaf_iterator
  public :: next_multiplicative_leaf
  public :: store_multiplicative_snapshot
  public :: capture_multiplicative_snapshot
  public :: set_multiplicative_real_configuration
  public :: multiplicative_leaf_has_snapshots
  public :: restore_multiplicative_leaf
  public :: add_multiplicative_block_density
  public :: complete_multiplicative_zero_branches
  public :: contract_multiplicative_leaf
  public :: invalidate_multiplicative_component_cache
  public :: store_multiplicative_component_cache
  public :: restore_multiplicative_component_cache

  interface
    integer function sdm_branch_component_count()
    end function sdm_branch_component_count

    integer function sdm_branch_corrected_count()
    end function sdm_branch_corrected_count

    integer function sdm_branch_component_id(position)
      integer, intent(in) :: position
    end function sdm_branch_component_id

    integer function sdm_branch_component_open_size(position)
      integer, intent(in) :: position
    end function sdm_branch_component_open_size

    integer function sdm_contribution_component_position(contribution)
      integer, intent(in) :: contribution
    end function sdm_contribution_component_position

    subroutine sdm_multiplicative_contraction( &
         branches, branch_choice, weight_count, result)
      import :: spin_density_branch_result
      type(spin_density_branch_result), intent(in) :: branches(*)
      integer, intent(in) :: branch_choice(*), weight_count
      complex(kind=8), intent(out) :: result(*)
    end subroutine sdm_multiplicative_contraction
  end interface

contains

  subroutine initialize_multiplicative_workspace( &
       workspace, component_ids, component_open_sizes, &
       contribution_positions, weight_count)
    type(multiplicative_nlo_workspace), intent(out) :: workspace
    integer, intent(in) :: component_ids(:), component_open_sizes(:)
    integer, intent(in) :: contribution_positions(:), weight_count
    integer :: component, contribution

    if (size(component_ids) < 1 .or. &
        size(component_open_sizes) /= size(component_ids)) then
      call fail_multiplicative_nlo('invalid component layout')
    end if
    if (size(contribution_positions) < 1 .or. weight_count < 0) then
      call fail_multiplicative_nlo('invalid corrected-block layout')
    end if
    if (any(component_open_sizes < 1)) then
      call fail_multiplicative_nlo('a component open size is invalid')
    end if
    if (any(contribution_positions < 1) .or. &
        any(contribution_positions > size(component_ids))) then
      call fail_multiplicative_nlo( &
           'a contribution position is out of range')
    end if
    do contribution = 1, size(contribution_positions)
      if (count(contribution_positions == &
                contribution_positions(contribution)) /= 1) then
        call fail_multiplicative_nlo( &
             'two NLO contributions correct the same physical block')
      end if
    end do

    workspace%component_count = size(component_ids)
    workspace%corrected_count = size(contribution_positions)
    workspace%weight_count = weight_count
    workspace%leaf_count = spin_density_branch_leaf_count( &
         workspace%corrected_count)
    allocate(workspace%component_ids(workspace%component_count))
    allocate(workspace%component_open_sizes(workspace%component_count))
    allocate(workspace%contribution_positions(workspace%corrected_count))
    allocate(workspace%branch_by_component(workspace%component_count))
    allocate(workspace%real_configuration_by_component( &
         workspace%component_count))
    allocate(workspace%has_snapshot(0:1, workspace%component_count))
    allocate(workspace%snapshots(0:1, workspace%component_count))
    workspace%component_ids = component_ids
    workspace%component_open_sizes = component_open_sizes
    workspace%contribution_positions = contribution_positions
    workspace%branch_by_component = spin_density_bornlike_branch
    workspace%real_configuration_by_component = 0
    workspace%has_snapshot = .false.
    if (weight_count > 0) then
      allocate(workspace%branches(workspace%component_count))
      do component = 1, workspace%component_count
        call initialize_spin_density_branches( &
             workspace%branches(component), component_ids(component), &
             component_open_sizes(component), weight_count)
      end do
    end if
    call reset_multiplicative_leaf_iterator(workspace)
  end subroutine initialize_multiplicative_workspace


  subroutine initialize_generated_multiplicative_workspace(workspace)
    type(multiplicative_nlo_workspace), intent(out) :: workspace
    integer, allocatable :: component_ids(:), component_open_sizes(:)
    integer, allocatable :: contribution_positions(:)
    integer :: component, contribution, component_count, corrected_count

    component_count = sdm_branch_component_count()
    corrected_count = sdm_branch_corrected_count()
    if (component_count < 1 .or. corrected_count < 1) then
      call fail_multiplicative_nlo( &
           'the generated multiplicative layout is empty')
    end if
    allocate(component_ids(component_count))
    allocate(component_open_sizes(component_count))
    allocate(contribution_positions(corrected_count))
    do component = 1, component_count
      component_ids(component) = sdm_branch_component_id(component)
      component_open_sizes(component) = &
           sdm_branch_component_open_size(component)
    end do
    do contribution = 1, corrected_count
      contribution_positions(contribution) = &
           sdm_contribution_component_position(contribution)
    end do
    call initialize_multiplicative_workspace( &
         workspace, component_ids, component_open_sizes, &
         contribution_positions, 0)
    deallocate(component_ids, component_open_sizes, &
               contribution_positions)
  end subroutine initialize_generated_multiplicative_workspace


  subroutine prepare_generated_multiplicative_workspace(workspace)
    type(multiplicative_nlo_workspace), intent(inout) :: workspace
    integer :: component, contribution

    if (workspace%component_count == 0) then
      call initialize_generated_multiplicative_workspace(workspace)
      return
    end if
    call validate_workspace_layout(workspace)
    if (workspace%component_count /= sdm_branch_component_count() .or. &
        workspace%corrected_count /= sdm_branch_corrected_count()) then
      call fail_multiplicative_nlo( &
           'the generated multiplicative layout changed during a run')
    end if
    do component = 1, workspace%component_count
      if (workspace%component_ids(component) /= &
          sdm_branch_component_id(component) .or. &
          workspace%component_open_sizes(component) /= &
          sdm_branch_component_open_size(component)) then
        call fail_multiplicative_nlo( &
             'the generated component layout changed during a run')
      end if
    end do
    do contribution = 1, workspace%corrected_count
      if (workspace%contribution_positions(contribution) /= &
          sdm_contribution_component_position(contribution)) then
        call fail_multiplicative_nlo( &
             'the generated contribution layout changed during a run')
      end if
    end do
    workspace%weight_count = 0
    workspace%real_configuration_by_component = 0
    workspace%has_snapshot = .false.
    call reset_multiplicative_leaf_iterator(workspace)
  end subroutine prepare_generated_multiplicative_workspace


  subroutine set_multiplicative_weight_count(workspace, weight_count)
    type(multiplicative_nlo_workspace), intent(inout) :: workspace
    integer, intent(in) :: weight_count
    integer :: component

    call validate_workspace_layout(workspace)
    if (weight_count < 1) then
      call fail_multiplicative_nlo('the multiplicative weight count is zero')
    end if
    if (workspace%weight_count /= 0) then
      call fail_multiplicative_nlo( &
           'the multiplicative weights were initialized twice')
    end if
    workspace%weight_count = weight_count
    if (.not. allocated(workspace%branches)) &
         allocate(workspace%branches(workspace%component_count))
    do component = 1, workspace%component_count
      call initialize_spin_density_branches( &
           workspace%branches(component), &
           workspace%component_ids(component), &
           workspace%component_open_sizes(component), weight_count)
    end do
  end subroutine set_multiplicative_weight_count


  subroutine reset_multiplicative_leaf_iterator(workspace)
    type(multiplicative_nlo_workspace), intent(inout) :: workspace

    call validate_workspace_layout(workspace)
    workspace%next_mask = 0_8
    workspace%branch_by_component = spin_density_bornlike_branch
  end subroutine reset_multiplicative_leaf_iterator


  subroutine next_multiplicative_leaf(workspace, mask, available)
    type(multiplicative_nlo_workspace), intent(inout) :: workspace
    integer(kind=8), intent(out) :: mask
    logical, intent(out) :: available
    integer :: contribution

    call validate_workspace_layout(workspace)
    available = workspace%next_mask < workspace%leaf_count
    if (.not. available) then
      mask = -1_8
      return
    end if
    mask = workspace%next_mask
    workspace%branch_by_component = spin_density_bornlike_branch
    do contribution = 1, workspace%corrected_count
      workspace%branch_by_component( &
           workspace%contribution_positions(contribution)) = &
           merge(spin_density_real_branch, &
                 spin_density_bornlike_branch, &
                 btest(mask, contribution - 1))
    end do
    workspace%next_mask = workspace%next_mask + 1_8
  end subroutine next_multiplicative_leaf


  subroutine store_multiplicative_snapshot( &
       workspace, component_position, branch, snapshot)
    type(multiplicative_nlo_workspace), intent(inout) :: workspace
    integer, intent(in) :: component_position, branch
    type(factorized_branch_snapshot), intent(in) :: snapshot

    call validate_multiplicative_snapshot( &
         workspace, component_position, branch, snapshot)
    workspace%snapshots(branch, component_position) = snapshot
    workspace%has_snapshot(branch, component_position) = .true.
  end subroutine store_multiplicative_snapshot


  subroutine capture_multiplicative_snapshot( &
       workspace, component_position, branch, event_slot)
    type(multiplicative_nlo_workspace), intent(inout) :: workspace
    integer, intent(in) :: component_position, branch, event_slot
    call validate_component_and_branch( &
         workspace, component_position, branch)
    call capture_factorized_branch_snapshot( &
         event_slot, workspace%component_ids(component_position), &
         workspace%snapshots(branch, component_position))
    call validate_multiplicative_snapshot( &
         workspace, component_position, branch, &
         workspace%snapshots(branch, component_position))
    workspace%has_snapshot(branch, component_position) = .true.
  end subroutine capture_multiplicative_snapshot


  subroutine validate_multiplicative_snapshot( &
       workspace, component_position, branch, snapshot)
    type(multiplicative_nlo_workspace), intent(in) :: workspace
    integer, intent(in) :: component_position, branch
    type(factorized_branch_snapshot), intent(in) :: snapshot

    call validate_component_and_branch( &
         workspace, component_position, branch)
    if (snapshot%block /= workspace%component_ids(component_position)) then
      call fail_multiplicative_nlo( &
           'a branch snapshot belongs to a different physical block')
    end if
    if (branch == spin_density_real_branch .and. &
        .not. any(workspace%contribution_positions == &
                 component_position)) then
      call fail_multiplicative_nlo( &
           'an uncorrected block received a real snapshot')
    end if
    if (branch == spin_density_bornlike_branch .and. &
        .not. snapshot%has_block) then
      call fail_multiplicative_nlo( &
           'a B branch snapshot has no canonical n-body block')
    end if
    if (branch == spin_density_real_branch .and. &
        .not. snapshot%has_embedded) then
      call fail_multiplicative_nlo( &
           'an R branch snapshot has no embedded real block')
    end if
  end subroutine validate_multiplicative_snapshot


  subroutine set_multiplicative_real_configuration( &
       workspace, contribution, configuration)
    type(multiplicative_nlo_workspace), intent(inout) :: workspace
    integer, intent(in) :: contribution, configuration
    integer :: component_position

    call validate_workspace_layout(workspace)
    if (contribution < 1 .or. contribution > workspace%corrected_count) then
      call fail_multiplicative_nlo( &
           'a real-branch contribution is out of range')
    end if
    if (configuration < 1) then
      call fail_multiplicative_nlo( &
           'a real-branch FKS configuration is invalid')
    end if
    component_position = workspace%contribution_positions(contribution)
    workspace%real_configuration_by_component(component_position) = &
         configuration
  end subroutine set_multiplicative_real_configuration


  subroutine restore_multiplicative_leaf(workspace, event_slot)
    type(multiplicative_nlo_workspace), intent(in) :: workspace
    integer, intent(in) :: event_slot
    integer :: branch, component, global_component

    call validate_workspace_layout(workspace)
    do component = 1, workspace%component_count
      branch = workspace%branch_by_component(component)
      if (.not. workspace%has_snapshot(branch, component)) then
        call fail_multiplicative_nlo( &
             'a selected B/R leaf has no block snapshot')
      end if
    end do
    global_component = 0
    do component = 1, workspace%component_count
      if (workspace%component_ids(component) == 0) then
        global_component = component
        exit
      end if
    end do
    if (global_component == 0) then
      call fail_multiplicative_nlo( &
           'the production component is absent from a branch workspace')
    end if

    call reset_factorized_phase_space()
    do component = 1, workspace%component_count
      branch = workspace%branch_by_component(component)
      call restore_factorized_branch_snapshot( &
           workspace%snapshots(branch, component), event_slot, &
           restore_global=component == global_component)
    end do
  end subroutine restore_multiplicative_leaf


  logical function multiplicative_leaf_has_snapshots(workspace)
    type(multiplicative_nlo_workspace), intent(in) :: workspace
    integer :: branch, component

    call validate_workspace_layout(workspace)
    multiplicative_leaf_has_snapshots = .false.
    do component = 1, workspace%component_count
      branch = workspace%branch_by_component(component)
      if (.not. workspace%has_snapshot(branch, component)) return
    end do
    multiplicative_leaf_has_snapshots = .true.
  end function multiplicative_leaf_has_snapshots


  subroutine add_multiplicative_block_density( &
       workspace, component_position, branch, density)
    type(multiplicative_nlo_workspace), intent(inout) :: workspace
    integer, intent(in) :: component_position, branch
    complex(kind=8), intent(in) :: density(:, :, :)

    call validate_component_and_branch( &
         workspace, component_position, branch)
    if (branch == spin_density_real_branch .and. &
        .not. any(workspace%contribution_positions == &
                 component_position)) then
      call fail_multiplicative_nlo( &
           'an uncorrected block received a real density')
    end if
    call add_spin_density_branch( &
         workspace%branches(component_position), branch, density)
  end subroutine add_multiplicative_block_density


  subroutine complete_multiplicative_zero_branches(workspace)
    ! A sampled map may omit one side of a local B/R pair before it records a
    ! weight line.  This includes a rejected real map and the second solution
    ! of a massive real-emission map, which intentionally has no Born-like
    ! counterevent.  Mark every missing branch as the exact zero matrix so
    ! contractions containing it vanish before they request a snapshot.
    type(multiplicative_nlo_workspace), intent(inout) :: workspace
    complex(kind=8), allocatable :: zero_density(:, :, :)
    integer :: component, contribution, open_size

    call validate_workspace(workspace)
    do contribution = 1, workspace%corrected_count
      component = workspace%contribution_positions(contribution)
      open_size = workspace%component_open_sizes(component)
      allocate(zero_density(workspace%weight_count, open_size, open_size))
      zero_density = (0d0, 0d0)
      if (.not. workspace%branches(component)%has_bornlike) then
        call add_multiplicative_block_density( &
             workspace, component, spin_density_bornlike_branch, &
             zero_density)
      end if
      if (.not. workspace%branches(component)%has_real) then
        call add_multiplicative_block_density( &
             workspace, component, spin_density_real_branch, zero_density)
      end if
      deallocate(zero_density)
    end do
  end subroutine complete_multiplicative_zero_branches


  subroutine contract_multiplicative_leaf(workspace, result)
    type(multiplicative_nlo_workspace), intent(in) :: workspace
    complex(kind=8), intent(out) :: result(:)

    call validate_workspace(workspace)
    if (size(result) /= workspace%weight_count) then
      call fail_multiplicative_nlo( &
           'a multiplicative contraction result has the wrong size')
    end if
    call sdm_multiplicative_contraction( &
         workspace%branches, workspace%branch_by_component, &
         workspace%weight_count, result)
  end subroutine contract_multiplicative_leaf


  subroutine invalidate_multiplicative_component_cache(cache)
    type(multiplicative_component_cache), intent(inout) :: cache

    cache%valid = .false.
    cache%component_position = 0
    cache%weight_count = 0
    cache%vegas_weight = 0d0
    cache%branch%has_bornlike = .false.
    cache%branch%has_real = .false.
  end subroutine invalidate_multiplicative_component_cache


  subroutine store_multiplicative_component_cache( &
       workspace, component_position, vegas_weight, cache)
    type(multiplicative_nlo_workspace), intent(in) :: workspace
    integer, intent(in) :: component_position
    double precision, intent(in) :: vegas_weight
    type(multiplicative_component_cache), intent(inout) :: cache

    call validate_workspace(workspace)
    if (component_position < 1 .or. &
        component_position > workspace%component_count .or. &
        vegas_weight <= 0d0) then
      call fail_multiplicative_nlo( &
           'cannot cache an invalid folded production component')
    end if
    cache%branch = workspace%branches(component_position)
    cache%component_position = component_position
    cache%weight_count = workspace%weight_count
    cache%vegas_weight = vegas_weight
    cache%valid = .true.
  end subroutine store_multiplicative_component_cache


  subroutine restore_multiplicative_component_cache( &
       workspace, component_position, vegas_weight, cache)
    type(multiplicative_nlo_workspace), intent(inout) :: workspace
    integer, intent(in) :: component_position
    double precision, intent(in) :: vegas_weight
    type(multiplicative_component_cache), intent(in) :: cache
    double precision :: rescaling

    call validate_workspace(workspace)
    if (.not. cache%valid .or. &
        cache%component_position /= component_position .or. &
        cache%weight_count /= workspace%weight_count .or. &
        cache%vegas_weight <= 0d0 .or. vegas_weight <= 0d0) then
      call fail_multiplicative_nlo( &
           'a folded production component cache has the wrong context')
    end if
    rescaling = vegas_weight/cache%vegas_weight
    workspace%branches(component_position) = cache%branch
    if (workspace%branches(component_position)%has_bornlike) then
      workspace%branches(component_position)%bornlike = &
           rescaling*workspace%branches(component_position)%bornlike
    end if
    if (workspace%branches(component_position)%has_real) then
      workspace%branches(component_position)%real = &
           rescaling*workspace%branches(component_position)%real
    end if
  end subroutine restore_multiplicative_component_cache


  subroutine validate_workspace(workspace)
    type(multiplicative_nlo_workspace), intent(in) :: workspace

    call validate_workspace_layout(workspace)
    if (workspace%weight_count < 1 .or. &
        .not. allocated(workspace%branches)) then
      call fail_multiplicative_nlo( &
           'the multiplicative branch weights are uninitialized')
    end if
  end subroutine validate_workspace


  subroutine validate_workspace_layout(workspace)
    type(multiplicative_nlo_workspace), intent(in) :: workspace

    if (workspace%component_count < 1 .or. &
        workspace%corrected_count < 1 .or. workspace%weight_count < 0 .or. &
        .not. allocated(workspace%component_ids) .or. &
        .not. allocated(workspace%component_open_sizes) .or. &
        .not. allocated(workspace%contribution_positions) .or. &
        .not. allocated(workspace%branch_by_component) .or. &
        .not. allocated(workspace%real_configuration_by_component) .or. &
        .not. allocated(workspace%has_snapshot) .or. &
        .not. allocated(workspace%snapshots)) then
      call fail_multiplicative_nlo( &
           'the multiplicative workspace is uninitialized')
    end if
  end subroutine validate_workspace_layout


  subroutine validate_component_and_branch( &
       workspace, component_position, branch)
    type(multiplicative_nlo_workspace), intent(in) :: workspace
    integer, intent(in) :: component_position, branch

    call validate_workspace_layout(workspace)
    if (component_position < 1 .or. &
        component_position > workspace%component_count) then
      call fail_multiplicative_nlo('a component position is out of range')
    end if
    if (branch /= spin_density_bornlike_branch .and. &
        branch /= spin_density_real_branch) then
      call fail_multiplicative_nlo('a B/R branch is invalid')
    end if
  end subroutine validate_component_and_branch


  subroutine fail_multiplicative_nlo(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') 'ERROR in multiplicative_nlo_decay: '//trim(message)
    stop 1
  end subroutine fail_multiplicative_nlo

end module multiplicative_nlo_decay
