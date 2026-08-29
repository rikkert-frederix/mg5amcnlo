module multiplicative_process_plan
  use nlo_contribution_bundle, only: nlo_contribution_count, &
       bundle_species_is_nlo
  use decay_chain_metadata, only: decay_node_count, node_pdg
  use decay_chain_parameters, only: decay_lo_width, decay_nlo_width, &
       decay_scale_species_count
  use multiplicative_generated_metadata, only: &
       multiplicative_block_count, multiplicative_physical_blocks, &
       multiplicative_block_pdgs, multiplicative_contribution_count, &
       multiplicative_contribution_positions
  implicit none
  private

  ! Process topology and ownership never change between integration points.
  ! Keep the generated lookup results and width layout in one immutable plan
  ! instead of reconstructing them in the integrand.
  type, public :: multiplicative_process_plan_type
    logical :: initialized = .false.
    integer :: contribution_count = 0
    integer :: component_count = 0
    integer :: production_position = 0
    integer :: corrected_width_count = 0
    integer :: decay_scale_factor_count = 0
    integer, allocatable :: component_blocks(:)
    integer, allocatable :: component_pdgs(:)
    integer, allocatable :: contribution_positions(:)
    integer, allocatable :: width_blocks(:)
    double precision, allocatable :: lo_widths(:)
    double precision, allocatable :: nlo_widths(:)
  end type multiplicative_process_plan_type

  type(multiplicative_process_plan_type), target, save :: process_plan

  public :: acquire_multiplicative_process_plan

contains

  subroutine acquire_multiplicative_process_plan(plan)
    type(multiplicative_process_plan_type), pointer, intent(out) :: plan

    if (.not. process_plan%initialized) call initialize_process_plan()
    plan => process_plan
  end subroutine acquire_multiplicative_process_plan


  subroutine initialize_process_plan()
    integer :: component, contribution, node, occurrence, owner, candidate
    integer :: block, pdg

    process_plan%contribution_count = nlo_contribution_count()
    process_plan%component_count = multiplicative_block_count
    if (process_plan%contribution_count < 1 .or. &
        process_plan%component_count < process_plan%contribution_count) then
      call fail_process_plan('the generated block graph is incomplete')
    end if
    if (multiplicative_contribution_count /= &
        process_plan%contribution_count) then
      call fail_process_plan( &
           'the generated contribution metadata has the wrong size')
    end if
    allocate(process_plan%component_blocks(process_plan%component_count))
    allocate(process_plan%component_pdgs(process_plan%component_count))
    allocate(process_plan%contribution_positions( &
         process_plan%contribution_count))

    process_plan%production_position = 0
    do component = 1, process_plan%component_count
      block = multiplicative_physical_blocks(component)
      process_plan%component_blocks(component) = block
      process_plan%component_pdgs(component) = &
           multiplicative_block_pdgs(component)
      if (block == 0) then
        if (process_plan%production_position /= 0) then
          call fail_process_plan('the block graph has two production roots')
        end if
        process_plan%production_position = component
      end if
      if (component > 1) then
        if (any(process_plan%component_blocks(1:component - 1) == block)) &
             call fail_process_plan('a physical block occurs twice')
      end if
    end do
    if (process_plan%production_position == 0) then
      call fail_process_plan('the block graph has no production root')
    end if

    do contribution = 1, process_plan%contribution_count
      component = multiplicative_contribution_positions(contribution)
      if (component < 1 .or. component > process_plan%component_count) then
        call fail_process_plan('an NLO contribution has no component')
      end if
      if (contribution > 1) then
        if (any(process_plan%contribution_positions(1:contribution - 1) == &
                component)) then
          call fail_process_plan('two NLO contributions own one component')
        end if
      end if
      process_plan%contribution_positions(contribution) = component
    end do
    if (process_plan%contribution_positions(1) /= &
        process_plan%production_position) then
      call fail_process_plan('the first NLO contribution is not production')
    end if

    process_plan%corrected_width_count = 0
    do node = 1, decay_node_count()
      if (bundle_species_is_nlo(node_pdg(node))) &
           process_plan%corrected_width_count = &
           process_plan%corrected_width_count + 1
    end do
    allocate(process_plan%width_blocks( &
         process_plan%corrected_width_count))
    allocate(process_plan%lo_widths(process_plan%corrected_width_count))
    allocate(process_plan%nlo_widths(process_plan%corrected_width_count))
    occurrence = 0
    do node = 1, decay_node_count()
      pdg = node_pdg(node)
      if (.not. bundle_species_is_nlo(pdg)) cycle
      owner = 0
      do candidate = 1, process_plan%component_count
        if (process_plan%component_blocks(candidate) /= node) cycle
        owner = candidate
        exit
      end do
      if (owner == 0) then
        call fail_process_plan('a corrected width has no decay block')
      end if
      occurrence = occurrence + 1
      process_plan%width_blocks(occurrence) = owner
      process_plan%lo_widths(occurrence) = decay_lo_width(pdg)
      process_plan%nlo_widths(occurrence) = decay_nlo_width(pdg)
    end do
    process_plan%decay_scale_factor_count = decay_scale_species_count()
    process_plan%initialized = .true.
  end subroutine initialize_process_plan


  subroutine fail_process_plan(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') 'ERROR in multiplicative_process_plan: '//trim(message)
    stop 1
  end subroutine fail_process_plan

end module multiplicative_process_plan
