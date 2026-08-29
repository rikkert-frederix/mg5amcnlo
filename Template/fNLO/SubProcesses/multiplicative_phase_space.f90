module multiplicative_phase_space
  use process_dimensions, only: nexternal, validate_process_dimensions
  use factorized_phase_space, only: factorized_block_snapshot, &
       reset_factorized_phase_space, capture_factorized_block, &
       restore_factorized_block
  use nlo_contribution_bundle, only: nlo_contribution_count, &
       contribution_is_nlo_decay, contribution_corrected_node
  implicit none
  private

  type, public :: multiplicative_phase_space_assembly
    logical :: initialized = .false.
    type(factorized_block_snapshot), allocatable :: baseline(:)
    type(factorized_block_snapshot), allocatable :: contribution(:)
    logical, allocatable :: contribution_valid(:)
    integer, allocatable :: contribution_block(:)
  end type multiplicative_phase_space_assembly

  public :: initialize_multiplicative_phase_space_assembly
  public :: capture_multiplicative_contribution
  public :: restore_multiplicative_phase_space_assembly

contains

  subroutine initialize_multiplicative_phase_space_assembly(assembly)
    type(multiplicative_phase_space_assembly), intent(inout) :: assembly
    integer :: block, contribution_count

    call validate_process_dimensions()
    contribution_count = nlo_contribution_count()
    if (allocated(assembly%baseline)) then
      if (lbound(assembly%baseline, 1) /= 0 .or. &
          ubound(assembly%baseline, 1) /= nexternal) &
           deallocate(assembly%baseline)
    end if
    if (.not. allocated(assembly%baseline)) &
         allocate(assembly%baseline(0:nexternal))
    if (allocated(assembly%contribution)) then
      if (size(assembly%contribution) /= contribution_count) then
        deallocate(assembly%contribution)
        deallocate(assembly%contribution_valid)
        deallocate(assembly%contribution_block)
      end if
    end if
    if (.not. allocated(assembly%contribution)) then
      allocate(assembly%contribution(contribution_count))
      allocate(assembly%contribution_valid(contribution_count))
      allocate(assembly%contribution_block(contribution_count))
    end if
    assembly%initialized = .false.
    assembly%contribution_valid = .false.
    assembly%contribution_block = -1
    do block = 0, nexternal
      call capture_factorized_block(block, assembly%baseline(block))
    end do
    assembly%initialized = .true.
  end subroutine initialize_multiplicative_phase_space_assembly


  subroutine capture_multiplicative_contribution(assembly, contribution)
    type(multiplicative_phase_space_assembly), intent(inout) :: assembly
    integer, intent(in) :: contribution
    integer :: block, previous

    call require_assembly(assembly)
    if (contribution < 1 .or. &
        contribution > size(assembly%contribution)) then
      call fail_multiplicative_phase_space( &
           'an NLO contribution index is out of range')
    end if
    block = contribution_physical_block(contribution)
    do previous = 1, size(assembly%contribution)
      if (.not. assembly%contribution_valid(previous)) cycle
      if (previous /= contribution .and. &
          assembly%contribution_block(previous) == block) then
        call fail_multiplicative_phase_space( &
             'two NLO contributions overwrite the same physical block')
      end if
    end do
    call capture_factorized_block( &
         block, assembly%contribution(contribution))
    assembly%contribution_block(contribution) = block
    assembly%contribution_valid(contribution) = .true.
  end subroutine capture_multiplicative_contribution


  subroutine restore_multiplicative_phase_space_assembly(assembly)
    type(multiplicative_phase_space_assembly), intent(in) :: assembly
    integer :: block, contribution

    call require_assembly(assembly)
    if (.not. all(assembly%contribution_valid)) then
      call fail_multiplicative_phase_space( &
           'not every NLO block has a generated phase-space family')
    end if
    call reset_factorized_phase_space()
    do block = 0, nexternal
      call restore_factorized_block(assembly%baseline(block))
    end do
    do contribution = 1, size(assembly%contribution)
      call restore_factorized_block(assembly%contribution(contribution))
    end do
  end subroutine restore_multiplicative_phase_space_assembly


  integer function contribution_physical_block(contribution)
    integer, intent(in) :: contribution

    if (contribution_is_nlo_decay(contribution)) then
      contribution_physical_block = &
           contribution_corrected_node(contribution)
    else
      contribution_physical_block = 0
    end if
  end function contribution_physical_block


  subroutine require_assembly(assembly)
    type(multiplicative_phase_space_assembly), intent(in) :: assembly

    if (.not. assembly%initialized) then
      call fail_multiplicative_phase_space( &
           'the phase-space assembly is uninitialized')
    end if
  end subroutine require_assembly


  subroutine fail_multiplicative_phase_space(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_phase_space: '//trim(message)
    stop 1
  end subroutine fail_multiplicative_phase_space
end module multiplicative_phase_space
