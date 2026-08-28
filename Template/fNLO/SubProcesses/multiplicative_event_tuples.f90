module multiplicative_event_tuples
  use process_dimensions, only: nexternal, validate_process_dimensions
  use fnlo_process_common, only: soft_counterevent, &
       collinear_counterevent, soft_collinear_counterevent, real_event
  implicit none
  private

  integer, parameter, public :: maximum_subtraction_terms = 4

  ! A family is kept as a list, not as a pre-summed block result.  Its terms
  ! generally have different momenta and therefore different cuts and
  ! observables.  Products of families are formed by the Cartesian decoder
  ! below before any signed terms are added.
  type, public :: block_subtraction_family
    integer :: block = -1
    integer :: term_count = 0
    integer :: event_slots(maximum_subtraction_terms) = soft_counterevent
    integer :: signs(maximum_subtraction_terms) = 0
  end type block_subtraction_family

  type, public :: block_event_tuple
    integer :: family_count = 0
    integer :: sign = 1
    integer, allocatable :: event_slots(:)
  end type block_event_tuple

  public :: initialize_nbody_family
  public :: initialize_real_subtraction_family
  public :: cartesian_subtraction_tuple_count
  public :: decode_cartesian_subtraction_tuple

contains

  subroutine initialize_nbody_family(family, block)
    type(block_subtraction_family), intent(out) :: family
    integer, intent(in) :: block

    call validate_block(block)
    family%block = block
    family%term_count = 1
    family%event_slots(1) = soft_counterevent
    family%signs(1) = 1
  end subroutine initialize_nbody_family


  subroutine initialize_real_subtraction_family(family, block, &
                                                 has_collinear_terms)
    type(block_subtraction_family), intent(out) :: family
    integer, intent(in) :: block
    logical, intent(in) :: has_collinear_terms

    call validate_block(block)
    family%block = block
    family%event_slots = soft_counterevent
    family%signs = 0
    if (has_collinear_terms) then
      family%term_count = 4
      family%event_slots(1:4) = (/real_event, soft_counterevent, &
           collinear_counterevent, soft_collinear_counterevent/)
      family%signs(1:4) = (/1, -1, -1, 1/)
    else
      family%term_count = 2
      family%event_slots(1:2) = (/real_event, soft_counterevent/)
      family%signs(1:2) = (/1, -1/)
    end if
  end subroutine initialize_real_subtraction_family


  integer function cartesian_subtraction_tuple_count(families)
    type(block_subtraction_family), intent(in) :: families(:)
    integer :: family

    call validate_families(families)
    cartesian_subtraction_tuple_count = 1
    do family = 1, size(families)
      if (cartesian_subtraction_tuple_count > &
          huge(cartesian_subtraction_tuple_count)/ &
          families(family)%term_count) then
        call fail_event_tuples('the Cartesian tuple count overflowed')
      end if
      cartesian_subtraction_tuple_count = &
           cartesian_subtraction_tuple_count*families(family)%term_count
    end do
  end function cartesian_subtraction_tuple_count


  subroutine decode_cartesian_subtraction_tuple(families, tuple_index, tuple)
    type(block_subtraction_family), intent(in) :: families(:)
    integer, intent(in) :: tuple_index
    type(block_event_tuple), intent(out) :: tuple
    integer :: family, term, remainder, tuple_count

    tuple_count = cartesian_subtraction_tuple_count(families)
    if (tuple_index < 1 .or. tuple_index > tuple_count) then
      call fail_event_tuples('a Cartesian tuple index is out of range')
    end if
    allocate(tuple%event_slots(0:nexternal))
    tuple%event_slots = soft_counterevent
    tuple%family_count = size(families)
    tuple%sign = 1
    remainder = tuple_index - 1
    do family = 1, size(families)
      term = mod(remainder, families(family)%term_count) + 1
      remainder = remainder/families(family)%term_count
      tuple%event_slots(families(family)%block) = &
           families(family)%event_slots(term)
      tuple%sign = tuple%sign*families(family)%signs(term)
    end do
  end subroutine decode_cartesian_subtraction_tuple


  subroutine validate_families(families)
    type(block_subtraction_family), intent(in) :: families(:)
    integer :: family, previous

    call validate_process_dimensions()
    do family = 1, size(families)
      call validate_block(families(family)%block)
      if (families(family)%term_count < 1 .or. &
          families(family)%term_count > maximum_subtraction_terms) then
        call fail_event_tuples('a subtraction family has an invalid size')
      end if
      if (any(abs(families(family)%signs( &
          1:families(family)%term_count)) /= 1)) then
        call fail_event_tuples('a subtraction-family sign is not plus/minus one')
      end if
      do previous = 1, family - 1
        if (families(previous)%block == families(family)%block) then
          call fail_event_tuples('two subtraction families own the same block')
        end if
      end do
    end do
  end subroutine validate_families


  subroutine validate_block(block)
    integer, intent(in) :: block

    call validate_process_dimensions()
    if (block < 0 .or. block > nexternal) then
      call fail_event_tuples('a factorized block index is out of range')
    end if
  end subroutine validate_block


  subroutine fail_event_tuples(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') 'ERROR in multiplicative_event_tuples: '//trim(message)
    stop 1
  end subroutine fail_event_tuples
end module multiplicative_event_tuples
