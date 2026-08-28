module multiplicative_density_terms
  use process_dimensions, only: nexternal, validate_process_dimensions
  use fnlo_process_common, only: soft_counterevent, real_event
  use spin_density_matrix_results, only: spin_density_no_insertion, &
       spin_density_color_insertion
  implicit none
  private

  integer, parameter, public :: density_scale_coefficient_count = 3

  ! A primitive is one separately generated density matrix multiplied by a
  ! block-local scalar kernel.  The three coefficients multiply respectively
  ! 1, log(mu_R^2/Q^2), and log(mu_F^2/Q^2).  Products of primitives therefore
  ! retain the full multivariate scale polynomial instead of projecting it
  ! back onto a single linear weight line.
  type, public :: density_primitive_descriptor
    integer :: insertion_kind = spin_density_no_insertion
    integer :: insertion_identifier = 0
    integer :: insertion_rank = 0
    integer :: correlation_leg = 0
    integer :: nlo_order = 0
    logical :: laurent_poles_cancelled = .false.
    complex(kind=8) :: scale_coefficients( &
         density_scale_coefficient_count) = (0d0, 0d0)
  end type density_primitive_descriptor

  ! Every term owns exactly one mapped momentum point.  Several primitives
  ! may be summed inside a term only because they share that point (for
  ! example the colour-linked pieces of one soft counterevent, or the finite
  ! B+V+I operator on one Born point).  Terms at R/S/C/SC momenta remain
  ! distinct and are never pre-added here.
  type, public :: block_distribution_term
    integer :: block = -1
    integer :: event_slot = -1
    integer :: sign = 0
    integer :: nlo_order = 0
    integer :: primitive_count = 0
    logical :: finalized = .false.
    type(density_primitive_descriptor), allocatable :: primitives(:)
  end type block_distribution_term

  type, public :: block_nlo_distribution
    integer :: block = -1
    integer :: term_count = 0
    logical :: finalized = .false.
    type(block_distribution_term), allocatable :: terms(:)
  end type block_nlo_distribution

  type, public :: multiplicative_density_tuple
    integer :: distribution_count = 0
    integer :: sign = 1
    integer :: nlo_order = 0
    integer, allocatable :: term_indices(:)
    integer, allocatable :: event_slots(:)
  end type multiplicative_density_tuple

  public :: initialize_block_distribution
  public :: initialize_block_distribution_term
  public :: set_density_primitive
  public :: finalize_block_distribution_term
  public :: finalize_block_distribution
  public :: density_cartesian_tuple_count
  public :: decode_density_cartesian_tuple
  public :: evaluate_density_primitive_coefficient

contains

  subroutine initialize_block_distribution(distribution, block, term_count)
    type(block_nlo_distribution), intent(out) :: distribution
    integer, intent(in) :: block, term_count

    call validate_block(block)
    if (term_count < 1) then
      call fail_density_terms('a block distribution has no terms')
    end if
    distribution%block = block
    distribution%term_count = term_count
    allocate(distribution%terms(term_count))
  end subroutine initialize_block_distribution


  subroutine initialize_block_distribution_term( &
       distribution, term_index, event_slot, sign, nlo_order, &
       primitive_count)
    type(block_nlo_distribution), intent(inout) :: distribution
    integer, intent(in) :: term_index, event_slot, sign, nlo_order
    integer, intent(in) :: primitive_count
    type(block_distribution_term) :: term

    call require_distribution_shape(distribution)
    if (distribution%finalized) then
      call fail_density_terms('a finalized distribution was modified')
    end if
    if (term_index < 1 .or. term_index > distribution%term_count) then
      call fail_density_terms('a block-term index is out of range')
    end if
    call validate_event_slot(event_slot)
    if (abs(sign) /= 1) then
      call fail_density_terms('a block-term sign is not plus/minus one')
    end if
    if (nlo_order < 0 .or. nlo_order > 1) then
      call fail_density_terms( &
           'one block term must be LO or one-NLO-order')
    end if
    if (primitive_count < 1) then
      call fail_density_terms('a block term has no density primitives')
    end if
    if (distribution%terms(term_index)%block /= -1) then
      call fail_density_terms('a block term was initialized twice')
    end if

    term%block = distribution%block
    term%event_slot = event_slot
    term%sign = sign
    term%nlo_order = nlo_order
    term%primitive_count = primitive_count
    allocate(term%primitives(primitive_count))
    distribution%terms(term_index) = term
  end subroutine initialize_block_distribution_term


  subroutine set_density_primitive( &
       term, primitive_index, insertion_kind, insertion_identifier, &
       insertion_rank, correlation_leg, nlo_order, scale_coefficients, &
       laurent_poles_cancelled)
    type(block_distribution_term), intent(inout) :: term
    integer, intent(in) :: primitive_index, insertion_kind
    integer, intent(in) :: insertion_identifier, insertion_rank
    integer, intent(in) :: correlation_leg, nlo_order
    complex(kind=8), intent(in) :: scale_coefficients( &
         density_scale_coefficient_count)
    logical, intent(in) :: laurent_poles_cancelled

    call require_term_shape(term)
    if (term%finalized) then
      call fail_density_terms('a finalized block term was modified')
    end if
    if (primitive_index < 1 .or. &
        primitive_index > term%primitive_count) then
      call fail_density_terms('a density-primitive index is out of range')
    end if
    if (insertion_kind < spin_density_no_insertion .or. &
        insertion_kind > spin_density_color_insertion) then
      call fail_density_terms('a density-primitive kind is invalid')
    end if
    if (insertion_identifier < 0 .or. insertion_rank < 0 .or. &
        correlation_leg < 0) then
      call fail_density_terms('a density-primitive identifier is negative')
    end if
    if (nlo_order < 0 .or. nlo_order > 1) then
      call fail_density_terms('a density primitive has invalid NLO order')
    end if
    if (insertion_kind == spin_density_no_insertion .and. &
        (insertion_identifier /= 0 .or. insertion_rank /= 0 .or. &
         correlation_leg /= 0)) then
      call fail_density_terms('an LO primitive carries insertion metadata')
    end if
    if (insertion_kind /= spin_density_no_insertion .and. &
        insertion_identifier == 0) then
      call fail_density_terms('an insertion primitive has no identifier')
    end if

    term%primitives(primitive_index)%insertion_kind = insertion_kind
    term%primitives(primitive_index)%insertion_identifier = &
         insertion_identifier
    term%primitives(primitive_index)%insertion_rank = insertion_rank
    term%primitives(primitive_index)%correlation_leg = correlation_leg
    term%primitives(primitive_index)%nlo_order = nlo_order
    term%primitives(primitive_index)%scale_coefficients = &
         scale_coefficients
    term%primitives(primitive_index)%laurent_poles_cancelled = &
         laurent_poles_cancelled
  end subroutine set_density_primitive


  subroutine finalize_block_distribution_term(term)
    type(block_distribution_term), intent(inout) :: term
    integer :: primitive, primitive_order

    call require_term_shape(term)
    primitive_order = -1
    do primitive = 1, term%primitive_count
      if (.not. term%primitives(primitive)%laurent_poles_cancelled) then
        call fail_density_terms( &
             'a raw Laurent-pole density cannot enter a product')
      end if
      if (all(term%primitives(primitive)%scale_coefficients == &
              (0d0, 0d0))) then
        call fail_density_terms('a density primitive has zero coefficient')
      end if
      if (primitive_order < 0) then
        primitive_order = term%primitives(primitive)%nlo_order
      else if (primitive_order /= &
               term%primitives(primitive)%nlo_order) then
        call fail_density_terms( &
             'one momentum term mixes different perturbative orders')
      end if
    end do
    if (primitive_order /= term%nlo_order) then
      call fail_density_terms( &
           'a block term disagrees with its density-primitive order')
    end if
    term%finalized = .true.
  end subroutine finalize_block_distribution_term


  subroutine finalize_block_distribution(distribution)
    type(block_nlo_distribution), intent(inout) :: distribution
    integer :: term

    call require_distribution_shape(distribution)
    do term = 1, distribution%term_count
      if (.not. distribution%terms(term)%finalized) then
        call fail_density_terms( &
             'a block distribution contains an unfinished term')
      end if
      if (distribution%terms(term)%block /= distribution%block) then
        call fail_density_terms('a block distribution mixes block owners')
      end if
    end do
    distribution%finalized = .true.
  end subroutine finalize_block_distribution


  integer function density_cartesian_tuple_count(distributions)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    integer :: distribution

    call validate_distributions(distributions)
    density_cartesian_tuple_count = 1
    do distribution = 1, size(distributions)
      if (density_cartesian_tuple_count > &
          huge(density_cartesian_tuple_count)/ &
          distributions(distribution)%term_count) then
        call fail_density_terms('the density tuple count overflowed')
      end if
      density_cartesian_tuple_count = density_cartesian_tuple_count* &
           distributions(distribution)%term_count
    end do
  end function density_cartesian_tuple_count


  subroutine decode_density_cartesian_tuple( &
       distributions, tuple_index, tuple)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    integer, intent(in) :: tuple_index
    type(multiplicative_density_tuple), intent(out) :: tuple
    integer :: distribution, term, remainder, tuple_count

    tuple_count = density_cartesian_tuple_count(distributions)
    if (tuple_index < 1 .or. tuple_index > tuple_count) then
      call fail_density_terms('a density tuple index is out of range')
    end if
    allocate(tuple%term_indices(size(distributions)))
    allocate(tuple%event_slots(0:nexternal))
    tuple%event_slots = soft_counterevent
    tuple%distribution_count = size(distributions)
    tuple%sign = 1
    tuple%nlo_order = 0
    remainder = tuple_index - 1
    do distribution = 1, size(distributions)
      term = mod(remainder, distributions(distribution)%term_count) + 1
      remainder = remainder/distributions(distribution)%term_count
      tuple%term_indices(distribution) = term
      tuple%event_slots(distributions(distribution)%block) = &
           distributions(distribution)%terms(term)%event_slot
      tuple%sign = tuple%sign*distributions(distribution)%terms(term)%sign
      tuple%nlo_order = tuple%nlo_order + &
           distributions(distribution)%terms(term)%nlo_order
    end do
  end subroutine decode_density_cartesian_tuple


  complex(kind=8) function evaluate_density_primitive_coefficient( &
       primitive, logarithmic_mu2_r, logarithmic_mu2_f)
    type(density_primitive_descriptor), intent(in) :: primitive
    double precision, intent(in) :: logarithmic_mu2_r, logarithmic_mu2_f

    if (.not. primitive%laurent_poles_cancelled) then
      call fail_density_terms( &
           'a raw Laurent-pole coefficient was evaluated')
    end if
    evaluate_density_primitive_coefficient = &
         primitive%scale_coefficients(1) + &
         primitive%scale_coefficients(2)*logarithmic_mu2_r + &
         primitive%scale_coefficients(3)*logarithmic_mu2_f
  end function evaluate_density_primitive_coefficient


  subroutine validate_distributions(distributions)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    integer :: distribution, previous

    call validate_process_dimensions()
    if (size(distributions) < 1) then
      call fail_density_terms('there are no block distributions')
    end if
    do distribution = 1, size(distributions)
      call require_distribution_shape(distributions(distribution))
      if (.not. distributions(distribution)%finalized) then
        call fail_density_terms('an unfinalized distribution was scheduled')
      end if
      do previous = 1, distribution - 1
        if (distributions(previous)%block == &
            distributions(distribution)%block) then
          call fail_density_terms( &
               'two distributions own the same physical block')
        end if
      end do
    end do
  end subroutine validate_distributions


  subroutine require_distribution_shape(distribution)
    type(block_nlo_distribution), intent(in) :: distribution

    call validate_block(distribution%block)
    if (distribution%term_count < 1 .or. &
        .not. allocated(distribution%terms) .or. &
        size(distribution%terms) /= distribution%term_count) then
      call fail_density_terms('a block distribution has invalid storage')
    end if
  end subroutine require_distribution_shape


  subroutine require_term_shape(term)
    type(block_distribution_term), intent(in) :: term

    call validate_block(term%block)
    call validate_event_slot(term%event_slot)
    if (term%primitive_count < 1 .or. .not. allocated(term%primitives) .or. &
        size(term%primitives) /= term%primitive_count) then
      call fail_density_terms('a block term has invalid primitive storage')
    end if
  end subroutine require_term_shape


  subroutine validate_event_slot(event_slot)
    integer, intent(in) :: event_slot

    if (event_slot < soft_counterevent .or. event_slot > real_event) then
      call fail_density_terms('a density term has an invalid event slot')
    end if
  end subroutine validate_event_slot


  subroutine validate_block(block)
    integer, intent(in) :: block

    call validate_process_dimensions()
    if (block < 0 .or. block > nexternal) then
      call fail_density_terms('a density block index is out of range')
    end if
  end subroutine validate_block


  subroutine fail_density_terms(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_density_terms: '//trim(message)
    stop 1
  end subroutine fail_density_terms
end module multiplicative_density_terms
