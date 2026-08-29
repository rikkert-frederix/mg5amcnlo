module multiplicative_density_terms
  use process_dimensions, only: nexternal, validate_process_dimensions
  use fnlo_process_common, only: soft_counterevent, real_event
  use spin_density_matrix_results, only: spin_density_no_insertion, &
       spin_density_fast_virtual_insertion
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
    integer :: radiation_group = 1
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
    ! The generated luminosity provider is configuration-specific.  This is
    ! deliberately independent of EVENT_SLOT: an n-body Born/virtual atom
    ! can be sampled through a qg real sector while still owning the
    ! underlying gg luminosity.  Zero is reserved for terms (such as an
    ! LO-only decay spectator) which can never determine incoming flavours.
    integer :: luminosity_configuration = 0
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
    logical :: kinematic_families_coalesced = .false.
    type(block_distribution_term), allocatable :: terms(:)
  end type block_nlo_distribution

  type, public :: multiplicative_density_tuple
    integer :: distribution_count = 0
    integer :: sign = 1
    integer :: nlo_order = 0
    integer, allocatable :: term_indices(:)
    integer, allocatable :: event_slots(:)
  end type multiplicative_density_tuple

  ! Immutable mixed-radix traversal metadata for one set of block
  ! distributions.  Physical block identifiers are topological: production
  ! is block zero and every decay parent precedes its children.  Storing the
  ! blocks from deepest to shallowest in FASTEST_DISTRIBUTIONS consequently
  ! keeps ancestor momenta fixed while descendants are varied.  Within one
  ! block, terms with an identical event slot and luminosity configuration
  ! are contiguous so that returning to an already visited boost is avoided.
  type, public :: density_tuple_schedule
    integer :: distribution_count = 0
    integer :: tuple_count = 0
    integer :: maximum_term_count = 0
    logical :: initialized = .false.
    integer, allocatable :: fastest_distributions(:)
    integer, allocatable :: ordered_terms(:, :)
  end type density_tuple_schedule

  public :: initialize_block_distribution
  public :: initialize_block_distribution_term
  public :: set_density_primitive
  public :: finalize_block_distribution_term
  public :: finalize_block_distribution
  public :: coalesce_block_kinematic_families
  public :: density_cartesian_tuple_count
  public :: initialize_density_tuple_schedule
  public :: prepare_scheduled_density_tuple
  public :: decode_scheduled_density_tuple
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

    ! Generated fNLO code is commonly compiled with -fno-automatic.  In
    ! that mode this local derived type persists between calls, including
    ! its allocatable component.
    if (allocated(term%primitives)) deallocate(term%primitives)
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
        insertion_kind > spin_density_fast_virtual_insertion) then
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
    integer :: primitive, primitive_order, minimum_order, maximum_order

    call require_term_shape(term)
    primitive_order = -1
    minimum_order = 1
    maximum_order = 0
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
        primitive_order = -2
      end if
      minimum_order = min(minimum_order, &
           term%primitives(primitive)%nlo_order)
      maximum_order = max(maximum_order, &
           term%primitives(primitive)%nlo_order)
    end do
    if (term%nlo_order == -1) then
      if (minimum_order == maximum_order) then
        call fail_density_terms( &
             'a mixed-order momentum family contains only one order')
      end if
    else if (primitive_order /= term%nlo_order) then
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


  subroutine coalesce_block_kinematic_families(distribution)
    type(block_nlo_distribution), intent(inout) :: distribution
    type(block_distribution_term), allocatable :: families(:)
    integer, allocatable :: family_for_term(:), primitive_counts(:)
    type(density_primitive_descriptor) :: candidate
    integer :: term, previous, family, family_count, primitive, target
    integer :: matching_primitive
    integer :: minimum_order, maximum_order

    call require_distribution_shape(distribution)
    if (.not. distribution%finalized) then
      call fail_density_terms( &
           'an unfinished distribution cannot be coalesced')
    end if
    if (distribution%kinematic_families_coalesced) return
    allocate(family_for_term(distribution%term_count))
    allocate(primitive_counts(distribution%term_count))
    family_for_term = 0
    primitive_counts = 0
    family_count = 0
    do term = 1, distribution%term_count
      family = 0
      do previous = 1, term - 1
        if (distribution%terms(previous)%event_slot == &
            distribution%terms(term)%event_slot .and. &
            distribution%terms(previous)%luminosity_configuration == &
            distribution%terms(term)%luminosity_configuration) then
          family = family_for_term(previous)
          exit
        end if
      end do
      if (family == 0) then
        family_count = family_count + 1
        family = family_count
      end if
      family_for_term(term) = family
      primitive_counts(family) = primitive_counts(family) + &
           distribution%terms(term)%primitive_count
    end do

    allocate(families(family_count))
    do family = 1, family_count
      do term = 1, distribution%term_count
        if (family_for_term(term) /= family) cycle
        families(family)%block = distribution%block
        families(family)%event_slot = distribution%terms(term)%event_slot
        families(family)%luminosity_configuration = &
             distribution%terms(term)%luminosity_configuration
        exit
      end do
      families(family)%sign = 1
      families(family)%primitive_count = primitive_counts(family)
      allocate(families(family)%primitives(primitive_counts(family)))
      target = 0
      minimum_order = 1
      maximum_order = 0
      do term = 1, distribution%term_count
        if (family_for_term(term) /= family) cycle
        do primitive = 1, distribution%terms(term)%primitive_count
          candidate = distribution%terms(term)%primitives(primitive)
          candidate%scale_coefficients = distribution%terms(term)%sign* &
               candidate%scale_coefficients
          matching_primitive = 0
          do previous = 1, target
            if (same_density_primitive_key( &
                families(family)%primitives(previous), candidate)) then
              matching_primitive = previous
              exit
            end if
          end do
          ! The contraction is linear in every block density.  Equal
          ! provider descriptors can therefore be summed before loading the
          ! matrix.  Keep an exactly cancelling pair separate: finalized
          ! descriptors deliberately reject zero coefficient vectors, and
          ! retaining the pair is both exact and harmless.
          if (matching_primitive > 0 .and. any( &
              families(family)%primitives(matching_primitive)% &
              scale_coefficients + candidate%scale_coefficients /= &
              (0d0, 0d0))) then
            families(family)%primitives(matching_primitive)% &
                 scale_coefficients = &
                 families(family)%primitives(matching_primitive)% &
                 scale_coefficients + candidate%scale_coefficients
          else
            target = target + 1
            families(family)%primitives(target) = candidate
          end if
          minimum_order = min(minimum_order, &
               candidate%nlo_order)
          maximum_order = max(maximum_order, &
               candidate%nlo_order)
        end do
      end do
      call resize_family_primitives(families(family), target)
      if (minimum_order == maximum_order) then
        families(family)%nlo_order = minimum_order
      else
        families(family)%nlo_order = -1
      end if
      call finalize_block_distribution_term(families(family))
    end do
    call move_alloc(families, distribution%terms)
    distribution%term_count = family_count
    distribution%kinematic_families_coalesced = .true.
    deallocate(family_for_term)
    deallocate(primitive_counts)
  end subroutine coalesce_block_kinematic_families


  logical function same_density_primitive_key(left, right)
    type(density_primitive_descriptor), intent(in) :: left, right

    same_density_primitive_key = &
         left%insertion_kind == right%insertion_kind .and. &
         left%insertion_identifier == right%insertion_identifier .and. &
         left%insertion_rank == right%insertion_rank .and. &
         left%correlation_leg == right%correlation_leg .and. &
         left%nlo_order == right%nlo_order .and. &
         left%radiation_group == right%radiation_group .and. &
         left%laurent_poles_cancelled .eqv. &
         right%laurent_poles_cancelled
  end function same_density_primitive_key


  subroutine resize_family_primitives(term, primitive_count)
    type(block_distribution_term), intent(inout) :: term
    integer, intent(in) :: primitive_count
    type(density_primitive_descriptor), allocatable :: compact(:)

    if (primitive_count < 1 .or. primitive_count > term%primitive_count) then
      call fail_density_terms( &
           'a coalesced family has an invalid primitive count')
    end if
    if (primitive_count == term%primitive_count) return
    allocate(compact(primitive_count))
    compact = term%primitives(1:primitive_count)
    call move_alloc(compact, term%primitives)
    term%primitive_count = primitive_count
  end subroutine resize_family_primitives


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


  subroutine initialize_density_tuple_schedule(distributions, schedule)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    type(density_tuple_schedule), intent(inout) :: schedule
    integer :: distribution, candidate, position
    integer :: selected, selected_term
    integer :: maximum_term_count

    call validate_distributions(distributions)
    maximum_term_count = 0
    do distribution = 1, size(distributions)
      maximum_term_count = max(maximum_term_count, &
           distributions(distribution)%term_count)
    end do
    if (allocated(schedule%fastest_distributions)) then
      if (size(schedule%fastest_distributions) /= size(distributions)) &
           deallocate(schedule%fastest_distributions)
    end if
    if (.not. allocated(schedule%fastest_distributions)) &
         allocate(schedule%fastest_distributions(size(distributions)))
    if (allocated(schedule%ordered_terms)) then
      if (size(schedule%ordered_terms, 1) /= maximum_term_count .or. &
          size(schedule%ordered_terms, 2) /= size(distributions)) &
           deallocate(schedule%ordered_terms)
    end if
    if (.not. allocated(schedule%ordered_terms)) &
         allocate(schedule%ordered_terms( &
         maximum_term_count, size(distributions)))
    schedule%distribution_count = size(distributions)
    schedule%tuple_count = density_cartesian_tuple_count(distributions)
    schedule%maximum_term_count = maximum_term_count
    schedule%fastest_distributions = 0
    schedule%ordered_terms = 0
    schedule%initialized = .false.

    ! Select physical blocks in descending topological order.  This is the
    ! fast-to-slow radix order; block zero (production) is therefore slowest.
    do position = 1, size(distributions)
      selected = 0
      do candidate = 1, size(distributions)
        if (any(schedule%fastest_distributions(1:position - 1) == &
                candidate)) cycle
        if (selected == 0 .or. distributions(candidate)%block > &
            distributions(selected)%block) selected = candidate
      end do
      if (selected == 0) then
        call fail_density_terms('cannot order the density distributions')
      end if
      schedule%fastest_distributions(position) = selected
    end do

    ! Stable selection sort by the exact kinematic-family key.  The original
    ! term index is the final tie breaker, keeping diagnostics deterministic.
    do distribution = 1, size(distributions)
      do position = 1, distributions(distribution)%term_count
        selected_term = 0
        do candidate = 1, distributions(distribution)%term_count
          if (any(schedule%ordered_terms(1:position - 1, distribution) == &
                  candidate)) cycle
          if (selected_term == 0 .or. term_precedes( &
              distributions(distribution)%terms(candidate), candidate, &
              distributions(distribution)%terms(selected_term), &
              selected_term)) selected_term = candidate
        end do
        if (selected_term == 0) then
          call fail_density_terms('cannot order a block distribution')
        end if
        schedule%ordered_terms(position, distribution) = selected_term
      end do
    end do
    schedule%initialized = .true.
  end subroutine initialize_density_tuple_schedule


  subroutine prepare_scheduled_density_tuple(schedule, tuple)
    type(density_tuple_schedule), intent(in) :: schedule
    type(multiplicative_density_tuple), intent(inout) :: tuple

    call validate_schedule(schedule)
    if (allocated(tuple%term_indices)) then
      if (size(tuple%term_indices) /= schedule%distribution_count) &
           deallocate(tuple%term_indices)
    end if
    if (.not. allocated(tuple%term_indices)) &
         allocate(tuple%term_indices(schedule%distribution_count))
    if (allocated(tuple%event_slots)) then
      if (lbound(tuple%event_slots, 1) /= 0 .or. &
          ubound(tuple%event_slots, 1) /= nexternal) &
           deallocate(tuple%event_slots)
    end if
    if (.not. allocated(tuple%event_slots)) &
         allocate(tuple%event_slots(0:nexternal))
    tuple%distribution_count = schedule%distribution_count
  end subroutine prepare_scheduled_density_tuple


  subroutine decode_scheduled_density_tuple( &
       distributions, schedule, tuple_index, tuple)
    type(block_nlo_distribution), intent(in) :: distributions(:)
    type(density_tuple_schedule), intent(in) :: schedule
    integer, intent(in) :: tuple_index
    type(multiplicative_density_tuple), intent(inout) :: tuple
    integer :: radix, distribution, ordered_term, term, remainder

    call validate_distributions(distributions)
    call validate_schedule(schedule)
    if (schedule%distribution_count /= size(distributions) .or. &
        schedule%tuple_count /= density_cartesian_tuple_count( &
        distributions)) then
      call fail_density_terms( &
           'a tuple schedule does not match its distributions')
    end if
    if (tuple_index < 1 .or. tuple_index > schedule%tuple_count) then
      call fail_density_terms('a density tuple index is out of range')
    end if
    call prepare_scheduled_density_tuple(schedule, tuple)
    tuple%event_slots = soft_counterevent
    tuple%term_indices = 0
    tuple%sign = 1
    tuple%nlo_order = 0
    remainder = tuple_index - 1
    do radix = 1, schedule%distribution_count
      distribution = schedule%fastest_distributions(radix)
      ordered_term = mod(remainder, &
           distributions(distribution)%term_count) + 1
      remainder = remainder/distributions(distribution)%term_count
      term = schedule%ordered_terms(ordered_term, distribution)
      tuple%term_indices(distribution) = term
      tuple%event_slots(distributions(distribution)%block) = &
           distributions(distribution)%terms(term)%event_slot
      tuple%sign = tuple%sign*distributions(distribution)%terms(term)%sign
      if (tuple%nlo_order >= 0) then
        if (distributions(distribution)%terms(term)%nlo_order < 0) then
          tuple%nlo_order = -1
        else
          tuple%nlo_order = tuple%nlo_order + &
               distributions(distribution)%terms(term)%nlo_order
        end if
      end if
    end do
  end subroutine decode_scheduled_density_tuple


  logical function term_precedes(left, left_index, right, right_index)
    type(block_distribution_term), intent(in) :: left, right
    integer, intent(in) :: left_index, right_index

    if (left%event_slot /= right%event_slot) then
      term_precedes = left%event_slot < right%event_slot
    else if (left%luminosity_configuration /= &
             right%luminosity_configuration) then
      term_precedes = left%luminosity_configuration < &
           right%luminosity_configuration
    else
      term_precedes = left_index < right_index
    end if
  end function term_precedes


  subroutine validate_schedule(schedule)
    type(density_tuple_schedule), intent(in) :: schedule

    if (.not. schedule%initialized .or. &
        schedule%distribution_count < 1 .or. schedule%tuple_count < 1 .or. &
        .not. allocated(schedule%fastest_distributions) .or. &
        .not. allocated(schedule%ordered_terms)) then
      call fail_density_terms('an uninitialized tuple schedule was used')
    end if
    if (size(schedule%fastest_distributions) /= &
        schedule%distribution_count .or. &
        size(schedule%ordered_terms, 2) /= &
        schedule%distribution_count .or. &
        size(schedule%ordered_terms, 1) /= &
        schedule%maximum_term_count) then
      call fail_density_terms('a tuple schedule has the wrong shape')
    end if
  end subroutine validate_schedule


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
