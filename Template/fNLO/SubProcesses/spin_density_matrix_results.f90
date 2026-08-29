module spin_density_matrix_results
  use process_dimensions, only: nexternal
  use fnlo_process_common, only: soft_counterevent, real_event
  use factorized_phase_space, only: &
       factorized_block_momentum_revision
  implicit none
  private

  integer, parameter, public :: spin_density_no_insertion = 0
  integer, parameter, public :: spin_density_born_insertion = 1
  integer, parameter, public :: spin_density_real_insertion = 2
  integer, parameter, public :: spin_density_virtual_insertion = 3
  integer, parameter, public :: spin_density_color_insertion = 4

  ! One physical production or decay block.  LO is the reusable spectator
  ! density.  INSERTION is the block-local object selected for one term of a
  ! subtraction tuple (an underlying Born, real, virtual, spin-correlated, or
  ! colour-correlated density).  Strict additive contractions select an
  ! insertion in at most one block; multiplicative contractions may select
  ! one insertion in every block, while retaining each block's own event slot.
  type, public :: spin_density_block_result
    integer :: block = -1
    integer :: event_slot = -1
    integer :: open_size = 0
    integer :: insertion_kind = spin_density_no_insertion
    integer :: insertion_order = 0
    logical :: has_lo = .false.
    logical :: has_insertion = .false.
    complex(kind=8), allocatable :: lo(:, :, :)
    complex(kind=8), allocatable :: insertion(:, :, :)
  end type spin_density_block_result

  type :: spin_density_cache_entry
    integer(kind=8) :: momentum_revision = 0_8
    integer :: open_size = 0
    logical :: valid = .false.
    complex(kind=8), allocatable :: value(:, :, :)
  end type spin_density_cache_entry

  type :: spin_density_insertion_cache_entry
    integer :: kind = spin_density_no_insertion
    integer :: order = 0
    integer :: identifier = 0
    integer :: correlation_leg = 0
    integer :: open_size = 0
    double precision :: precision_asked = 0d0
    double precision :: precision_found = 0d0
    integer :: return_code = 0
    complex(kind=8), allocatable :: value(:, :, :)
  end type spin_density_insertion_cache_entry

  type :: spin_density_insertion_cache
    integer(kind=8) :: momentum_revision = 0_8
    type(spin_density_insertion_cache_entry), allocatable :: entries(:)
  end type spin_density_insertion_cache

  type(spin_density_cache_entry), allocatable, save :: lo_cache(:, :)
  type(spin_density_insertion_cache), allocatable, save :: &
       insertion_cache(:, :)

  public :: initialize_spin_density_block
  public :: reset_spin_density_caches
  public :: load_cached_lo_density, record_lo_density
  public :: load_cached_spin_density_insertion
  public :: record_spin_density_insertion
  public :: set_spin_density_insertion
  public :: strict_spin_density_product
  public :: multiplicative_spin_density_product
  public :: spin_density_product_order

contains

  subroutine reset_spin_density_caches()
    ! A diagnostic or channel-partition contraction can evaluate several
    ! blocks while one global model coupling is active.  Those matrices are
    ! valid for that contraction, but they must not leak into the subsequent
    ! block-factorized evaluation, where every block is evaluated at its own
    ! immutable reference coupling.  Momentum revisions alone cannot
    ! distinguish the two coupling contexts.
    if (allocated(lo_cache)) deallocate(lo_cache)
    if (allocated(insertion_cache)) deallocate(insertion_cache)
  end subroutine reset_spin_density_caches

  subroutine initialize_spin_density_block(result, event_slot, block, &
                                             open_size)
    type(spin_density_block_result), intent(out) :: result
    integer, intent(in) :: event_slot, block, open_size

    call validate_identity(event_slot, block, open_size)
    result%block = block
    result%event_slot = event_slot
    result%open_size = open_size
  end subroutine initialize_spin_density_block


  subroutine load_cached_lo_density(result, available)
    type(spin_density_block_result), intent(inout) :: result
    logical, intent(out) :: available
    integer(kind=8) :: revision

    call ensure_cache()
    call validate_result_identity(result)
    revision = factorized_block_momentum_revision( &
         result%event_slot, result%block)
    available = revision > 0_8 .and. &
         lo_cache(result%block, result%event_slot)%valid .and. &
         lo_cache(result%block, result%event_slot)%momentum_revision == &
         revision .and. &
         lo_cache(result%block, result%event_slot)%open_size == &
         result%open_size
    if (.not. available) return
    allocate(result%lo(1, result%open_size, result%open_size))
    result%lo = lo_cache(result%block, result%event_slot)%value
    result%has_lo = .true.
  end subroutine load_cached_lo_density


  subroutine record_lo_density(result, density)
    type(spin_density_block_result), intent(inout) :: result
    complex(kind=8), intent(in) :: density(:, :, :)
    integer(kind=8) :: revision

    call ensure_cache()
    call validate_result_identity(result)
    if (size(density, 1) /= 1 .or. &
        size(density, 2) /= result%open_size .or. &
        size(density, 3) /= result%open_size) then
      call fail_spin_density_results('an LO density has the wrong shape')
    end if
    revision = factorized_block_momentum_revision( &
         result%event_slot, result%block)
    if (revision <= 0_8) then
      call fail_spin_density_results( &
           'cannot cache an LO density without boosted block momenta')
    end if
    if (allocated(result%lo)) deallocate(result%lo)
    allocate(result%lo(1, result%open_size, result%open_size))
    result%lo = density
    result%has_lo = .true.

    associate(entry => lo_cache(result%block, result%event_slot))
      if (allocated(entry%value)) deallocate(entry%value)
      allocate(entry%value(1, result%open_size, result%open_size))
      entry%value = density
      entry%momentum_revision = revision
      entry%open_size = result%open_size
      entry%valid = .true.
    end associate
  end subroutine record_lo_density


  subroutine load_cached_spin_density_insertion( &
       result, kind, identifier, correlation_leg, precision_asked, &
       available, precision_found, return_code)
    type(spin_density_block_result), intent(inout) :: result
    integer, intent(in) :: kind, identifier, correlation_leg
    double precision, intent(in) :: precision_asked
    logical, intent(out) :: available
    double precision, intent(out) :: precision_found
    integer, intent(out) :: return_code
    integer(kind=8) :: revision
    integer :: entry_index

    call validate_result_identity(result)
    call validate_insertion_key(kind, identifier, correlation_leg)
    call prepare_insertion_cache(result, revision)
    available = .false.
    precision_found = 0d0
    return_code = 0
    associate(cache => insertion_cache(result%block, result%event_slot))
      if (.not. allocated(cache%entries)) return
      do entry_index = 1, size(cache%entries)
        associate(entry => cache%entries(entry_index))
          if (entry%kind /= kind .or. &
              entry%identifier /= identifier .or. &
              entry%correlation_leg /= correlation_leg .or. &
              entry%open_size /= result%open_size .or. &
              entry%precision_asked /= precision_asked) cycle
          call set_spin_density_insertion( &
               result, entry%kind, entry%order, entry%value)
          precision_found = entry%precision_found
          return_code = entry%return_code
          available = .true.
          return
        end associate
      end do
    end associate
  end subroutine load_cached_spin_density_insertion


  subroutine record_spin_density_insertion( &
       result, kind, order, identifier, correlation_leg, &
       precision_asked, precision_found, return_code, density)
    type(spin_density_block_result), intent(inout) :: result
    integer, intent(in) :: kind, order, identifier, correlation_leg
    double precision, intent(in) :: precision_asked, precision_found
    integer, intent(in) :: return_code
    complex(kind=8), intent(in) :: density(:, :, :)
    type(spin_density_insertion_cache_entry), allocatable :: grown(:)
    integer(kind=8) :: revision
    integer :: entry_count

    call validate_insertion_key(kind, identifier, correlation_leg)
    call prepare_insertion_cache(result, revision)
    call set_spin_density_insertion(result, kind, order, density)
    associate(cache => insertion_cache(result%block, result%event_slot))
      entry_count = 0
      if (allocated(cache%entries)) entry_count = size(cache%entries)
      allocate(grown(entry_count + 1))
      if (entry_count > 0) grown(1:entry_count) = cache%entries
      associate(entry => grown(entry_count + 1))
        entry%kind = kind
        entry%order = order
        entry%identifier = identifier
        entry%correlation_leg = correlation_leg
        entry%open_size = result%open_size
        entry%precision_asked = precision_asked
        entry%precision_found = precision_found
        entry%return_code = return_code
        allocate(entry%value(size(density, 1), size(density, 2), &
                             size(density, 3)))
        entry%value = density
      end associate
      call move_alloc(grown, cache%entries)
    end associate
  end subroutine record_spin_density_insertion


  subroutine set_spin_density_insertion(result, kind, order, density)
    type(spin_density_block_result), intent(inout) :: result
    integer, intent(in) :: kind, order
    complex(kind=8), intent(in) :: density(:, :, :)

    call validate_result_identity(result)
    if (result%has_insertion) then
      call fail_spin_density_results( &
           'a block result received more than one insertion')
    end if
    if (kind <= spin_density_no_insertion .or. &
        kind > spin_density_color_insertion) then
      call fail_spin_density_results('an insertion kind is invalid')
    end if
    if (order < 0 .or. order > 1) then
      call fail_spin_density_results( &
           'only LO or one-NLO-order insertions are supported')
    end if
    if (size(density, 1) < 1 .or. &
        size(density, 2) /= result%open_size .or. &
        size(density, 3) /= result%open_size) then
      call fail_spin_density_results( &
           'an insertion density has the wrong shape')
    end if
    allocate(result%insertion(size(density, 1), result%open_size, &
                              result%open_size))
    result%insertion = density
    result%insertion_kind = kind
    result%insertion_order = order
    result%has_insertion = .true.
  end subroutine set_spin_density_insertion


  complex(kind=8) function strict_spin_density_product( &
       results, active_position, insertion_rank, left_indices, &
       right_indices)
    type(spin_density_block_result), intent(in) :: results(:)
    integer, intent(in) :: active_position, insertion_rank
    integer, intent(in) :: left_indices(:), right_indices(:)
    integer :: insertion_ranks(size(results))

    if (size(left_indices) /= size(results) .or. &
        size(right_indices) /= size(results)) then
      call fail_spin_density_results( &
           'a contraction index vector has the wrong size')
    end if
    if (active_position < 0 .or. active_position > size(results)) then
      call fail_spin_density_results( &
           'a contraction active position is out of range')
    end if
    if (active_position == 0 .and. insertion_rank /= 0) then
      call fail_spin_density_results( &
           'an all-LO contraction requested an insertion rank')
    end if

    insertion_ranks = 0
    if (active_position > 0) then
      if (.not. results(active_position)%has_insertion) then
        call fail_spin_density_results( &
             'the active block has no density insertion')
      end if
      if (results(active_position)%insertion_order > 1) then
        call fail_spin_density_results( &
             'a higher-order density insertion was requested')
      end if
      insertion_ranks(active_position) = insertion_rank
    end if
    strict_spin_density_product = multiplicative_spin_density_product( &
         results, insertion_ranks, left_indices, right_indices)
  end function strict_spin_density_product


  complex(kind=8) function multiplicative_spin_density_product( &
       results, insertion_ranks, left_indices, right_indices)
    type(spin_density_block_result), intent(in) :: results(:)
    integer, intent(in) :: insertion_ranks(:)
    integer, intent(in) :: left_indices(:), right_indices(:)
    integer :: position, rank

    call validate_contraction_vectors(results, insertion_ranks, &
                                      left_indices, right_indices)
    multiplicative_spin_density_product = (1d0, 0d0)
    do position = 1, size(results)
      rank = insertion_ranks(position)
      if (rank == 0) then
        if (.not. results(position)%has_lo) then
          call fail_spin_density_results( &
               'a selected LO block has no cached LO density')
        end if
        multiplicative_spin_density_product = &
             multiplicative_spin_density_product* &
             results(position)%lo(1, left_indices(position), &
                                  right_indices(position))
      else
        if (.not. results(position)%has_insertion) then
          call fail_spin_density_results( &
               'a selected block has no density insertion')
        end if
        if (rank < 1 .or. rank > size(results(position)%insertion, 1)) then
          call fail_spin_density_results( &
               'an insertion rank is out of range')
        end if
        multiplicative_spin_density_product = &
             multiplicative_spin_density_product* &
             results(position)%insertion(rank, left_indices(position), &
                                         right_indices(position))
      end if
    end do
  end function multiplicative_spin_density_product


  integer function spin_density_product_order(results, insertion_ranks)
    type(spin_density_block_result), intent(in) :: results(:)
    integer, intent(in) :: insertion_ranks(:)
    integer :: position

    if (size(insertion_ranks) /= size(results)) then
      call fail_spin_density_results( &
           'a contraction insertion-rank vector has the wrong size')
    end if
    spin_density_product_order = 0
    do position = 1, size(results)
      if (insertion_ranks(position) < 0) then
        call fail_spin_density_results('an insertion rank is negative')
      end if
      if (insertion_ranks(position) == 0) cycle
      if (.not. results(position)%has_insertion) then
        call fail_spin_density_results( &
             'a selected block has no density insertion')
      end if
      if (insertion_ranks(position) > &
          size(results(position)%insertion, 1)) then
        call fail_spin_density_results( &
             'an insertion rank is out of range')
      end if
      spin_density_product_order = spin_density_product_order + &
           results(position)%insertion_order
    end do
  end function spin_density_product_order


  subroutine validate_contraction_vectors(results, insertion_ranks, &
                                           left_indices, right_indices)
    type(spin_density_block_result), intent(in) :: results(:)
    integer, intent(in) :: insertion_ranks(:)
    integer, intent(in) :: left_indices(:), right_indices(:)
    integer :: position

    if (size(insertion_ranks) /= size(results) .or. &
        size(left_indices) /= size(results) .or. &
        size(right_indices) /= size(results)) then
      call fail_spin_density_results( &
           'a contraction vector has the wrong size')
    end if
    do position = 1, size(results)
      if (insertion_ranks(position) < 0) then
        call fail_spin_density_results('an insertion rank is negative')
      end if
      if (left_indices(position) < 1 .or. &
          left_indices(position) > results(position)%open_size .or. &
          right_indices(position) < 1 .or. &
          right_indices(position) > results(position)%open_size) then
        call fail_spin_density_results( &
             'a contraction spin index is out of range')
      end if
    end do
  end subroutine validate_contraction_vectors


  subroutine ensure_cache()
    if (.not. allocated(lo_cache)) then
      allocate(lo_cache(0:nexternal, soft_counterevent:real_event))
    end if
    if (.not. allocated(insertion_cache)) then
      allocate(insertion_cache(0:nexternal, &
                               soft_counterevent:real_event))
    end if
  end subroutine ensure_cache


  subroutine prepare_insertion_cache(result, revision)
    type(spin_density_block_result), intent(in) :: result
    integer(kind=8), intent(out) :: revision

    call ensure_cache()
    call validate_result_identity(result)
    revision = factorized_block_momentum_revision( &
         result%event_slot, result%block)
    if (revision <= 0_8) then
      call fail_spin_density_results( &
           'cannot cache an insertion without boosted block momenta')
    end if
    associate(cache => insertion_cache(result%block, result%event_slot))
      if (cache%momentum_revision /= revision) then
        if (allocated(cache%entries)) deallocate(cache%entries)
        cache%momentum_revision = revision
      end if
    end associate
  end subroutine prepare_insertion_cache


  subroutine validate_insertion_key(kind, identifier, correlation_leg)
    integer, intent(in) :: kind, identifier, correlation_leg

    if (kind <= spin_density_no_insertion .or. &
        kind > spin_density_color_insertion) then
      call fail_spin_density_results('an insertion kind is invalid')
    end if
    if (identifier < 0) then
      call fail_spin_density_results('an insertion identifier is negative')
    end if
    if (correlation_leg < 0) then
      call fail_spin_density_results('a correlation leg is negative')
    end if
  end subroutine validate_insertion_key


  subroutine validate_result_identity(result)
    type(spin_density_block_result), intent(in) :: result

    call validate_identity(result%event_slot, result%block, &
                           result%open_size)
  end subroutine validate_result_identity


  subroutine validate_identity(event_slot, block, open_size)
    integer, intent(in) :: event_slot, block, open_size

    if (event_slot < soft_counterevent .or. event_slot > real_event) then
      call fail_spin_density_results('an event slot is out of range')
    end if
    if (block < 0 .or. block > nexternal) then
      call fail_spin_density_results('a block index is out of range')
    end if
    if (open_size < 1) then
      call fail_spin_density_results('an open-spin size is invalid')
    end if
  end subroutine validate_identity


  subroutine fail_spin_density_results(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') 'ERROR in spin_density_matrix_results: '// &
         trim(message)
    stop 1
  end subroutine fail_spin_density_results

end module spin_density_matrix_results
