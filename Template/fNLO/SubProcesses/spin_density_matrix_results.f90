module spin_density_matrix_results
  use process_dimensions, only: nexternal
  use fnlo_process_common, only: soft_counterevent, real_event
  use factorized_phase_space, only: &
       factorized_block_momentum_revision
  use multiplicative_scale_state, only: &
       multiplicative_active_coupling_context
  implicit none
  private

  integer, parameter, public :: spin_density_no_insertion = 0
  integer, parameter, public :: spin_density_born_insertion = 1
  integer, parameter, public :: spin_density_real_insertion = 2
  integer, parameter, public :: spin_density_virtual_insertion = 3
  integer, parameter, public :: spin_density_color_insertion = 4
  integer, parameter, public :: spin_density_fast_virtual_insertion = 5

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
    integer(kind=8) :: coupling_context = 0_8
    logical :: has_lo = .false.
    logical :: has_insertion = .false.
    complex(kind=8), allocatable :: lo(:, :, :)
    complex(kind=8), allocatable :: insertion(:, :, :)
  end type spin_density_block_result

  type :: spin_density_cache_entry
    integer(kind=8) :: momentum_revision = 0_8
    integer(kind=8) :: coupling_context = 0_8
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
    integer(kind=8) :: coupling_context = 0_8
    type(spin_density_insertion_cache_entry), allocatable :: entries(:)
  end type spin_density_insertion_cache

  type, public :: spin_density_cache_statistics
    integer(kind=8) :: raw_amplitude_hits = 0_8
    integer(kind=8) :: raw_amplitude_misses = 0_8
    integer(kind=8) :: raw_amplitude_provider_evaluations = 0_8
    integer(kind=8) :: lo_hits = 0_8
    integer(kind=8) :: lo_misses = 0_8
    integer(kind=8) :: insertion_hits = 0_8
    integer(kind=8) :: insertion_misses = 0_8
    integer(kind=8) :: lo_provider_evaluations = 0_8
    integer(kind=8) :: insertion_provider_evaluations = 0_8
    integer(kind=8) :: virtual_madloop_evaluations = 0_8
    integer(kind=8) :: virtual_direct_reconstructions = 0_8
    integer(kind=8) :: virtual_direct_fallbacks = 0_8
    integer(kind=8) :: virtual_tomography_reconstructions = 0_8
    integer(kind=8) :: exact_family_candidates = 0_8
    integer(kind=8) :: exact_family_cut_rejections = 0_8
    integer(kind=8) :: exact_family_acceptances = 0_8
    integer(kind=8) :: density_contractions = 0_8
    integer(kind=8) :: scale_reweight_evaluations = 0_8
    integer(kind=8) :: scale_luminosity_evaluations = 0_8
    integer(kind=8) :: scale_luminosity_cache_hits = 0_8
    integer(kind=8) :: pdf_member_initializations = 0_8
    integer(kind=8) :: pdf_luminosity_evaluations = 0_8
    integer(kind=8) :: pdf_luminosity_cache_hits = 0_8
    integer(kind=8) :: histogram_family_fills = 0_8
  end type spin_density_cache_statistics

  type(spin_density_cache_entry), allocatable, save :: lo_cache(:, :)
  type(spin_density_insertion_cache), allocatable, save :: &
       insertion_cache(:, :)
  type(spin_density_cache_statistics), save :: cache_statistics

  public :: initialize_spin_density_block
  public :: reset_spin_density_caches
  public :: reset_spin_density_cache_statistics
  public :: fetch_spin_density_cache_statistics
  public :: record_raw_amplitude_cache_hit
  public :: record_raw_amplitude_cache_miss
  public :: record_direct_virtual_reconstruction
  public :: record_virtual_tomography_reconstruction
  public :: record_exact_family_candidate
  public :: record_exact_family_cut_rejection
  public :: record_exact_family_acceptance
  public :: record_density_contraction
  public :: record_scale_reweight_evaluation
  public :: record_scale_luminosity_evaluation
  public :: record_scale_luminosity_cache_hit
  public :: record_pdf_member_initialization
  public :: record_pdf_luminosity_evaluation
  public :: record_pdf_luminosity_cache_hit
  public :: record_histogram_family_fill
  public :: load_cached_lo_density, record_lo_density
  public :: load_cached_spin_density_insertion
  public :: record_spin_density_insertion
  public :: set_spin_density_insertion
  public :: strict_spin_density_product
  public :: multiplicative_spin_density_product
  public :: spin_density_product_order

contains

  subroutine record_raw_amplitude_cache_hit()
    cache_statistics%raw_amplitude_hits = &
         cache_statistics%raw_amplitude_hits + 1_8
  end subroutine record_raw_amplitude_cache_hit


  subroutine record_raw_amplitude_cache_miss()
    cache_statistics%raw_amplitude_misses = &
         cache_statistics%raw_amplitude_misses + 1_8
    cache_statistics%raw_amplitude_provider_evaluations = &
         cache_statistics%raw_amplitude_provider_evaluations + 1_8
  end subroutine record_raw_amplitude_cache_miss


  subroutine record_direct_virtual_reconstruction(success, evaluations)
    logical, intent(in) :: success
    integer, intent(in) :: evaluations

    if (evaluations < 1) then
      call fail_spin_density_results( &
           'a direct virtual reconstruction recorded no evaluations')
    end if
    cache_statistics%virtual_madloop_evaluations = &
         cache_statistics%virtual_madloop_evaluations + int(evaluations, 8)
    if (success) then
      cache_statistics%virtual_direct_reconstructions = &
           cache_statistics%virtual_direct_reconstructions + 1_8
    else
      cache_statistics%virtual_direct_fallbacks = &
           cache_statistics%virtual_direct_fallbacks + 1_8
    end if
  end subroutine record_direct_virtual_reconstruction


  subroutine record_virtual_tomography_reconstruction(evaluations)
    integer, intent(in) :: evaluations

    if (evaluations < 1) then
      call fail_spin_density_results( &
           'a virtual tomography recorded no evaluations')
    end if
    cache_statistics%virtual_madloop_evaluations = &
         cache_statistics%virtual_madloop_evaluations + int(evaluations, 8)
    cache_statistics%virtual_tomography_reconstructions = &
         cache_statistics%virtual_tomography_reconstructions + 1_8
  end subroutine record_virtual_tomography_reconstruction


  subroutine record_exact_family_candidate()
    cache_statistics%exact_family_candidates = &
         cache_statistics%exact_family_candidates + 1_8
  end subroutine record_exact_family_candidate


  subroutine record_exact_family_cut_rejection()
    cache_statistics%exact_family_cut_rejections = &
         cache_statistics%exact_family_cut_rejections + 1_8
  end subroutine record_exact_family_cut_rejection


  subroutine record_exact_family_acceptance()
    cache_statistics%exact_family_acceptances = &
         cache_statistics%exact_family_acceptances + 1_8
  end subroutine record_exact_family_acceptance


  subroutine record_density_contraction()
    cache_statistics%density_contractions = &
         cache_statistics%density_contractions + 1_8
  end subroutine record_density_contraction


  subroutine record_scale_reweight_evaluation()
    cache_statistics%scale_reweight_evaluations = &
         cache_statistics%scale_reweight_evaluations + 1_8
  end subroutine record_scale_reweight_evaluation


  subroutine record_scale_luminosity_evaluation()
    cache_statistics%scale_luminosity_evaluations = &
         cache_statistics%scale_luminosity_evaluations + 1_8
  end subroutine record_scale_luminosity_evaluation


  subroutine record_scale_luminosity_cache_hit()
    cache_statistics%scale_luminosity_cache_hits = &
         cache_statistics%scale_luminosity_cache_hits + 1_8
  end subroutine record_scale_luminosity_cache_hit


  subroutine record_pdf_member_initialization()
    cache_statistics%pdf_member_initializations = &
         cache_statistics%pdf_member_initializations + 1_8
  end subroutine record_pdf_member_initialization


  subroutine record_pdf_luminosity_evaluation()
    cache_statistics%pdf_luminosity_evaluations = &
         cache_statistics%pdf_luminosity_evaluations + 1_8
  end subroutine record_pdf_luminosity_evaluation


  subroutine record_pdf_luminosity_cache_hit()
    cache_statistics%pdf_luminosity_cache_hits = &
         cache_statistics%pdf_luminosity_cache_hits + 1_8
  end subroutine record_pdf_luminosity_cache_hit


  subroutine record_histogram_family_fill()
    cache_statistics%histogram_family_fills = &
         cache_statistics%histogram_family_fills + 1_8
  end subroutine record_histogram_family_fill

  subroutine reset_spin_density_cache_statistics()
    cache_statistics = spin_density_cache_statistics()
  end subroutine reset_spin_density_cache_statistics


  subroutine fetch_spin_density_cache_statistics(statistics)
    type(spin_density_cache_statistics), intent(out) :: statistics

    statistics = cache_statistics
  end subroutine fetch_spin_density_cache_statistics

  subroutine reset_spin_density_caches()
    ! Explicit reset remains available for process changes and diagnostics.
    ! Ordinary block-reference changes are distinguished by the cache key.
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
    result%coupling_context = multiplicative_active_coupling_context()
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
         lo_cache(result%block, result%event_slot)%coupling_context == &
         result%coupling_context .and. &
         lo_cache(result%block, result%event_slot)%open_size == &
         result%open_size
    if (.not. available) then
      cache_statistics%lo_misses = cache_statistics%lo_misses + 1_8
      return
    end if
    cache_statistics%lo_hits = cache_statistics%lo_hits + 1_8
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
    cache_statistics%lo_provider_evaluations = &
         cache_statistics%lo_provider_evaluations + 1_8

    associate(entry => lo_cache(result%block, result%event_slot))
      if (allocated(entry%value)) deallocate(entry%value)
      allocate(entry%value(1, result%open_size, result%open_size))
      entry%value = density
      entry%momentum_revision = revision
      entry%coupling_context = result%coupling_context
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
      if (allocated(cache%entries)) then
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
            cache_statistics%insertion_hits = &
                 cache_statistics%insertion_hits + 1_8
            return
          end associate
        end do
      end if
    end associate
    cache_statistics%insertion_misses = &
         cache_statistics%insertion_misses + 1_8
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
    cache_statistics%insertion_provider_evaluations = &
         cache_statistics%insertion_provider_evaluations + 1_8
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
        kind > spin_density_fast_virtual_insertion) then
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
      if (cache%momentum_revision /= revision .or. &
          cache%coupling_context /= result%coupling_context) then
        if (allocated(cache%entries)) deallocate(cache%entries)
        cache%momentum_revision = revision
        cache%coupling_context = result%coupling_context
      end if
    end associate
  end subroutine prepare_insertion_cache


  subroutine validate_insertion_key(kind, identifier, correlation_leg)
    integer, intent(in) :: kind, identifier, correlation_leg

    if (kind <= spin_density_no_insertion .or. &
        kind > spin_density_fast_virtual_insertion) then
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
