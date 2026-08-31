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

  ! The subtraction slots remain an implementation detail of one FKS block.
  ! Once their own matrix-element and phase-space factors have been applied,
  ! a block exposes only the canonical n-body (B) and resolved-real (R)
  ! branches to the decay-tree contraction.
  integer, parameter, public :: spin_density_bornlike_branch = 0
  integer, parameter, public :: spin_density_real_branch = 1

  ! One physical production or decay block.  LO is the reusable spectator
  ! density.  INSERTION is the sole block-local object that may replace it in
  ! a fixed-order contraction (a local underlying Born, real, virtual,
  ! spin-correlated, or colour-correlated density).  The contraction routine
  ! below never reads INSERTION from more than one block.
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

  ! A fully weighted result for one production or decay block.  The first
  ! index is deliberately generic: it may label the central weight and its
  ! reweighting variations, or any other set of weights which must remain
  ! separate until after the global contraction.  BORNLIKE contains the sum
  ! LO+SV+I+S+C+SC on one canonical n-body event; REAL contains the resolved
  ! real contribution on the local (n+1)-body event.
  type, public :: spin_density_branch_result
    integer :: block = -1
    integer :: open_size = 0
    integer :: weight_count = 0
    integer :: bornlike_event_slot = soft_counterevent
    integer :: real_event_slot = real_event
    logical :: has_bornlike = .false.
    logical :: has_real = .false.
    complex(kind=8), allocatable :: bornlike(:, :, :)
    complex(kind=8), allocatable :: real(:, :, :)
  end type spin_density_branch_result

  type :: spin_density_cache_entry
    integer(kind=8) :: momentum_revision = 0_8
    integer :: open_size = 0
    integer :: qcd_power = 0
    double precision :: strong_coupling = 0d0
    logical :: valid = .false.
    complex(kind=8), allocatable :: value(:, :, :)
  end type spin_density_cache_entry

  type(spin_density_cache_entry), allocatable, save :: lo_cache(:, :)

  public :: initialize_spin_density_block
  public :: load_cached_lo_density, record_lo_density
  public :: fetch_cached_lo_density
  public :: set_spin_density_insertion
  public :: selected_spin_density_product
  public :: strict_spin_density_product
  public :: initialize_spin_density_branches
  public :: add_spin_density_branch
  public :: spin_density_branch_product
  public :: spin_density_branch_leaf_count
  public :: decode_spin_density_branch_mask

contains

  subroutine initialize_spin_density_block(result, event_slot, block, &
                                             open_size)
    type(spin_density_block_result), intent(out) :: result
    integer, intent(in) :: event_slot, block, open_size

    call validate_identity(event_slot, block, open_size)
    result%block = block
    result%event_slot = event_slot
    result%open_size = open_size
  end subroutine initialize_spin_density_block


  subroutine load_cached_lo_density( &
       result, qcd_power, strong_coupling, available)
    type(spin_density_block_result), intent(inout) :: result
    integer, intent(in) :: qcd_power
    double precision, intent(in) :: strong_coupling
    logical, intent(out) :: available
    double precision :: coupling_rescaling
    integer(kind=8) :: revision
    integer :: cached_slot

    call ensure_cache()
    call validate_result_identity(result)
    call validate_coupling_key(qcd_power, strong_coupling)
    revision = factorized_block_momentum_revision( &
         result%event_slot, result%block)
    call find_cached_lo_slot( &
         result%block, result%open_size, qcd_power, revision, &
         cached_slot, available)
    if (.not. available) return
    coupling_rescaling = coupling_cache_rescaling( &
         qcd_power, lo_cache(result%block, cached_slot)%strong_coupling, &
         strong_coupling)
    allocate(result%lo(1, result%open_size, result%open_size))
    result%lo = coupling_rescaling* &
         lo_cache(result%block, cached_slot)%value
    result%has_lo = .true.
  end subroutine load_cached_lo_density


  subroutine record_lo_density( &
       result, density, qcd_power, strong_coupling)
    type(spin_density_block_result), intent(inout) :: result
    complex(kind=8), intent(in) :: density(:, :, :)
    integer, intent(in) :: qcd_power
    double precision, intent(in) :: strong_coupling
    integer(kind=8) :: revision

    call ensure_cache()
    call validate_result_identity(result)
    call validate_coupling_key(qcd_power, strong_coupling)
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
      entry%qcd_power = qcd_power
      entry%strong_coupling = strong_coupling
      entry%valid = .true.
    end associate
  end subroutine record_lo_density


  subroutine fetch_cached_lo_density( &
       event_slot, block, open_size, qcd_power, strong_coupling, &
       density, available)
    integer, intent(in) :: event_slot, block, open_size, qcd_power
    double precision, intent(in) :: strong_coupling
    complex(kind=8), intent(out) :: density(1, open_size, open_size)
    logical, intent(out) :: available
    double precision :: coupling_rescaling
    integer(kind=8) :: revision
    integer :: cached_slot

    call ensure_cache()
    call validate_identity(event_slot, block, open_size)
    call validate_coupling_key(qcd_power, strong_coupling)
    revision = factorized_block_momentum_revision(event_slot, block)
    call find_cached_lo_slot( &
         block, open_size, qcd_power, revision, cached_slot, available)
    density = (0d0, 0d0)
    if (available) then
      coupling_rescaling = coupling_cache_rescaling( &
           qcd_power, lo_cache(block, cached_slot)%strong_coupling, &
           strong_coupling)
      density = coupling_rescaling*lo_cache(block, cached_slot)%value
    end if
  end subroutine fetch_cached_lo_density


  subroutine find_cached_lo_slot( &
       block, open_size, qcd_power, revision, cached_slot, available)
    integer, intent(in) :: block, open_size, qcd_power
    integer(kind=8), intent(in) :: revision
    integer, intent(out) :: cached_slot
    logical, intent(out) :: available
    integer :: event_slot

    cached_slot = soft_counterevent
    available = .false.
    if (revision <= 0_8) return
    do event_slot = soft_counterevent, real_event
      if (.not. lo_cache(block, event_slot)%valid) cycle
      if (lo_cache(block, event_slot)%momentum_revision /= revision) cycle
      if (lo_cache(block, event_slot)%open_size /= open_size) cycle
      if (lo_cache(block, event_slot)%qcd_power /= qcd_power) cycle
      cached_slot = event_slot
      available = .true.
      return
    end do
  end subroutine find_cached_lo_slot


  pure double precision function coupling_cache_rescaling( &
       qcd_power, cached_coupling, requested_coupling)
    integer, intent(in) :: qcd_power
    double precision, intent(in) :: cached_coupling, requested_coupling

    if (qcd_power == 0) then
      coupling_cache_rescaling = 1d0
    else
      coupling_cache_rescaling = &
           (requested_coupling/cached_coupling)**qcd_power
    end if
  end function coupling_cache_rescaling


  subroutine validate_coupling_key(qcd_power, strong_coupling)
    integer, intent(in) :: qcd_power
    double precision, intent(in) :: strong_coupling

    if (qcd_power < 0) then
      call fail_spin_density_results( &
           'an LO density has a negative QCD coupling power')
    end if
    if (qcd_power > 0 .and. strong_coupling <= 0d0) then
      call fail_spin_density_results( &
           'an LO density has an invalid strong coupling')
    end if
  end subroutine validate_coupling_key


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


  complex(kind=8) function selected_spin_density_product( &
       results, insertion_ranks, left_indices, right_indices)
    type(spin_density_block_result), intent(in) :: results(:)
    integer, intent(in) :: insertion_ranks(:)
    integer, intent(in) :: left_indices(:), right_indices(:)
    integer :: position, rank

    if (size(insertion_ranks) /= size(results) .or. &
        size(left_indices) /= size(results) .or. &
        size(right_indices) /= size(results)) then
      call fail_spin_density_results( &
           'a selected contraction vector has the wrong size')
    end if

    selected_spin_density_product = (1d0, 0d0)
    do position = 1, size(results)
      call validate_contraction_indices( &
           results(position)%open_size, left_indices(position), &
           right_indices(position))
      rank = insertion_ranks(position)
      if (rank == 0) then
        if (.not. results(position)%has_lo) then
          call fail_spin_density_results( &
               'a selected LO block has no cached LO density')
        end if
        selected_spin_density_product = selected_spin_density_product* &
             results(position)%lo(1, left_indices(position), &
                                  right_indices(position))
      else
        if (.not. results(position)%has_insertion) then
          call fail_spin_density_results( &
               'a selected block has no density insertion')
        end if
        if (results(position)%insertion_order > 1) then
          call fail_spin_density_results( &
               'a higher-order density insertion was requested')
        end if
        if (rank < 1 .or. rank > &
            size(results(position)%insertion, 1)) then
          call fail_spin_density_results( &
               'a selected insertion rank is out of range')
        end if
        selected_spin_density_product = selected_spin_density_product* &
             results(position)%insertion( &
             rank, left_indices(position), right_indices(position))
      end if
    end do
  end function selected_spin_density_product


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
    if (active_position > 0) then
      if (.not. results(active_position)%has_insertion) then
        call fail_spin_density_results( &
             'the active block has no density insertion')
      end if
      if (insertion_rank < 1 .or. insertion_rank > &
          size(results(active_position)%insertion, 1)) then
        call fail_spin_density_results('an insertion rank is out of range')
      end if
    end if

    insertion_ranks = 0
    if (active_position > 0) insertion_ranks(active_position) = &
         insertion_rank
    strict_spin_density_product = selected_spin_density_product( &
         results, insertion_ranks, left_indices, right_indices)
  end function strict_spin_density_product


  subroutine initialize_spin_density_branches(result, block, open_size, &
                                               weight_count)
    type(spin_density_branch_result), intent(inout) :: result
    integer, intent(in) :: block, open_size, weight_count

    call validate_identity(soft_counterevent, block, open_size)
    if (weight_count < 1) then
      call fail_spin_density_results('a branch weight count is invalid')
    end if
    result%block = block
    result%open_size = open_size
    result%weight_count = weight_count
    result%bornlike_event_slot = soft_counterevent
    result%real_event_slot = real_event
    result%has_bornlike = .false.
    result%has_real = .false.
    if (allocated(result%bornlike)) then
      if (any(shape(result%bornlike) /= &
              [weight_count, open_size, open_size])) then
        deallocate(result%bornlike)
      end if
    end if
    if (allocated(result%real)) then
      if (any(shape(result%real) /= &
              [weight_count, open_size, open_size])) then
        deallocate(result%real)
      end if
    end if
    if (.not. allocated(result%bornlike)) &
         allocate(result%bornlike(weight_count, open_size, open_size))
    if (.not. allocated(result%real)) &
         allocate(result%real(weight_count, open_size, open_size))
    result%bornlike = (0d0, 0d0)
    result%real = (0d0, 0d0)
  end subroutine initialize_spin_density_branches


  subroutine add_spin_density_branch(result, branch, density)
    type(spin_density_branch_result), intent(inout) :: result
    integer, intent(in) :: branch
    complex(kind=8), intent(in) :: density(:, :, :)

    call validate_branch_result(result)
    if (size(density, 1) /= result%weight_count .or. &
        size(density, 2) /= result%open_size .or. &
        size(density, 3) /= result%open_size) then
      call fail_spin_density_results('a branch density has the wrong shape')
    end if
    select case (branch)
    case (spin_density_bornlike_branch)
      result%bornlike = result%bornlike + density
      result%has_bornlike = .true.
    case (spin_density_real_branch)
      result%real = result%real + density
      result%has_real = .true.
    case default
      call fail_spin_density_results('a spin-density branch is invalid')
    end select
  end subroutine add_spin_density_branch


  complex(kind=8) function spin_density_branch_product( &
       results, branch_by_position, weight_rank, left_indices, &
       right_indices)
    type(spin_density_branch_result), intent(in) :: results(:)
    integer, intent(in) :: branch_by_position(:), weight_rank
    integer, intent(in) :: left_indices(:), right_indices(:)
    integer :: position

    if (size(branch_by_position) /= size(results) .or. &
        size(left_indices) /= size(results) .or. &
        size(right_indices) /= size(results)) then
      call fail_spin_density_results( &
           'a branch contraction vector has the wrong size')
    end if
    spin_density_branch_product = (1d0, 0d0)
    do position = 1, size(results)
      call validate_branch_result(results(position))
      call validate_contraction_indices( &
           results(position)%open_size, left_indices(position), &
           right_indices(position))
      if (weight_rank < 1 .or. &
          weight_rank > results(position)%weight_count) then
        call fail_spin_density_results( &
             'a branch contraction weight rank is out of range')
      end if
      select case (branch_by_position(position))
      case (spin_density_bornlike_branch)
        if (.not. results(position)%has_bornlike) then
          call fail_spin_density_results( &
               'a selected block has no Born-like branch')
        end if
        spin_density_branch_product = spin_density_branch_product* &
             results(position)%bornlike( &
             weight_rank, left_indices(position), right_indices(position))
      case (spin_density_real_branch)
        if (.not. results(position)%has_real) then
          call fail_spin_density_results( &
               'a selected block has no real branch')
        end if
        spin_density_branch_product = spin_density_branch_product* &
             results(position)%real( &
             weight_rank, left_indices(position), right_indices(position))
      case default
        call fail_spin_density_results( &
             'a branch contraction selected an invalid branch')
      end select
    end do
  end function spin_density_branch_product


  integer(kind=8) function spin_density_branch_leaf_count(block_count)
    integer, intent(in) :: block_count

    if (block_count < 0 .or. block_count >= bit_size(0_8) - 1) then
      call fail_spin_density_results( &
           'a branch leaf count would overflow its bit mask')
    end if
    spin_density_branch_leaf_count = shiftl(1_8, block_count)
  end function spin_density_branch_leaf_count


  subroutine decode_spin_density_branch_mask(mask, branch_by_position)
    integer(kind=8), intent(in) :: mask
    integer, intent(out) :: branch_by_position(:)
    integer(kind=8) :: leaf_count
    integer :: position

    leaf_count = spin_density_branch_leaf_count(size(branch_by_position))
    if (mask < 0_8 .or. mask >= leaf_count) then
      call fail_spin_density_results('a branch mask is out of range')
    end if
    do position = 1, size(branch_by_position)
      branch_by_position(position) = merge( &
           spin_density_real_branch, spin_density_bornlike_branch, &
           btest(mask, position - 1))
    end do
  end subroutine decode_spin_density_branch_mask


  subroutine ensure_cache()
    if (allocated(lo_cache)) return
    allocate(lo_cache(0:nexternal, soft_counterevent:real_event))
  end subroutine ensure_cache


  subroutine validate_result_identity(result)
    type(spin_density_block_result), intent(in) :: result

    call validate_identity(result%event_slot, result%block, &
                           result%open_size)
  end subroutine validate_result_identity


  subroutine validate_branch_result(result)
    type(spin_density_branch_result), intent(in) :: result

    call validate_identity(soft_counterevent, result%block, &
                           result%open_size)
    if (result%weight_count < 1 .or. &
        .not. allocated(result%bornlike) .or. &
        .not. allocated(result%real)) then
      call fail_spin_density_results( &
           'a spin-density branch result is not initialized')
    end if
    if (any(shape(result%bornlike) /= &
            [result%weight_count, result%open_size, result%open_size]) .or. &
        any(shape(result%real) /= &
            [result%weight_count, result%open_size, result%open_size])) then
      call fail_spin_density_results( &
           'a spin-density branch result has inconsistent storage')
    end if
  end subroutine validate_branch_result


  subroutine validate_contraction_indices(open_size, left_index, &
                                          right_index)
    integer, intent(in) :: open_size, left_index, right_index

    if (left_index < 1 .or. left_index > open_size .or. &
        right_index < 1 .or. right_index > open_size) then
      call fail_spin_density_results( &
           'a spin-density contraction index is out of range')
    end if
  end subroutine validate_contraction_indices


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
