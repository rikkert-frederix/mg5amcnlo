module factorized_phase_space
  use process_dimensions, only: nexternal, validate_process_dimensions
  use fnlo_process_common, only: soft_counterevent, real_event
  implicit none
  private

  integer(kind=8), save :: next_block_momentum_revision = 0_8

  ! Particle identities are part of a matrix-element block, just as much as
  ! its four-vectors.  Keep them beside the boosted block cache.  The live
  ! local layout may subsequently change to a real-emission context (an ISR
  ! soft projection can even change an incoming flavour), while the reduced
  ! Born block in slot zero must retain its own immutable layout.
  ! Embedding storage is deliberately distinct from the matrix-element
  ! block cache above.  In particular, a soft projection can have a
  ! real-context particle layout while the Born density matrix in slot zero
  ! must retain its Born layout.  Event materialization consumes this cache;
  ! matrix elements never do.
  ! Matrix elements and subtraction kernels deliberately use different
  ! representations of the same block.  The former needs the block boosted
  ! into the event frame, whereas the latter is most naturally evaluated in
  ! the frame in which its three FKS variables were generated.
  ! Tuple realization starts from an immutable block-local representation.
  ! For a decay this is the parent rest frame; for production it is the
  ! generated partonic event frame.  It must remain distinct from both the
  ! boosted matrix-element cache and the subtraction-kernel cache: selecting
  ! a different ancestor term changes the boost of every descendant without
  ! changing the descendant's independently generated local configuration.
  integer, parameter, public :: factorized_no_target = 0
  integer, parameter, public :: factorized_visible_target = 1
  integer, parameter, public :: factorized_block_target = 2

  type, public :: factorized_radiation_state
    double precision :: jacobian = -1d0
    double precision :: radiation_jacobian = -1d0
    double precision :: radiation_weight = -1d0
    double precision :: xi = -1d0
    double precision :: y = -2d0
    double precision :: xi_hat = -1d0
    double precision :: xi_max = -1d0
    double precision :: xi_norm = -1d0
    double precision :: fks_momentum(0:3) = -1d0
    double precision :: shat = -1d0
    double precision :: sqrt_shat = -1d0
    double precision :: y_to_cm = 0d0
    double precision :: bjorken_x(2) = -1d0
    double precision :: y_to_lab = 0d0
  end type factorized_radiation_state

  ! The underlying-Born measure is a product of one independently generated
  ! production block and one factor for every decay node.  Event measures
  ! contain only the event-dependent radiation map, flux and incoming MINT
  ! weight.  Keeping the two layers separate lets a block be regenerated or
  ! replaced without passing a cumulative Jacobian through the decay tree.
  type, public :: factorized_measure_state
    double precision :: jacobian = 1d0
    double precision :: phase_space_weight = 1d0
  end type factorized_measure_state

  ! All data belonging to one (block,event-slot) identity live in one value.
  ! This replaces the former two dozen parallel global arrays and makes a
  ! snapshot an intrinsic assignment.  The process dimensions are initialized
  ! at run time, so the particle buffers are allocatable components which are
  ! allocated once when the persistent block store is created.
  type :: factorized_slot_state
    double precision, allocatable :: momenta(:,:)
    integer :: particle_count = 0
    logical :: is_valid = .false.
    integer(kind=8) :: momentum_revision = 0_8
    integer :: matrix_count = 0
    integer, allocatable :: matrix_pdg(:)
    logical, allocatable :: matrix_final(:)
    logical :: matrix_layout_valid = .false.
    double precision, allocatable :: embedded(:,:)
    integer :: embedded_count = 0
    logical :: embedded_valid = .false.
    double precision, allocatable :: kernel(:,:)
    integer :: kernel_count = 0
    logical :: kernel_valid = .false.
    double precision, allocatable :: local(:,:)
    integer :: local_count = 0
    logical :: local_valid = .false.
    integer, allocatable :: local_pdg(:)
    logical, allocatable :: local_final(:)
    integer, allocatable :: local_target_kind(:)
    integer, allocatable :: local_target_id(:)
    logical :: local_layout_valid = .false.
    type(factorized_radiation_state) :: radiation
    logical :: radiation_valid = .false.
    type(factorized_measure_state) :: measure
    logical :: measure_valid = .false.
  end type factorized_slot_state

  type :: factorized_block_state
    integer :: block = -1
    logical :: initialized = .false.
    type(factorized_slot_state) :: slot( &
         soft_counterevent:real_event)
    type(factorized_measure_state) :: base
    logical :: base_valid = .false.
  end type factorized_block_state

  ! Sequential generation of independent radiative blocks must not let the
  ! later block overwrite an earlier block's real/counterevent family.  A
  ! snapshot owns one physical block in every event slot and can therefore be
  ! overlaid on a common Born baseline without copying or pre-summing events.
  type, public :: factorized_block_snapshot
    logical :: initialized = .false.
    type(factorized_block_state) :: state
  end type factorized_block_snapshot

  type(factorized_block_state), allocatable, save :: block_state(:)

  public :: reset_factorized_phase_space
  public :: store_factorized_block_momenta
  public :: fetch_factorized_block_momenta
  public :: fetch_factorized_matrix_momenta
  public :: fetch_factorized_ordered_matrix_momenta
  public :: factorized_block_momentum_revision
  public :: store_factorized_embedded_momenta
  public :: fetch_factorized_embedded_momenta
  public :: store_factorized_kernel_momenta
  public :: fetch_factorized_kernel_momenta
  public :: store_factorized_local_momenta
  public :: fetch_factorized_local_momenta
  public :: factorized_local_count
  public :: store_factorized_local_layout
  public :: fetch_factorized_local_layout
  public :: store_factorized_radiation_state
  public :: fetch_factorized_radiation_state
  public :: scale_factorized_radiation_jacobians
  public :: store_factorized_base_measure
  public :: multiply_factorized_base_measure
  public :: fetch_factorized_base_measure
  public :: store_factorized_event_measure
  public :: multiply_factorized_event_measure
  public :: fetch_factorized_event_measure
  public :: compose_factorized_base_measure
  public :: compose_factorized_event_measure
  public :: compose_factorized_tuple_measure
  public :: scale_factorized_event_measures
  public :: capture_factorized_block, restore_factorized_block

contains

  subroutine reset_factorized_phase_space()
    integer :: block

    call ensure_storage()
    do block = 0, nexternal
      call initialize_block_state(block_state(block), block)
    end do
  end subroutine reset_factorized_phase_space


  subroutine capture_factorized_block(block, snapshot)
    integer, intent(in) :: block
    type(factorized_block_snapshot), intent(out) :: snapshot

    call ensure_storage()
    call validate_block(block)
    snapshot%state = block_state(block)
    snapshot%initialized = .true.
  end subroutine capture_factorized_block


  subroutine restore_factorized_block(snapshot)
    type(factorized_block_snapshot), intent(in) :: snapshot
    integer :: block

    call ensure_storage()
    if (.not. snapshot%initialized) then
      call fail_factorized_phase_space('a block snapshot is uninitialized')
    end if
    block = snapshot%state%block
    call validate_block(block)
    if (.not. snapshot%state%initialized) then
      call fail_factorized_phase_space('a block snapshot has no state')
    end if
    block_state(block) = snapshot%state
  end subroutine restore_factorized_block


  subroutine store_factorized_block_momenta(event_slot, block, &
                                             particle_count, momenta)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(in) :: momenta(0:, :)
    integer :: cached_block, cached_slot

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    if (size(momenta, 1) < 4 .or. size(momenta, 2) < particle_count) then
      call fail_factorized_phase_space( &
           'a block momentum array has inconsistent bounds')
    end if

    ! Tuple realization is deliberately repeated for every Cartesian NLO
    ! atom, but many atoms differ only in a sibling block or in the density
    ! insertion selected for an otherwise identical block.  Preserve the
    ! momentum identity in that case.  Besides making the revision describe
    ! the matrix-element input rather than the number of realization calls,
    ! this lets all block-local Born/real/virtual densities share one cache.
    ! The matrix layout is immutable once these momenta have been captured;
    ! a subsequently installed real-context layout must not invalidate the
    ! reduced-Born identity of an unchanged soft block.
    if (block_state(block)%slot(event_slot)%is_valid .and. &
        block_state(block)%slot(event_slot)%particle_count == &
        particle_count) then
      if (all(block_state(block)%slot(event_slot)% &
              momenta(:, 1:particle_count) == &
              momenta(0:3, 1:particle_count))) return
    end if
    block_state(block)%slot(event_slot)%momenta = 0d0
    block_state(block)%slot(event_slot)%momenta(:, 1:particle_count) = &
         momenta(0:3, 1:particle_count)
    block_state(block)%slot(event_slot)%particle_count = particle_count
    block_state(block)%slot(event_slot)%is_valid = .true.
    block_state(block)%slot(event_slot)%matrix_count = 0
    block_state(block)%slot(event_slot)%matrix_pdg = 0
    block_state(block)%slot(event_slot)%matrix_final = .false.
    block_state(block)%slot(event_slot)%matrix_layout_valid = .false.
    if (block_state(block)%slot(event_slot)%local_layout_valid .and. &
        block_state(block)%slot(event_slot)%local_count >= &
        particle_count) then
      call capture_matrix_layout(event_slot, block, particle_count)
    end if
    next_block_momentum_revision = next_block_momentum_revision + 1_8
    if (next_block_momentum_revision <= 0_8) then
      ! A wrap is fantastically unlikely, but resetting all identities is
      ! safer than allowing a stale density-matrix cache entry to match.
      next_block_momentum_revision = 1_8
      do cached_block = 0, nexternal
        do cached_slot = soft_counterevent, real_event
          block_state(cached_block)%slot(cached_slot)% &
               momentum_revision = 0_8
        end do
      end do
    end if
    block_state(block)%slot(event_slot)%momentum_revision = &
         next_block_momentum_revision
  end subroutine store_factorized_block_momenta


  subroutine fetch_factorized_block_momenta(event_slot, block, &
                                             particle_count, momenta, &
                                             available)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(out) :: momenta(0:3, particle_count)
    logical, intent(out) :: available

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    available = block_state(block)%slot(event_slot)%is_valid .and. &
         block_state(block)%slot(event_slot)%particle_count == particle_count
    momenta = 0d0
    if (available) then
      momenta = block_state(block)%slot(event_slot)% &
           momenta(:, 1:particle_count)
    end if
  end subroutine fetch_factorized_block_momenta


  subroutine fetch_factorized_matrix_momenta(event_slot, block, &
                                              particle_count, momenta, &
                                              available)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(out) :: momenta(0:3, particle_count)
    logical, intent(out) :: available

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    ! A real or counterevent local layout retains the emitted leg so that
    ! cuts and observables see that atom's own kinematics.  Born, virtual,
    ! soft and collinear density insertions use the reduced Born prefix of
    ! precisely that already-boosted layout.  The exporter guarantees that
    ! the extra FKS leg is appended after all Born legs.
    available = block_state(block)%slot(event_slot)%is_valid .and. &
         block_state(block)%slot(event_slot)%particle_count >= particle_count
    momenta = 0d0
    if (available) then
      momenta = block_state(block)%slot(event_slot)% &
           momenta(:, 1:particle_count)
    end if
  end subroutine fetch_factorized_matrix_momenta


  subroutine fetch_factorized_ordered_matrix_momenta( &
       event_slot, block, particle_count, expected_pdgs, expected_final, &
       momenta, available)
    integer, intent(in) :: event_slot, block, particle_count
    integer, intent(in) :: expected_pdgs(particle_count)
    integer, intent(in) :: expected_final(particle_count)
    double precision, intent(out) :: momenta(0:3, particle_count)
    logical, intent(out) :: available

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    momenta = 0d0
    available = block_state(block)%slot(event_slot)%is_valid
    if (.not. available) return
    if (any(expected_final /= 0 .and. expected_final /= 1)) then
      call fail_factorized_phase_space( &
           'an ordered matrix-element layout has an invalid state')
    end if

    ! First try the event-specific identity.  Real atoms need precisely this
    ! layout.  A reduced collinear/soft-collinear atom, however, can retain a
    ! real-context local layout whose incoming flavour differs from its
    ! underlying Born matrix element.  Spectator blocks likewise have no
    ! event-specific identity at all.  In both cases retry with the immutable
    ! underlying-Born identity while keeping every returned four-vector in
    ! the requested, already-boosted event slot.
    call match_ordered_matrix_layout( &
         event_slot, block, particle_count, expected_pdgs, expected_final, &
         event_slot, momenta, available)
    if (available .or. event_slot == soft_counterevent) return
    call match_ordered_matrix_layout( &
         event_slot, block, particle_count, expected_pdgs, expected_final, &
         soft_counterevent, momenta, available)
  end subroutine fetch_factorized_ordered_matrix_momenta


  subroutine match_ordered_matrix_layout( &
       event_slot, block, particle_count, expected_pdgs, expected_final, &
       layout_slot, momenta, available)
    integer, intent(in) :: event_slot, block, particle_count, layout_slot
    integer, intent(in) :: expected_pdgs(particle_count)
    integer, intent(in) :: expected_final(particle_count)
    double precision, intent(out) :: momenta(0:3, particle_count)
    logical, intent(out) :: available
    logical :: used(nexternal)
    integer :: source_count, expected, source

    momenta = 0d0
    available = block_state(block)%slot(layout_slot)%matrix_layout_valid
    if (.not. available) return
    source_count = min( &
         block_state(block)%slot(event_slot)%particle_count, &
         block_state(block)%slot(layout_slot)%matrix_count)
    if (source_count < particle_count) then
      available = .false.
      return
    end if

    ! Match by signed PDG and incoming/final state.  Stable first-match
    ! ordering preserves repeated-particle identities and also converts
    ! between user-process and FKS leg orderings such as t -> W b / t -> b W.
    used = .false.
    do expected = 1, particle_count
      available = .false.
      do source = 1, source_count
        if (used(source)) cycle
        if (block_state(block)%slot(layout_slot)%matrix_pdg(source) /= &
            expected_pdgs(expected)) cycle
        if (block_state(block)%slot(layout_slot)%matrix_final(source) .neqv. &
            (expected_final(expected) == 1)) cycle
        momenta(:, expected) = &
             block_state(block)%slot(event_slot)%momenta(:, source)
        used(source) = .true.
        available = .true.
        exit
      end do
      if (.not. available) then
        momenta = 0d0
        return
      end if
    end do
  end subroutine match_ordered_matrix_layout


  integer(kind=8) function factorized_block_momentum_revision( &
       event_slot, block)
    integer, intent(in) :: event_slot, block

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    if (block_state(block)%slot(event_slot)%is_valid) then
      factorized_block_momentum_revision = &
           block_state(block)%slot(event_slot)%momentum_revision
    else
      factorized_block_momentum_revision = 0_8
    end if
  end function factorized_block_momentum_revision


  subroutine store_factorized_embedded_momenta(event_slot, block, &
                                                particle_count, momenta)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(in) :: momenta(0:, :)

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    if (size(momenta, 1) < 4 .or. size(momenta, 2) < particle_count) then
      call fail_factorized_phase_space( &
           'an embedded momentum array has inconsistent bounds')
    end if
    block_state(block)%slot(event_slot)%embedded = 0d0
    block_state(block)%slot(event_slot)% &
         embedded(:, 1:particle_count) = &
         momenta(0:3, 1:particle_count)
    block_state(block)%slot(event_slot)%embedded_count = particle_count
    block_state(block)%slot(event_slot)%embedded_valid = .true.
  end subroutine store_factorized_embedded_momenta


  subroutine fetch_factorized_embedded_momenta(event_slot, block, &
                                                particle_count, momenta, &
                                                available)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(out) :: momenta(0:3, particle_count)
    logical, intent(out) :: available

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    available = block_state(block)%slot(event_slot)%embedded_valid .and. &
         block_state(block)%slot(event_slot)%embedded_count == particle_count
    momenta = 0d0
    if (available) then
      momenta = block_state(block)%slot(event_slot)% &
           embedded(:, 1:particle_count)
    end if
  end subroutine fetch_factorized_embedded_momenta


  subroutine store_factorized_kernel_momenta(event_slot, block, &
                                              particle_count, momenta)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(in) :: momenta(0:, :)

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    if (size(momenta, 1) < 4 .or. size(momenta, 2) < particle_count) then
      call fail_factorized_phase_space( &
           'a kernel momentum array has inconsistent bounds')
    end if
    block_state(block)%slot(event_slot)%kernel = 0d0
    block_state(block)%slot(event_slot)%kernel(:, 1:particle_count) = &
         momenta(0:3, 1:particle_count)
    block_state(block)%slot(event_slot)%kernel_count = particle_count
    block_state(block)%slot(event_slot)%kernel_valid = .true.
  end subroutine store_factorized_kernel_momenta


  subroutine fetch_factorized_kernel_momenta(event_slot, block, &
                                              particle_count, momenta, &
                                              available)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(out) :: momenta(0:3, particle_count)
    logical, intent(out) :: available

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    available = block_state(block)%slot(event_slot)%kernel_valid .and. &
         block_state(block)%slot(event_slot)%kernel_count == particle_count
    momenta = 0d0
    if (available) then
      momenta = block_state(block)%slot(event_slot)% &
           kernel(:, 1:particle_count)
    end if
  end subroutine fetch_factorized_kernel_momenta


  subroutine store_factorized_local_momenta(event_slot, block, &
                                             particle_count, momenta)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(in) :: momenta(0:, :)

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    if (size(momenta, 1) < 4 .or. size(momenta, 2) < particle_count) then
      call fail_factorized_phase_space( &
           'a block-local momentum array has inconsistent bounds')
    end if
    block_state(block)%slot(event_slot)%local = 0d0
    block_state(block)%slot(event_slot)%local(:, 1:particle_count) = &
         momenta(0:3, 1:particle_count)
    block_state(block)%slot(event_slot)%local_count = particle_count
    block_state(block)%slot(event_slot)%local_valid = .true.
  end subroutine store_factorized_local_momenta


  subroutine fetch_factorized_local_momenta(event_slot, block, &
                                             particle_count, momenta, &
                                             available)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(out) :: momenta(0:3, particle_count)
    logical, intent(out) :: available

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    available = block_state(block)%slot(event_slot)%local_valid .and. &
         block_state(block)%slot(event_slot)%local_count == particle_count
    momenta = 0d0
    if (available) then
      momenta = block_state(block)%slot(event_slot)% &
           local(:, 1:particle_count)
    end if
  end subroutine fetch_factorized_local_momenta


  integer function factorized_local_count(event_slot, block)
    integer, intent(in) :: event_slot, block

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    if (block_state(block)%slot(event_slot)%local_valid) then
      factorized_local_count = &
           block_state(block)%slot(event_slot)%local_count
    else
      factorized_local_count = 0
    end if
  end function factorized_local_count


  subroutine store_factorized_local_layout( &
       event_slot, block, particle_count, pdgs, particle_is_final, &
       target_kinds, target_ids)
    integer, intent(in) :: event_slot, block, particle_count
    integer, intent(in) :: pdgs(:), target_kinds(:), target_ids(:)
    logical, intent(in) :: particle_is_final(:)
    integer :: particle

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    if (size(pdgs) < particle_count .or. &
        size(particle_is_final) < particle_count .or. &
        size(target_kinds) < particle_count .or. &
        size(target_ids) < particle_count) then
      call fail_factorized_phase_space( &
           'a block-local layout has inconsistent bounds')
    end if
    do particle = 1, particle_count
      if (target_kinds(particle) < factorized_no_target .or. &
          target_kinds(particle) > factorized_block_target) then
        call fail_factorized_phase_space( &
             'a block-local target kind is invalid')
      end if
      if (particle_is_final(particle)) then
        if (target_kinds(particle) == factorized_no_target .or. &
            target_ids(particle) < 1) then
          call fail_factorized_phase_space( &
               'a final block-local particle has no target')
        end if
      else if (target_kinds(particle) /= factorized_no_target .or. &
               target_ids(particle) /= 0) then
        call fail_factorized_phase_space( &
             'an incoming block-local particle has a target')
      end if
      if (target_kinds(particle) == factorized_block_target) then
        call validate_block(target_ids(particle))
        if (target_ids(particle) == block) then
          call fail_factorized_phase_space( &
               'a block-local particle targets its own block')
        end if
      end if
    end do

    block_state(block)%slot(event_slot)%local_pdg = 0
    block_state(block)%slot(event_slot)%local_final = .false.
    block_state(block)%slot(event_slot)%local_target_kind = &
         factorized_no_target
    block_state(block)%slot(event_slot)%local_target_id = 0
    block_state(block)%slot(event_slot)%local_pdg(1:particle_count) = &
         pdgs(1:particle_count)
    block_state(block)%slot(event_slot)% &
         local_final(1:particle_count) = &
         particle_is_final(1:particle_count)
    block_state(block)%slot(event_slot)% &
         local_target_kind(1:particle_count) = &
         target_kinds(1:particle_count)
    block_state(block)%slot(event_slot)% &
         local_target_id(1:particle_count) = &
         target_ids(1:particle_count)
    block_state(block)%slot(event_slot)%local_layout_valid = .true.
    if (block_state(block)%slot(event_slot)%is_valid .and. &
        .not. block_state(block)%slot(event_slot)%matrix_layout_valid .and. &
        block_state(block)%slot(event_slot)%particle_count <= &
        particle_count) then
      call capture_matrix_layout( &
           event_slot, block, &
           block_state(block)%slot(event_slot)%particle_count)
    end if
  end subroutine store_factorized_local_layout


  subroutine capture_matrix_layout(event_slot, block, particle_count)
    integer, intent(in) :: event_slot, block, particle_count

    block_state(block)%slot(event_slot)%matrix_count = particle_count
    block_state(block)%slot(event_slot)%matrix_pdg = 0
    block_state(block)%slot(event_slot)%matrix_final = .false.
    block_state(block)%slot(event_slot)% &
         matrix_pdg(1:particle_count) = &
         block_state(block)%slot(event_slot)% &
         local_pdg(1:particle_count)
    block_state(block)%slot(event_slot)% &
         matrix_final(1:particle_count) = &
         block_state(block)%slot(event_slot)% &
         local_final(1:particle_count)
    block_state(block)%slot(event_slot)%matrix_layout_valid = .true.
  end subroutine capture_matrix_layout


  subroutine fetch_factorized_local_layout( &
       event_slot, block, particle_count, pdgs, particle_is_final, &
       target_kinds, target_ids, available)
    integer, intent(in) :: event_slot, block, particle_count
    integer, intent(out) :: pdgs(particle_count)
    logical, intent(out) :: particle_is_final(particle_count)
    integer, intent(out) :: target_kinds(particle_count)
    integer, intent(out) :: target_ids(particle_count)
    logical, intent(out) :: available

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    available = block_state(block)%slot(event_slot)% &
         local_layout_valid .and. &
         block_state(block)%slot(event_slot)%local_count == particle_count
    pdgs = 0
    particle_is_final = .false.
    target_kinds = factorized_no_target
    target_ids = 0
    if (available) then
      pdgs = block_state(block)%slot(event_slot)% &
           local_pdg(1:particle_count)
      particle_is_final = &
           block_state(block)%slot(event_slot)% &
           local_final(1:particle_count)
      target_kinds = &
           block_state(block)%slot(event_slot)% &
           local_target_kind(1:particle_count)
      target_ids = block_state(block)%slot(event_slot)% &
           local_target_id(1:particle_count)
    end if
  end subroutine fetch_factorized_local_layout


  subroutine store_factorized_radiation_state(event_slot, block, state)
    integer, intent(in) :: event_slot, block
    type(factorized_radiation_state), intent(in) :: state

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    block_state(block)%slot(event_slot)%radiation = state
    block_state(block)%slot(event_slot)%radiation_valid = .true.
  end subroutine store_factorized_radiation_state


  subroutine fetch_factorized_radiation_state(event_slot, block, state, &
                                               available)
    integer, intent(in) :: event_slot, block
    type(factorized_radiation_state), intent(out) :: state
    logical, intent(out) :: available

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    available = block_state(block)%slot(event_slot)%radiation_valid
    state = factorized_radiation_state()
    if (available) state = block_state(block)%slot(event_slot)%radiation
  end subroutine fetch_factorized_radiation_state


  subroutine scale_factorized_radiation_jacobians(weight)
    double precision, intent(in) :: weight
    integer :: block, event_slot

    call ensure_storage()
    do event_slot = soft_counterevent, real_event
      do block = 0, nexternal
        if (block_state(block)%slot(event_slot)%radiation_valid) then
          block_state(block)%slot(event_slot)%radiation%jacobian = &
               block_state(block)%slot(event_slot)% &
               radiation%jacobian*weight
        end if
      end do
    end do
  end subroutine scale_factorized_radiation_jacobians


  subroutine store_factorized_base_measure(block, measure)
    integer, intent(in) :: block
    type(factorized_measure_state), intent(in) :: measure

    call ensure_storage()
    call validate_block(block)
    call validate_measure(measure)
    block_state(block)%base = measure
    block_state(block)%base_valid = .true.
  end subroutine store_factorized_base_measure


  subroutine multiply_factorized_base_measure(block, measure)
    integer, intent(in) :: block
    type(factorized_measure_state), intent(in) :: measure

    call ensure_storage()
    call validate_block(block)
    call validate_measure(measure)
    if (.not. block_state(block)%base_valid) then
      block_state(block)%base = factorized_measure_state()
      block_state(block)%base_valid = .true.
    end if
    call multiply_measure(block_state(block)%base, measure)
  end subroutine multiply_factorized_base_measure


  subroutine fetch_factorized_base_measure(block, measure, available)
    integer, intent(in) :: block
    type(factorized_measure_state), intent(out) :: measure
    logical, intent(out) :: available

    call ensure_storage()
    call validate_block(block)
    available = block_state(block)%base_valid
    measure = factorized_measure_state()
    if (available) measure = block_state(block)%base
  end subroutine fetch_factorized_base_measure


  subroutine store_factorized_event_measure(event_slot, block, measure)
    integer, intent(in) :: event_slot, block
    type(factorized_measure_state), intent(in) :: measure

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    call validate_measure(measure)
    block_state(block)%slot(event_slot)%measure = measure
    block_state(block)%slot(event_slot)%measure_valid = .true.
  end subroutine store_factorized_event_measure


  subroutine multiply_factorized_event_measure(event_slot, block, measure)
    integer, intent(in) :: event_slot, block
    type(factorized_measure_state), intent(in) :: measure

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    call validate_measure(measure)
    if (.not. block_state(block)%slot(event_slot)%measure_valid) then
      block_state(block)%slot(event_slot)%measure = &
           factorized_measure_state()
      block_state(block)%slot(event_slot)%measure_valid = .true.
    end if
    call multiply_measure(block_state(block)%slot(event_slot)%measure, measure)
  end subroutine multiply_factorized_event_measure


  subroutine fetch_factorized_event_measure(event_slot, block, measure, &
                                             available)
    integer, intent(in) :: event_slot, block
    type(factorized_measure_state), intent(out) :: measure
    logical, intent(out) :: available

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    available = block_state(block)%slot(event_slot)%measure_valid
    measure = factorized_measure_state()
    if (available) measure = block_state(block)%slot(event_slot)%measure
  end subroutine fetch_factorized_event_measure


  subroutine compose_factorized_base_measure(jacobian, phase_space_weight, &
                                              available)
    double precision, intent(out) :: jacobian, phase_space_weight
    logical, intent(out) :: available
    integer :: block

    call ensure_storage()
    jacobian = 1d0
    phase_space_weight = 1d0
    available = block_state(0)%base_valid
    if (.not. available) return
    do block = 0, nexternal
      if (.not. block_state(block)%base_valid) cycle
      jacobian = jacobian*block_state(block)%base%jacobian
      phase_space_weight = phase_space_weight* &
           block_state(block)%base%phase_space_weight
    end do
  end subroutine compose_factorized_base_measure


  subroutine compose_factorized_event_measure(event_slot, jacobian, &
                                               phase_space_weight, available)
    integer, intent(in) :: event_slot
    double precision, intent(out) :: jacobian, phase_space_weight
    logical, intent(out) :: available
    integer :: event_slots(0:nexternal)

    call validate_event_and_block(event_slot, 0)
    event_slots = event_slot
    call compose_factorized_tuple_measure( &
         event_slots, jacobian, phase_space_weight, available)
  end subroutine compose_factorized_event_measure


  subroutine compose_factorized_tuple_measure(event_slots, jacobian, &
                                               phase_space_weight, available)
    integer, intent(in) :: event_slots(0:)
    double precision, intent(out) :: jacobian, phase_space_weight
    logical, intent(out) :: available
    integer :: block, selected_slot

    call ensure_storage()
    if (ubound(event_slots, 1) < nexternal) then
      call fail_factorized_phase_space( &
           'an event-slot tuple has the wrong size')
    end if
    call compose_factorized_base_measure( &
         jacobian, phase_space_weight, available)
    if (.not. available) return
    available = .false.
    do block = 0, nexternal
      selected_slot = event_slots(block)
      call validate_event_and_block(selected_slot, block)
      if (block_state(block)%slot(selected_slot)%measure_valid) then
        available = .true.
        jacobian = jacobian* &
             block_state(block)%slot(selected_slot)%measure%jacobian
        phase_space_weight = phase_space_weight* &
             block_state(block)%slot(selected_slot)% &
             measure%phase_space_weight
      end if
    end do
    if (.not. available) return
  end subroutine compose_factorized_tuple_measure


  subroutine scale_factorized_event_measures(weight, block)
    double precision, intent(in) :: weight
    integer, intent(in), optional :: block
    type(factorized_measure_state) :: incoming_measure
    integer :: event_slot, target_block

    if (weight < 0d0) then
      call fail_factorized_phase_space('an incoming weight is negative')
    end if
    target_block = 0
    if (present(block)) target_block = block
    call validate_block(target_block)
    incoming_measure = factorized_measure_state()
    incoming_measure%jacobian = weight
    call ensure_storage()
    do event_slot = soft_counterevent, real_event
      if (block_state(target_block)%slot(event_slot)%measure_valid) then
        call multiply_factorized_event_measure( &
             event_slot, target_block, incoming_measure)
      end if
    end do
  end subroutine scale_factorized_event_measures


  subroutine multiply_measure(target, factor)
    type(factorized_measure_state), intent(inout) :: target
    type(factorized_measure_state), intent(in) :: factor

    target%jacobian = target%jacobian*factor%jacobian
    target%phase_space_weight = target%phase_space_weight* &
         factor%phase_space_weight
  end subroutine multiply_measure


  subroutine ensure_storage()
    integer :: block

    call validate_process_dimensions()
    if (allocated(block_state)) return
    allocate(block_state(0:nexternal))
    do block = 0, nexternal
      call initialize_block_state(block_state(block), block)
    end do
  end subroutine ensure_storage


  subroutine initialize_block_state(state, block)
    type(factorized_block_state), intent(inout) :: state
    integer, intent(in) :: block
    integer :: event_slot

    state%block = block
    state%initialized = .true.
    state%base = factorized_measure_state()
    state%base_valid = .false.
    do event_slot = soft_counterevent, real_event
      call initialize_slot_state(state%slot(event_slot))
    end do
  end subroutine initialize_block_state


  subroutine initialize_slot_state(state)
    type(factorized_slot_state), intent(inout) :: state

    if (.not. allocated(state%momenta)) &
         allocate(state%momenta(0:3, nexternal))
    if (.not. allocated(state%matrix_pdg)) &
         allocate(state%matrix_pdg(nexternal))
    if (.not. allocated(state%matrix_final)) &
         allocate(state%matrix_final(nexternal))
    if (.not. allocated(state%embedded)) &
         allocate(state%embedded(0:3, nexternal))
    if (.not. allocated(state%kernel)) &
         allocate(state%kernel(0:3, nexternal))
    if (.not. allocated(state%local)) &
         allocate(state%local(0:3, nexternal))
    if (.not. allocated(state%local_pdg)) &
         allocate(state%local_pdg(nexternal))
    if (.not. allocated(state%local_final)) &
         allocate(state%local_final(nexternal))
    if (.not. allocated(state%local_target_kind)) &
         allocate(state%local_target_kind(nexternal))
    if (.not. allocated(state%local_target_id)) &
         allocate(state%local_target_id(nexternal))

    state%particle_count = 0
    state%is_valid = .false.
    state%momentum_revision = 0_8
    state%matrix_count = 0
    state%matrix_layout_valid = .false.
    state%embedded_count = 0
    state%embedded_valid = .false.
    state%kernel_count = 0
    state%kernel_valid = .false.
    state%local_count = 0
    state%local_valid = .false.
    state%local_layout_valid = .false.
    state%radiation = factorized_radiation_state()
    state%radiation_valid = .false.
    state%measure = factorized_measure_state()
    state%measure_valid = .false.
    state%momenta = 0d0
    state%matrix_pdg = 0
    state%matrix_final = .false.
    state%embedded = 0d0
    state%kernel = 0d0
    state%local = 0d0
    state%local_pdg = 0
    state%local_final = .false.
    state%local_target_kind = factorized_no_target
    state%local_target_id = 0
  end subroutine initialize_slot_state


  subroutine validate_indices(event_slot, block, particle_count)
    integer, intent(in) :: event_slot, block, particle_count

    call validate_event_and_block(event_slot, block)
    if (particle_count < 1 .or. particle_count > nexternal) then
      call fail_factorized_phase_space( &
           'a block particle count is out of range')
    end if
  end subroutine validate_indices


  subroutine validate_event_and_block(event_slot, block)
    integer, intent(in) :: event_slot, block

    if (event_slot < soft_counterevent .or. event_slot > real_event) then
      call fail_factorized_phase_space('an event slot is out of range')
    end if
    call validate_block(block)
  end subroutine validate_event_and_block


  subroutine validate_block(block)
    integer, intent(in) :: block

    if (block < 0 .or. block > nexternal) then
      call fail_factorized_phase_space('a phase-space block is out of range')
    end if
  end subroutine validate_block


  subroutine validate_measure(measure)
    type(factorized_measure_state), intent(in) :: measure

    if (measure%jacobian < 0d0 .or. measure%phase_space_weight < 0d0) then
      call fail_factorized_phase_space('a block measure is negative')
    end if
  end subroutine validate_measure


  subroutine fail_factorized_phase_space(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in factorized_phase_space: '//trim(message)
    stop 1
  end subroutine fail_factorized_phase_space

end module factorized_phase_space


subroutine get_factorized_block_momenta(event_slot, block, &
                                        particle_count, momenta)
  use factorized_phase_space, only: fetch_factorized_matrix_momenta
  implicit none
  integer, intent(in) :: event_slot, block, particle_count
  double precision, intent(out) :: momenta(0:3, particle_count)
  logical :: available

  call fetch_factorized_matrix_momenta(event_slot, block, particle_count, &
                                       momenta, available)
  if (.not. available) then
    write (*, '(a,i0,a,i0)') &
         'ERROR: boosted factorized momenta are unavailable for block ', &
         block, ' in event slot ', event_slot
    stop 1
  end if
end subroutine get_factorized_block_momenta


subroutine get_factorized_block_momenta_ordered( &
     event_slot, block, particle_count, expected_pdgs, expected_final, &
     momenta)
  use factorized_phase_space, only: &
       fetch_factorized_ordered_matrix_momenta
  implicit none
  integer, intent(in) :: event_slot, block, particle_count
  integer, intent(in) :: expected_pdgs(particle_count)
  integer, intent(in) :: expected_final(particle_count)
  double precision, intent(out) :: momenta(0:3, particle_count)
  logical :: available

  call fetch_factorized_ordered_matrix_momenta( &
       event_slot, block, particle_count, expected_pdgs, expected_final, &
       momenta, available)
  if (.not. available) then
    write (*, '(a,i0,a,i0)') &
         'ERROR: ordered boosted momenta are unavailable for block ', &
         block, ' in event slot ', event_slot
    stop 1
  end if
end subroutine get_factorized_block_momenta_ordered
