module factorized_phase_space
  use process_dimensions, only: nexternal, validate_process_dimensions
  use fnlo_process_common, only: soft_counterevent, real_event
  implicit none
  private

  double precision, allocatable, save :: block_momenta(:, :, :, :)
  integer, allocatable, save :: block_particle_count(:, :)
  logical, allocatable, save :: block_is_valid(:, :)
  integer(kind=8), allocatable, save :: block_momentum_revision(:, :)
  integer(kind=8), save :: next_block_momentum_revision = 0_8

  ! Embedding storage is deliberately distinct from the matrix-element
  ! block cache above.  In particular, a soft projection can have a
  ! real-context particle layout while the Born density matrix in slot zero
  ! must retain its Born layout.  Event materialization consumes this cache;
  ! matrix elements never do.
  double precision, allocatable, save :: embedded_momenta(:, :, :, :)
  integer, allocatable, save :: embedded_particle_count(:, :)
  logical, allocatable, save :: embedded_is_valid(:, :)

  ! Matrix elements and subtraction kernels deliberately use different
  ! representations of the same block.  The former needs the block boosted
  ! into the event frame, whereas the latter is most naturally evaluated in
  ! the frame in which its three FKS variables were generated.
  double precision, allocatable, save :: kernel_momenta(:, :, :, :)
  integer, allocatable, save :: kernel_particle_count(:, :)
  logical, allocatable, save :: kernel_is_valid(:, :)

  ! Tuple realization starts from an immutable block-local representation.
  ! For a decay this is the parent rest frame; for production it is the
  ! generated partonic event frame.  It must remain distinct from both the
  ! boosted matrix-element cache and the subtraction-kernel cache: selecting
  ! a different ancestor term changes the boost of every descendant without
  ! changing the descendant's independently generated local configuration.
  double precision, allocatable, save :: local_momenta(:, :, :, :)
  integer, allocatable, save :: local_particle_count(:, :)
  logical, allocatable, save :: local_is_valid(:, :)
  integer, parameter, public :: factorized_no_target = 0
  integer, parameter, public :: factorized_visible_target = 1
  integer, parameter, public :: factorized_block_target = 2
  integer, allocatable, save :: local_pdg(:, :, :)
  logical, allocatable, save :: local_particle_is_final(:, :, :)
  integer, allocatable, save :: local_target_kind(:, :, :)
  integer, allocatable, save :: local_target_id(:, :, :)
  logical, allocatable, save :: local_layout_is_valid(:, :)

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

  type(factorized_radiation_state), allocatable, save :: block_radiation(:, :)
  logical, allocatable, save :: block_radiation_is_valid(:, :)

  ! The underlying-Born measure is a product of one independently generated
  ! production block and one factor for every decay node.  Event measures
  ! contain only the event-dependent radiation map, flux and incoming MINT
  ! weight.  Keeping the two layers separate lets a block be regenerated or
  ! replaced without passing a cumulative Jacobian through the decay tree.
  type, public :: factorized_measure_state
    double precision :: jacobian = 1d0
    double precision :: phase_space_weight = 1d0
  end type factorized_measure_state

  type(factorized_measure_state), allocatable, save :: base_measure(:)
  logical, allocatable, save :: base_measure_is_valid(:)
  type(factorized_measure_state), allocatable, save :: event_measure(:, :)
  logical, allocatable, save :: event_measure_is_valid(:, :)

  ! Sequential generation of independent radiative blocks must not let the
  ! later block overwrite an earlier block's real/counterevent family.  A
  ! snapshot owns one physical block in every event slot and can therefore be
  ! overlaid on a common Born baseline without copying or pre-summing events.
  type, public :: factorized_block_snapshot
    integer :: block = -1
    logical :: initialized = .false.
    double precision, allocatable :: momenta(:, :, :)
    integer, allocatable :: particle_count(:)
    logical, allocatable :: is_valid(:)
    integer(kind=8), allocatable :: momentum_revision(:)
    double precision, allocatable :: embedded(:, :, :)
    integer, allocatable :: embedded_count(:)
    logical, allocatable :: embedded_valid(:)
    double precision, allocatable :: kernel(:, :, :)
    integer, allocatable :: kernel_count(:)
    logical, allocatable :: kernel_valid(:)
    double precision, allocatable :: local(:, :, :)
    integer, allocatable :: local_count(:)
    logical, allocatable :: local_valid(:)
    integer, allocatable :: local_pdg(:, :)
    logical, allocatable :: local_final(:, :)
    integer, allocatable :: local_target_kind(:, :)
    integer, allocatable :: local_target_id(:, :)
    logical, allocatable :: local_layout_valid(:)
    type(factorized_radiation_state), allocatable :: radiation(:)
    logical, allocatable :: radiation_valid(:)
    type(factorized_measure_state) :: base
    logical :: base_valid = .false.
    type(factorized_measure_state), allocatable :: measure(:)
    logical, allocatable :: measure_valid(:)
  end type factorized_block_snapshot

  public :: reset_factorized_phase_space
  public :: store_factorized_block_momenta
  public :: fetch_factorized_block_momenta
  public :: fetch_factorized_matrix_momenta
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
    call ensure_storage()
    block_momenta = 0d0
    block_particle_count = 0
    block_is_valid = .false.
    block_momentum_revision = 0_8
    embedded_momenta = 0d0
    embedded_particle_count = 0
    embedded_is_valid = .false.
    kernel_momenta = 0d0
    kernel_particle_count = 0
    kernel_is_valid = .false.
    local_momenta = 0d0
    local_particle_count = 0
    local_is_valid = .false.
    local_pdg = 0
    local_particle_is_final = .false.
    local_target_kind = factorized_no_target
    local_target_id = 0
    local_layout_is_valid = .false.
    block_radiation = factorized_radiation_state()
    block_radiation_is_valid = .false.
    base_measure = factorized_measure_state()
    base_measure_is_valid = .false.
    event_measure = factorized_measure_state()
    event_measure_is_valid = .false.
  end subroutine reset_factorized_phase_space


  subroutine capture_factorized_block(block, snapshot)
    integer, intent(in) :: block
    type(factorized_block_snapshot), intent(out) :: snapshot

    call ensure_storage()
    call validate_block(block)
    allocate(snapshot%momenta(0:3, nexternal, &
                              soft_counterevent:real_event))
    allocate(snapshot%particle_count(soft_counterevent:real_event))
    allocate(snapshot%is_valid(soft_counterevent:real_event))
    allocate(snapshot%momentum_revision(soft_counterevent:real_event))
    allocate(snapshot%embedded(0:3, nexternal, &
                               soft_counterevent:real_event))
    allocate(snapshot%embedded_count(soft_counterevent:real_event))
    allocate(snapshot%embedded_valid(soft_counterevent:real_event))
    allocate(snapshot%kernel(0:3, nexternal, &
                             soft_counterevent:real_event))
    allocate(snapshot%kernel_count(soft_counterevent:real_event))
    allocate(snapshot%kernel_valid(soft_counterevent:real_event))
    allocate(snapshot%local(0:3, nexternal, &
                            soft_counterevent:real_event))
    allocate(snapshot%local_count(soft_counterevent:real_event))
    allocate(snapshot%local_valid(soft_counterevent:real_event))
    allocate(snapshot%local_pdg(nexternal, &
                                soft_counterevent:real_event))
    allocate(snapshot%local_final(nexternal, &
                                  soft_counterevent:real_event))
    allocate(snapshot%local_target_kind(nexternal, &
                                        soft_counterevent:real_event))
    allocate(snapshot%local_target_id(nexternal, &
                                      soft_counterevent:real_event))
    allocate(snapshot%local_layout_valid(soft_counterevent:real_event))
    allocate(snapshot%radiation(soft_counterevent:real_event))
    allocate(snapshot%radiation_valid(soft_counterevent:real_event))
    allocate(snapshot%measure(soft_counterevent:real_event))
    allocate(snapshot%measure_valid(soft_counterevent:real_event))

    snapshot%block = block
    snapshot%momenta = block_momenta(:, :, block, :)
    snapshot%particle_count = block_particle_count(block, :)
    snapshot%is_valid = block_is_valid(block, :)
    snapshot%momentum_revision = block_momentum_revision(block, :)
    snapshot%embedded = embedded_momenta(:, :, block, :)
    snapshot%embedded_count = embedded_particle_count(block, :)
    snapshot%embedded_valid = embedded_is_valid(block, :)
    snapshot%kernel = kernel_momenta(:, :, block, :)
    snapshot%kernel_count = kernel_particle_count(block, :)
    snapshot%kernel_valid = kernel_is_valid(block, :)
    snapshot%local = local_momenta(:, :, block, :)
    snapshot%local_count = local_particle_count(block, :)
    snapshot%local_valid = local_is_valid(block, :)
    snapshot%local_pdg = local_pdg(:, block, :)
    snapshot%local_final = local_particle_is_final(:, block, :)
    snapshot%local_target_kind = local_target_kind(:, block, :)
    snapshot%local_target_id = local_target_id(:, block, :)
    snapshot%local_layout_valid = local_layout_is_valid(block, :)
    snapshot%radiation = block_radiation(block, :)
    snapshot%radiation_valid = block_radiation_is_valid(block, :)
    snapshot%base = base_measure(block)
    snapshot%base_valid = base_measure_is_valid(block)
    snapshot%measure = event_measure(block, :)
    snapshot%measure_valid = event_measure_is_valid(block, :)
    snapshot%initialized = .true.
  end subroutine capture_factorized_block


  subroutine restore_factorized_block(snapshot)
    type(factorized_block_snapshot), intent(in) :: snapshot
    integer :: block

    call ensure_storage()
    if (.not. snapshot%initialized) then
      call fail_factorized_phase_space('a block snapshot is uninitialized')
    end if
    block = snapshot%block
    call validate_block(block)
    if (.not. allocated(snapshot%momenta) .or. &
        size(snapshot%momenta, 1) /= 4 .or. &
        size(snapshot%momenta, 2) /= nexternal .or. &
        size(snapshot%momenta, 3) /= &
        real_event - soft_counterevent + 1) then
      call fail_factorized_phase_space('a block snapshot has the wrong shape')
    end if

    block_momenta(:, :, block, :) = snapshot%momenta
    block_particle_count(block, :) = snapshot%particle_count
    block_is_valid(block, :) = snapshot%is_valid
    block_momentum_revision(block, :) = snapshot%momentum_revision
    embedded_momenta(:, :, block, :) = snapshot%embedded
    embedded_particle_count(block, :) = snapshot%embedded_count
    embedded_is_valid(block, :) = snapshot%embedded_valid
    kernel_momenta(:, :, block, :) = snapshot%kernel
    kernel_particle_count(block, :) = snapshot%kernel_count
    kernel_is_valid(block, :) = snapshot%kernel_valid
    local_momenta(:, :, block, :) = snapshot%local
    local_particle_count(block, :) = snapshot%local_count
    local_is_valid(block, :) = snapshot%local_valid
    local_pdg(:, block, :) = snapshot%local_pdg
    local_particle_is_final(:, block, :) = snapshot%local_final
    local_target_kind(:, block, :) = snapshot%local_target_kind
    local_target_id(:, block, :) = snapshot%local_target_id
    local_layout_is_valid(block, :) = snapshot%local_layout_valid
    block_radiation(block, :) = snapshot%radiation
    block_radiation_is_valid(block, :) = snapshot%radiation_valid
    base_measure(block) = snapshot%base
    base_measure_is_valid(block) = snapshot%base_valid
    event_measure(block, :) = snapshot%measure
    event_measure_is_valid(block, :) = snapshot%measure_valid
  end subroutine restore_factorized_block


  subroutine store_factorized_block_momenta(event_slot, block, &
                                             particle_count, momenta)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(in) :: momenta(0:, :)

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    if (size(momenta, 1) < 4 .or. size(momenta, 2) < particle_count) then
      call fail_factorized_phase_space( &
           'a block momentum array has inconsistent bounds')
    end if
    block_momenta(:, :, block, event_slot) = 0d0
    block_momenta(:, 1:particle_count, block, event_slot) = &
         momenta(0:3, 1:particle_count)
    block_particle_count(block, event_slot) = particle_count
    block_is_valid(block, event_slot) = .true.
    next_block_momentum_revision = next_block_momentum_revision + 1_8
    if (next_block_momentum_revision <= 0_8) then
      ! A wrap is fantastically unlikely, but resetting all identities is
      ! safer than allowing a stale density-matrix cache entry to match.
      next_block_momentum_revision = 1_8
      block_momentum_revision = 0_8
    end if
    block_momentum_revision(block, event_slot) = &
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
    available = block_is_valid(block, event_slot) .and. &
         block_particle_count(block, event_slot) == particle_count
    momenta = 0d0
    if (available) then
      momenta = block_momenta(:, 1:particle_count, block, event_slot)
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
    available = block_is_valid(block, event_slot) .and. &
         block_particle_count(block, event_slot) >= particle_count
    momenta = 0d0
    if (available) then
      momenta = block_momenta(:, 1:particle_count, block, event_slot)
    end if
  end subroutine fetch_factorized_matrix_momenta


  integer(kind=8) function factorized_block_momentum_revision( &
       event_slot, block)
    integer, intent(in) :: event_slot, block

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    if (block_is_valid(block, event_slot)) then
      factorized_block_momentum_revision = &
           block_momentum_revision(block, event_slot)
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
    embedded_momenta(:, :, block, event_slot) = 0d0
    embedded_momenta(:, 1:particle_count, block, event_slot) = &
         momenta(0:3, 1:particle_count)
    embedded_particle_count(block, event_slot) = particle_count
    embedded_is_valid(block, event_slot) = .true.
  end subroutine store_factorized_embedded_momenta


  subroutine fetch_factorized_embedded_momenta(event_slot, block, &
                                                particle_count, momenta, &
                                                available)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(out) :: momenta(0:3, particle_count)
    logical, intent(out) :: available

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    available = embedded_is_valid(block, event_slot) .and. &
         embedded_particle_count(block, event_slot) == particle_count
    momenta = 0d0
    if (available) then
      momenta = embedded_momenta(:, 1:particle_count, block, event_slot)
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
    kernel_momenta(:, :, block, event_slot) = 0d0
    kernel_momenta(:, 1:particle_count, block, event_slot) = &
         momenta(0:3, 1:particle_count)
    kernel_particle_count(block, event_slot) = particle_count
    kernel_is_valid(block, event_slot) = .true.
  end subroutine store_factorized_kernel_momenta


  subroutine fetch_factorized_kernel_momenta(event_slot, block, &
                                              particle_count, momenta, &
                                              available)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(out) :: momenta(0:3, particle_count)
    logical, intent(out) :: available

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    available = kernel_is_valid(block, event_slot) .and. &
         kernel_particle_count(block, event_slot) == particle_count
    momenta = 0d0
    if (available) then
      momenta = kernel_momenta(:, 1:particle_count, block, event_slot)
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
    local_momenta(:, :, block, event_slot) = 0d0
    local_momenta(:, 1:particle_count, block, event_slot) = &
         momenta(0:3, 1:particle_count)
    local_particle_count(block, event_slot) = particle_count
    local_is_valid(block, event_slot) = .true.
  end subroutine store_factorized_local_momenta


  subroutine fetch_factorized_local_momenta(event_slot, block, &
                                             particle_count, momenta, &
                                             available)
    integer, intent(in) :: event_slot, block, particle_count
    double precision, intent(out) :: momenta(0:3, particle_count)
    logical, intent(out) :: available

    call ensure_storage()
    call validate_indices(event_slot, block, particle_count)
    available = local_is_valid(block, event_slot) .and. &
         local_particle_count(block, event_slot) == particle_count
    momenta = 0d0
    if (available) then
      momenta = local_momenta(:, 1:particle_count, block, event_slot)
    end if
  end subroutine fetch_factorized_local_momenta


  integer function factorized_local_count(event_slot, block)
    integer, intent(in) :: event_slot, block

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    if (local_is_valid(block, event_slot)) then
      factorized_local_count = local_particle_count(block, event_slot)
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

    local_pdg(:, block, event_slot) = 0
    local_particle_is_final(:, block, event_slot) = .false.
    local_target_kind(:, block, event_slot) = factorized_no_target
    local_target_id(:, block, event_slot) = 0
    local_pdg(1:particle_count, block, event_slot) = pdgs(1:particle_count)
    local_particle_is_final(1:particle_count, block, event_slot) = &
         particle_is_final(1:particle_count)
    local_target_kind(1:particle_count, block, event_slot) = &
         target_kinds(1:particle_count)
    local_target_id(1:particle_count, block, event_slot) = &
         target_ids(1:particle_count)
    local_layout_is_valid(block, event_slot) = .true.
  end subroutine store_factorized_local_layout


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
    available = local_layout_is_valid(block, event_slot) .and. &
         local_particle_count(block, event_slot) == particle_count
    pdgs = 0
    particle_is_final = .false.
    target_kinds = factorized_no_target
    target_ids = 0
    if (available) then
      pdgs = local_pdg(1:particle_count, block, event_slot)
      particle_is_final = &
           local_particle_is_final(1:particle_count, block, event_slot)
      target_kinds = &
           local_target_kind(1:particle_count, block, event_slot)
      target_ids = local_target_id(1:particle_count, block, event_slot)
    end if
  end subroutine fetch_factorized_local_layout


  subroutine store_factorized_radiation_state(event_slot, block, state)
    integer, intent(in) :: event_slot, block
    type(factorized_radiation_state), intent(in) :: state

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    block_radiation(block, event_slot) = state
    block_radiation_is_valid(block, event_slot) = .true.
  end subroutine store_factorized_radiation_state


  subroutine fetch_factorized_radiation_state(event_slot, block, state, &
                                               available)
    integer, intent(in) :: event_slot, block
    type(factorized_radiation_state), intent(out) :: state
    logical, intent(out) :: available

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    available = block_radiation_is_valid(block, event_slot)
    state = factorized_radiation_state()
    if (available) state = block_radiation(block, event_slot)
  end subroutine fetch_factorized_radiation_state


  subroutine scale_factorized_radiation_jacobians(weight)
    double precision, intent(in) :: weight
    integer :: block, event_slot

    call ensure_storage()
    do event_slot = soft_counterevent, real_event
      do block = 0, nexternal
        if (block_radiation_is_valid(block, event_slot)) then
          block_radiation(block, event_slot)%jacobian = &
               block_radiation(block, event_slot)%jacobian*weight
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
    base_measure(block) = measure
    base_measure_is_valid(block) = .true.
  end subroutine store_factorized_base_measure


  subroutine multiply_factorized_base_measure(block, measure)
    integer, intent(in) :: block
    type(factorized_measure_state), intent(in) :: measure

    call ensure_storage()
    call validate_block(block)
    call validate_measure(measure)
    if (.not. base_measure_is_valid(block)) then
      base_measure(block) = factorized_measure_state()
      base_measure_is_valid(block) = .true.
    end if
    call multiply_measure(base_measure(block), measure)
  end subroutine multiply_factorized_base_measure


  subroutine fetch_factorized_base_measure(block, measure, available)
    integer, intent(in) :: block
    type(factorized_measure_state), intent(out) :: measure
    logical, intent(out) :: available

    call ensure_storage()
    call validate_block(block)
    available = base_measure_is_valid(block)
    measure = factorized_measure_state()
    if (available) measure = base_measure(block)
  end subroutine fetch_factorized_base_measure


  subroutine store_factorized_event_measure(event_slot, block, measure)
    integer, intent(in) :: event_slot, block
    type(factorized_measure_state), intent(in) :: measure

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    call validate_measure(measure)
    event_measure(block, event_slot) = measure
    event_measure_is_valid(block, event_slot) = .true.
  end subroutine store_factorized_event_measure


  subroutine multiply_factorized_event_measure(event_slot, block, measure)
    integer, intent(in) :: event_slot, block
    type(factorized_measure_state), intent(in) :: measure

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    call validate_measure(measure)
    if (.not. event_measure_is_valid(block, event_slot)) then
      event_measure(block, event_slot) = factorized_measure_state()
      event_measure_is_valid(block, event_slot) = .true.
    end if
    call multiply_measure(event_measure(block, event_slot), measure)
  end subroutine multiply_factorized_event_measure


  subroutine fetch_factorized_event_measure(event_slot, block, measure, &
                                             available)
    integer, intent(in) :: event_slot, block
    type(factorized_measure_state), intent(out) :: measure
    logical, intent(out) :: available

    call ensure_storage()
    call validate_event_and_block(event_slot, block)
    available = event_measure_is_valid(block, event_slot)
    measure = factorized_measure_state()
    if (available) measure = event_measure(block, event_slot)
  end subroutine fetch_factorized_event_measure


  subroutine compose_factorized_base_measure(jacobian, phase_space_weight, &
                                              available)
    double precision, intent(out) :: jacobian, phase_space_weight
    logical, intent(out) :: available
    integer :: block

    call ensure_storage()
    jacobian = 1d0
    phase_space_weight = 1d0
    available = base_measure_is_valid(0)
    if (.not. available) return
    do block = 0, nexternal
      if (.not. base_measure_is_valid(block)) cycle
      jacobian = jacobian*base_measure(block)%jacobian
      phase_space_weight = phase_space_weight* &
           base_measure(block)%phase_space_weight
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
      if (event_measure_is_valid(block, selected_slot)) then
        available = .true.
        jacobian = jacobian* &
             event_measure(block, selected_slot)%jacobian
        phase_space_weight = phase_space_weight* &
             event_measure(block, selected_slot)%phase_space_weight
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
      if (event_measure_is_valid(target_block, event_slot)) then
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
    call validate_process_dimensions()
    if (allocated(block_momenta)) return
    allocate(block_momenta(0:3, nexternal, 0:nexternal, &
                           soft_counterevent:real_event))
    allocate(block_particle_count(0:nexternal, &
                                  soft_counterevent:real_event))
    allocate(block_is_valid(0:nexternal, &
                            soft_counterevent:real_event))
    allocate(block_momentum_revision(0:nexternal, &
                                     soft_counterevent:real_event))
    allocate(embedded_momenta(0:3, nexternal, 0:nexternal, &
                              soft_counterevent:real_event))
    allocate(embedded_particle_count(0:nexternal, &
                                     soft_counterevent:real_event))
    allocate(embedded_is_valid(0:nexternal, &
                               soft_counterevent:real_event))
    allocate(kernel_momenta(0:3, nexternal, 0:nexternal, &
                            soft_counterevent:real_event))
    allocate(kernel_particle_count(0:nexternal, &
                                   soft_counterevent:real_event))
    allocate(kernel_is_valid(0:nexternal, &
                             soft_counterevent:real_event))
    allocate(local_momenta(0:3, nexternal, 0:nexternal, &
                           soft_counterevent:real_event))
    allocate(local_particle_count(0:nexternal, &
                                  soft_counterevent:real_event))
    allocate(local_is_valid(0:nexternal, &
                            soft_counterevent:real_event))
    allocate(local_pdg(nexternal, 0:nexternal, &
                       soft_counterevent:real_event))
    allocate(local_particle_is_final(nexternal, 0:nexternal, &
                                     soft_counterevent:real_event))
    allocate(local_target_kind(nexternal, 0:nexternal, &
                               soft_counterevent:real_event))
    allocate(local_target_id(nexternal, 0:nexternal, &
                             soft_counterevent:real_event))
    allocate(local_layout_is_valid(0:nexternal, &
                                   soft_counterevent:real_event))
    allocate(block_radiation(0:nexternal, &
                             soft_counterevent:real_event))
    allocate(block_radiation_is_valid(0:nexternal, &
                                      soft_counterevent:real_event))
    allocate(base_measure(0:nexternal))
    allocate(base_measure_is_valid(0:nexternal))
    allocate(event_measure(0:nexternal, &
                           soft_counterevent:real_event))
    allocate(event_measure_is_valid(0:nexternal, &
                                    soft_counterevent:real_event))
    block_momenta = 0d0
    block_particle_count = 0
    block_is_valid = .false.
    block_momentum_revision = 0_8
    embedded_momenta = 0d0
    embedded_particle_count = 0
    embedded_is_valid = .false.
    kernel_momenta = 0d0
    kernel_particle_count = 0
    kernel_is_valid = .false.
    local_momenta = 0d0
    local_particle_count = 0
    local_is_valid = .false.
    local_pdg = 0
    local_particle_is_final = .false.
    local_target_kind = factorized_no_target
    local_target_id = 0
    local_layout_is_valid = .false.
    block_radiation = factorized_radiation_state()
    block_radiation_is_valid = .false.
    base_measure = factorized_measure_state()
    base_measure_is_valid = .false.
    event_measure = factorized_measure_state()
    event_measure_is_valid = .false.
  end subroutine ensure_storage


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
