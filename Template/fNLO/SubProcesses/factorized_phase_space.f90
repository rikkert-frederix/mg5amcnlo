module factorized_phase_space
  use process_dimensions, only: nexternal, validate_process_dimensions
  use fnlo_process_common, only: soft_counterevent, real_event
  implicit none
  private

  double precision, allocatable, save :: block_momenta(:, :, :, :)
  integer, allocatable, save :: block_particle_count(:, :)
  logical, allocatable, save :: block_is_valid(:, :)

  ! Matrix elements and subtraction kernels deliberately use different
  ! representations of the same block.  The former needs the block boosted
  ! into the event frame, whereas the latter is most naturally evaluated in
  ! the frame in which its three FKS variables were generated.
  double precision, allocatable, save :: kernel_momenta(:, :, :, :)
  integer, allocatable, save :: kernel_particle_count(:, :)
  logical, allocatable, save :: kernel_is_valid(:, :)

  type, public :: factorized_radiation_state
    double precision :: jacobian = -1d0
    double precision :: xi = -1d0
    double precision :: y = -2d0
    double precision :: xi_hat = -1d0
    double precision :: xi_max = -1d0
    double precision :: xi_norm = -1d0
    double precision :: fks_momentum(0:3) = -1d0
    double precision :: shat = -1d0
    double precision :: sqrt_shat = -1d0
    double precision :: y_to_cm = 0d0
  end type factorized_radiation_state

  type(factorized_radiation_state), allocatable, save :: block_radiation(:, :)
  logical, allocatable, save :: block_radiation_is_valid(:, :)

  public :: reset_factorized_phase_space
  public :: store_factorized_block_momenta
  public :: fetch_factorized_block_momenta
  public :: store_factorized_kernel_momenta
  public :: fetch_factorized_kernel_momenta
  public :: store_factorized_radiation_state
  public :: fetch_factorized_radiation_state
  public :: scale_factorized_radiation_jacobians

contains

  subroutine reset_factorized_phase_space()
    call ensure_storage()
    block_momenta = 0d0
    block_particle_count = 0
    block_is_valid = .false.
    kernel_momenta = 0d0
    kernel_particle_count = 0
    kernel_is_valid = .false.
    block_radiation = factorized_radiation_state()
    block_radiation_is_valid = .false.
  end subroutine reset_factorized_phase_space


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


  subroutine ensure_storage()
    call validate_process_dimensions()
    if (allocated(block_momenta)) return
    allocate(block_momenta(0:3, nexternal, 0:nexternal, &
                           soft_counterevent:real_event))
    allocate(block_particle_count(0:nexternal, &
                                  soft_counterevent:real_event))
    allocate(block_is_valid(0:nexternal, &
                            soft_counterevent:real_event))
    allocate(kernel_momenta(0:3, nexternal, 0:nexternal, &
                            soft_counterevent:real_event))
    allocate(kernel_particle_count(0:nexternal, &
                                   soft_counterevent:real_event))
    allocate(kernel_is_valid(0:nexternal, &
                             soft_counterevent:real_event))
    allocate(block_radiation(0:nexternal, &
                             soft_counterevent:real_event))
    allocate(block_radiation_is_valid(0:nexternal, &
                                      soft_counterevent:real_event))
    block_momenta = 0d0
    block_particle_count = 0
    block_is_valid = .false.
    kernel_momenta = 0d0
    kernel_particle_count = 0
    kernel_is_valid = .false.
    block_radiation = factorized_radiation_state()
    block_radiation_is_valid = .false.
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
    if (block < 0 .or. block > nexternal) then
      call fail_factorized_phase_space('a phase-space block is out of range')
    end if
  end subroutine validate_event_and_block


  subroutine fail_factorized_phase_space(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in factorized_phase_space: '//trim(message)
    stop 1
  end subroutine fail_factorized_phase_space

end module factorized_phase_space


subroutine get_factorized_block_momenta(event_slot, block, &
                                        particle_count, momenta)
  use factorized_phase_space, only: fetch_factorized_block_momenta
  implicit none
  integer, intent(in) :: event_slot, block, particle_count
  double precision, intent(out) :: momenta(0:3, particle_count)
  logical :: available

  call fetch_factorized_block_momenta(event_slot, block, particle_count, &
                                      momenta, available)
  if (.not. available) then
    write (*, '(a,i0,a,i0)') &
         'ERROR: boosted factorized momenta are unavailable for block ', &
         block, ' in event slot ', event_slot
    stop 1
  end if
end subroutine get_factorized_block_momenta
