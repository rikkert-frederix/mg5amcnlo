! Monte Carlo integration over integer-valued channels with importance
! sampling.  The integration state is owned by this module and grows on
! demand; the initial capacities retain the layout limits of the original
! fixed-form implementation.
module mc_integer_module
  implicit none
  private

  integer, parameter :: default_dimension_capacity = 50
  integer, parameter :: default_interval_capacity = 200
  double precision, parameter :: minimum_cell_size = 1d-3

  integer, save :: dimension_capacity = 0
  integer, save :: interval_capacity = 0
  integer, allocatable, save :: nintervals(:)
  integer, allocatable, save :: ncall(:, :)
  double precision, allocatable, save :: grid(:, :)
  double precision, allocatable, save :: acc(:, :)
  logical, allocatable, save :: firsttime(:)

  public :: initialize_mc_integer
  public :: finalize_mc_integer
  public :: get_mc_integer
  public :: reset_mc_grid
  public :: fill_mc_integer
  public :: empty_mc_integer
  public :: regrid_mc_integer

contains

  ! Allocate the module state, or grow it without discarding an existing
  ! integration grid.  Repeated calls that already fit are no-ops.
  subroutine initialize_mc_integer(required_dimensions, required_intervals)
    integer, intent(in), optional :: required_dimensions
    integer, intent(in), optional :: required_intervals
    integer :: requested_dimensions, requested_intervals
    integer :: old_dimension_capacity, old_interval_capacity
    integer, allocatable :: old_nintervals(:), old_ncall(:, :)
    double precision, allocatable :: old_grid(:, :), old_acc(:, :)
    logical, allocatable :: old_firsttime(:)

    requested_dimensions = default_dimension_capacity
    requested_intervals = default_interval_capacity
    if (present(required_dimensions)) then
      if (required_dimensions < 1) then
        write (*, *) 'ERROR in initialize_MC_integer: invalid dimension', &
             required_dimensions
        stop 1
      end if
      requested_dimensions = max(requested_dimensions, required_dimensions)
    end if
    if (present(required_intervals)) then
      if (required_intervals < 0) then
        write (*, *) 'ERROR in initialize_MC_integer: invalid interval count', &
             required_intervals
        stop 1
      end if
      requested_intervals = max(requested_intervals, required_intervals)
    end if

    if (.not. allocated(grid)) then
      dimension_capacity = requested_dimensions
      interval_capacity = requested_intervals
      allocate(grid(0:interval_capacity, dimension_capacity))
      allocate(acc(0:interval_capacity, dimension_capacity))
      allocate(ncall(0:interval_capacity, dimension_capacity))
      allocate(nintervals(dimension_capacity))
      allocate(firsttime(dimension_capacity))
      grid = 0d0
      acc = 0d0
      ncall = 0
      nintervals = 0
      firsttime = .true.
      return
    end if

    if (requested_dimensions <= dimension_capacity .and. &
        requested_intervals <= interval_capacity) return

    old_dimension_capacity = dimension_capacity
    old_interval_capacity = interval_capacity
    allocate(old_grid(0:old_interval_capacity, old_dimension_capacity))
    allocate(old_acc(0:old_interval_capacity, old_dimension_capacity))
    allocate(old_ncall(0:old_interval_capacity, old_dimension_capacity))
    allocate(old_nintervals(old_dimension_capacity))
    allocate(old_firsttime(old_dimension_capacity))
    old_grid = grid
    old_acc = acc
    old_ncall = ncall
    old_nintervals = nintervals
    old_firsttime = firsttime

    dimension_capacity = max(requested_dimensions, &
         2 * old_dimension_capacity)
    interval_capacity = max(requested_intervals, 2 * old_interval_capacity)
    deallocate(grid, acc, ncall, nintervals, firsttime)
    allocate(grid(0:interval_capacity, dimension_capacity))
    allocate(acc(0:interval_capacity, dimension_capacity))
    allocate(ncall(0:interval_capacity, dimension_capacity))
    allocate(nintervals(dimension_capacity))
    allocate(firsttime(dimension_capacity))
    grid = 0d0
    acc = 0d0
    ncall = 0
    nintervals = 0
    firsttime = .true.
    grid(0:old_interval_capacity, 1:old_dimension_capacity) = old_grid
    acc(0:old_interval_capacity, 1:old_dimension_capacity) = old_acc
    ncall(0:old_interval_capacity, 1:old_dimension_capacity) = old_ncall
    nintervals(1:old_dimension_capacity) = old_nintervals
    firsttime(1:old_dimension_capacity) = old_firsttime
    deallocate(old_grid, old_acc, old_ncall, old_nintervals, old_firsttime)
  end subroutine initialize_mc_integer


  subroutine finalize_mc_integer
    if (allocated(grid)) then
      deallocate(grid, acc, ncall, nintervals, firsttime)
    end if
    dimension_capacity = 0
    interval_capacity = 0
  end subroutine finalize_mc_integer


  subroutine get_mc_integer(this_dim, niint_thisd, iint, vol)
    integer, intent(in) :: this_dim, niint_thisd
    integer, intent(out) :: iint
    double precision, intent(out) :: vol
    integer :: i, io_status
    double precision :: rnd
    double precision, external :: ran2
    character(len=1) :: cdum
    character(len=3) :: action
    logical :: grid_file_open
    logical :: flat_grid
    common /to_readgrid/ flat_grid

    if (niint_thisd < 1) then
      write (*, *) 'ERROR in get_MC_integer: invalid interval count', &
           niint_thisd
      stop 1
    end if
    call initialize_mc_integer(this_dim, niint_thisd)

    if (firsttime(this_dim)) then
      firsttime(this_dim) = .false.
      nintervals(this_dim) = niint_thisd
      if (flat_grid) then
        call set_flat_grid(this_dim)
      else
        grid_file_open = .false.
        open(unit=52, file='grid.MC_integer', status='old', &
             iostat=io_status)
        if (io_status /= 0) then
          call warn_and_set_flat_grid(this_dim, grid_file_open)
        else
          grid_file_open = .true.
          do i = 1, this_dim - 1
            read(52, *, iostat=io_status) cdum
            if (io_status /= 0) exit
          end do
          if (io_status == 0) then
            read(52, *, iostat=io_status) action, &
                 (grid(i, this_dim), i=0, nintervals(this_dim))
          end if
          if (io_status == 0) then
            do i = this_dim + 1, dimension_capacity
              read(52, *, iostat=io_status) cdum
              if (io_status /= 0) exit
            end do
          end if
          if (io_status /= 0) then
            call warn_and_set_flat_grid(this_dim, grid_file_open)
          else
            close(52)
            grid_file_open = .false.
          end if
        end if
      end if
      acc(0:nintervals(this_dim), this_dim) = 0d0
      ncall(0:nintervals(this_dim), this_dim) = 0
    else if (nintervals(this_dim) /= niint_thisd) then
      write (*, *) 'ERROR in get_MC_integer: interval count changed for', &
           this_dim, nintervals(this_dim), niint_thisd
      stop 1
    end if

    rnd = ran2()
    iint = 0
    do while (rnd > grid(iint, this_dim))
      iint = iint + 1
    end do
    if (iint == 0 .or. iint > nintervals(this_dim)) then
      write (*, *) 'ERROR in get_MC_integer', iint, nintervals(this_dim), &
           (grid(i, this_dim), i=1, nintervals(this_dim))
      stop
    end if
    vol = grid(iint, this_dim) - grid(iint - 1, this_dim)
  end subroutine get_mc_integer


  subroutine warn_and_set_flat_grid(this_dim, grid_file_open)
    integer, intent(in) :: this_dim
    logical, intent(inout) :: grid_file_open

    if (grid_file_open) then
      close(52)
      grid_file_open = .false.
    end if
    write (*, *) 'WARNING: File "grid.MC_integer" not found. ' // &
         'Using flat grid to start for', this_dim
    call set_flat_grid(this_dim)
  end subroutine warn_and_set_flat_grid


  subroutine set_flat_grid(this_dim)
    integer, intent(in) :: this_dim
    integer :: i

    do i = 0, nintervals(this_dim)
      grid(i, this_dim) = dble(i) / nintervals(this_dim)
    end do
  end subroutine set_flat_grid


  subroutine reset_mc_grid
    integer :: this_dim

    call initialize_mc_integer
    do this_dim = 1, dimension_capacity
      if (nintervals(this_dim) /= 0) call set_flat_grid(this_dim)
      acc(0:nintervals(this_dim), this_dim) = 0d0
      ncall(0:nintervals(this_dim), this_dim) = 0
    end do
  end subroutine reset_mc_grid


  subroutine fill_mc_integer(this_dim, iint, f_abs)
    integer, intent(in) :: this_dim, iint
    double precision, intent(in) :: f_abs

    call initialize_mc_integer(this_dim, max(0, iint))
    if (nintervals(this_dim) == 0 .or. iint < 1 .or. &
        iint > nintervals(this_dim)) then
      write (*, *) 'ERROR in fill_MC_integer', this_dim, iint, &
           nintervals(this_dim)
      stop 1
    end if
    acc(iint, this_dim) = acc(iint, this_dim) + f_abs
    ncall(iint, this_dim) = ncall(iint, this_dim) + 1
  end subroutine fill_mc_integer


  subroutine empty_mc_integer
    integer :: this_dim

    call initialize_mc_integer
    do this_dim = 1, dimension_capacity
      acc(0:nintervals(this_dim), this_dim) = 0d0
      ncall(0:nintervals(this_dim), this_dim) = 0
    end do
  end subroutine empty_mc_integer


  subroutine regrid_mc_integer
    integer :: i, ib, this_dim, io_status
    character(len=101) :: buff

    call initialize_mc_integer
    do this_dim = 1, dimension_capacity
      if (nintervals(this_dim) == 0) cycle

      ncall(0, this_dim) = 0
      do i = 1, nintervals(this_dim)
        if (ncall(i, this_dim) /= 0) then
          acc(i, this_dim) = acc(i - 1, this_dim) + &
               acc(i, this_dim) / ncall(i, this_dim)
          ncall(0, this_dim) = ncall(0, this_dim) + ncall(i, this_dim)
        else
          acc(i, this_dim) = acc(i - 1, this_dim)
        end if
      end do
      if (ncall(0, this_dim) <= max(nintervals(this_dim), 10)) then
        acc(0:nintervals(this_dim), this_dim) = 0d0
        ncall(0:nintervals(this_dim), this_dim) = 0
        cycle
      end if

      if (acc(nintervals(this_dim), this_dim) /= 0d0) then
        do i = 0, nintervals(this_dim)
          grid(i, this_dim) = acc(i, this_dim) / &
               acc(nintervals(this_dim), this_dim)
        end do
      end if

      do i = 1, nintervals(this_dim)
        if (grid(i, this_dim) <= grid(i - 1, this_dim) + &
            minimum_cell_size) then
          grid(i, this_dim) = grid(i - 1, this_dim) + minimum_cell_size
        end if
      end do
      grid(nintervals(this_dim), this_dim) = 1d0
      do i = 1, nintervals(this_dim)
        if (grid(nintervals(this_dim) - i, this_dim) >= &
            grid(nintervals(this_dim) - i + 1, this_dim) - &
            minimum_cell_size) then
          grid(nintervals(this_dim) - i, this_dim) = &
               1d0 - dble(i) * minimum_cell_size
        else
          exit
        end if
      end do
    end do

    open(unit=52, file='grid.MC_integer', status='unknown', &
         iostat=io_status)
    if (io_status /= 0) then
      write (*, *) 'Cannot open "grid.MC_integer" file'
      stop
    end if
    do this_dim = 1, dimension_capacity
      write(52, *) 'AVE', &
           (grid(i, this_dim), i=0, nintervals(this_dim))
    end do
    close(52)

    do this_dim = 1, dimension_capacity
      if (nintervals(this_dim) == 0) cycle
      buff = ' '
      do i = 0, nintervals(this_dim)
        ib = 1 + int(grid(i, this_dim) * 100)
        write(buff(ib:ib), '(i1)') mod(i, 10)
      end do
      write (*, '(i3,a,a)') this_dim, ':  ', buff
    end do

    call empty_mc_integer
  end subroutine regrid_mc_integer

end module mc_integer_module
