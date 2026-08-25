module orderstags_glob_module
  implicit none
  private

  integer, allocatable :: global_order_tags(:)
  integer :: n_orderstags = 0
  logical :: orderstags_initialized = .false.

  public :: initialize_orderstags_glob
  public :: validate_orderstags_glob
  public :: finalize_orderstags_glob
  public :: orderstags_glob_is_initialized
  public :: orderstags_glob_count
  public :: orderstags_glob_position

contains

  subroutine initialize_orderstags_glob()
    integer :: file_n_orderstags, j

    if (orderstags_initialized) then
      call validate_orderstags_glob()
      return
    end if

    open(unit=78, file='orderstags_glob.dat', status='old', err=101)
    go to 99
101 open(unit=78, file='../orderstags_glob.dat', status='old')

99  read(78, *) file_n_orderstags
    write (*, *) 'get_orderstags_glob_infos: n_orderstags=', &
         file_n_orderstags
    if (file_n_orderstags < 1) then
      write (*, *) 'ERROR, get_orderstags_glob_infos, invalid count', &
           file_n_orderstags
      stop 1
    end if

    allocate(global_order_tags(file_n_orderstags))
    read(78, *) (global_order_tags(j), j=1, file_n_orderstags)
    close(78)
    n_orderstags = file_n_orderstags
    orderstags_initialized = .true.

    write (*, *) 'get_orderstags_glob_infos: orderstags_glob', &
         (global_order_tags(j), j=1, n_orderstags)
    call validate_orderstags_glob()
  end subroutine initialize_orderstags_glob


  subroutine validate_orderstags_glob()
    if (.not. orderstags_initialized) then
      write (*, *) 'ERROR, orderstags_glob is not initialized'
      stop 1
    end if
    if (.not. allocated(global_order_tags)) then
      write (*, *) 'ERROR, orderstags_glob storage is not allocated'
      stop 1
    end if
    if (n_orderstags < 1 .or. size(global_order_tags) /= n_orderstags) then
      write (*, *) 'ERROR, orderstags_glob storage is inconsistent', &
           n_orderstags, size(global_order_tags)
      stop 1
    end if
  end subroutine validate_orderstags_glob


  subroutine finalize_orderstags_glob()
    if (allocated(global_order_tags)) deallocate(global_order_tags)
    n_orderstags = 0
    orderstags_initialized = .false.
  end subroutine finalize_orderstags_glob


  logical function orderstags_glob_is_initialized()
    orderstags_glob_is_initialized = orderstags_initialized
  end function orderstags_glob_is_initialized


  integer function orderstags_glob_count()
    call validate_orderstags_glob()
    orderstags_glob_count = n_orderstags
  end function orderstags_glob_count


  integer function orderstags_glob_position(tag)
    integer, intent(in) :: tag
    integer :: j

    call validate_orderstags_glob()
    do j = 1, n_orderstags
      if (global_order_tags(j) == tag) then
        orderstags_glob_position = j
        return
      end if
    end do
    write (*, *) 'ERROR, get_orderstags_glob_pos, not found', tag
    stop 1
  end function orderstags_glob_position

end module orderstags_glob_module
