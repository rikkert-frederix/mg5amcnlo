module fks_channel_map
  use process_dimensions, only: nexternal, nincoming, fks_configs, &
                                validate_process_dimensions
  use fks_metadata, only: fks_j_d, need_color_links_d, &
                          validate_fks_metadata
  implicit none
  private

  integer, allocatable, save :: channel_map(:, :)
  integer, allocatable, save :: born_process_map(:)
  logical, save :: channel_map_printed = .false.
  logical, save :: born_map_printed = .false.

  public :: fks_channel_count, fks_channel_configuration
  public :: print_fks_channel_map, get_born_fks_process

contains

  integer function fks_channel_count(category)
    implicit none
    integer, intent(in) :: category

    call ensure_channel_map()
    call validate_category(category)
    fks_channel_count = channel_map(category, 0)
  end function fks_channel_count


  integer function fks_channel_configuration(category, position)
    implicit none
    integer, intent(in) :: category, position

    call ensure_channel_map()
    call validate_category(category)
    if (position < 1 .or. position > channel_map(category, 0)) then
      call fail_channel_map('FKS channel position is out of range')
    end if
    fks_channel_configuration = channel_map(category, position)
  end function fks_channel_configuration


  subroutine print_fks_channel_map()
    implicit none

    call ensure_channel_map()
    if (channel_map_printed) return
    channel_map_printed = .true.
    write (*, *) 'initial-final FKS maps:'
    write (*, *) 0, ':', channel_map(0, :)
    write (*, *) 1, ':', channel_map(1, :)
    write (*, *) 2, ':', channel_map(2, :)
  end subroutine print_fks_channel_map


  subroutine get_born_fks_process(nfks_in, nfks_out)
    implicit none
    integer, intent(in) :: nfks_in
    integer, intent(out) :: nfks_out

    call ensure_born_map()
    if (nfks_in < 1 .or. nfks_in > fks_configs) then
      call fail_channel_map('FKS configuration is out of range')
    end if
    if (.not. born_map_printed) then
      born_map_printed = .true.
      write (*, *) 'Total number of FKS directories is', fks_configs
      write (*, *) 'For the Born we use nFKSprocesses:'
      write (*, *) born_process_map
    end if
    if (born_process_map(nfks_in) == 0) then
      write (*, *) 'Could not find the correct map to Born '// &
        'FKS configuration for the NLO FKS '// &
        'configuration', nfks_in
      stop 1
    end if
    nfks_out = born_process_map(nfks_in)
  end subroutine get_born_fks_process


  subroutine ensure_channel_map()
    implicit none
    integer :: ifks, emitter, emitter_position

    if (allocated(channel_map)) return
    call validate_process_dimensions()
    call validate_fks_metadata()
    allocate (channel_map(0:2, 0:fks_configs))
    channel_map = 0
    do ifks = 1, fks_configs
      emitter = fks_j_d(ifks)
      channel_map(0, 0) = channel_map(0, 0) + 1
      emitter_position = channel_map(0, 0)
      channel_map(0, emitter_position) = ifks
      if (emitter <= nincoming .and. emitter > 0) then
        channel_map(2, 0) = channel_map(2, 0) + 1
        emitter_position = channel_map(2, 0)
        channel_map(2, emitter_position) = ifks
      else if (emitter > nincoming .and. emitter <= nexternal) then
        channel_map(1, 0) = channel_map(1, 0) + 1
        emitter_position = channel_map(1, 0)
        channel_map(1, emitter_position) = ifks
      else
        write (*, *) 'ERROR in fks_channel_map', emitter, &
          nincoming, ifks
        stop 1
      end if
    end do
  end subroutine ensure_channel_map


  subroutine ensure_born_map()
    implicit none
    integer :: ifks, candidate, emitter, candidate_emitter

    if (allocated(born_process_map)) return
    call validate_process_dimensions()
    call validate_fks_metadata()
    allocate (born_process_map(fks_configs))
    born_process_map = 0

    do ifks = 1, fks_configs
      emitter = fks_j_d(ifks)
      if (need_color_links_d(ifks)) born_process_map(ifks) = ifks
      if (born_process_map(ifks) == 0) then
        do candidate = 1, fks_configs
          if (.not. need_color_links_d(candidate)) cycle
          candidate_emitter = fks_j_d(candidate)
          if (emitter == candidate_emitter) then
            born_process_map(ifks) = candidate
            exit
          end if
        end do
      end if
      if (born_process_map(ifks) == 0) then
        do candidate = 1, fks_configs
          if (.not. need_color_links_d(candidate)) cycle
          candidate_emitter = fks_j_d(candidate)
          if (candidate_emitter <= nincoming .and. &
              emitter <= nincoming) then
            born_process_map(ifks) = candidate
            exit
          else if (candidate_emitter > nincoming .and. &
                   emitter > nincoming) then
            born_process_map(ifks) = candidate
            exit
          end if
        end do
      end if
      if (born_process_map(ifks) == 0) then
        do candidate = 1, fks_configs
          if (need_color_links_d(candidate)) then
            born_process_map(ifks) = candidate
          end if
        end do
      end if
      if (born_process_map(ifks) == 0) born_process_map(ifks) = ifks
    end do
  end subroutine ensure_born_map


  subroutine validate_category(category)
    implicit none
    integer, intent(in) :: category

    if (category < 0 .or. category > 2) then
      call fail_channel_map('FKS initial/final category is out of range')
    end if
  end subroutine validate_category


  subroutine fail_channel_map(message)
    implicit none
    character(len=*), intent(in) :: message

    write (*, '(a)') 'ERROR in fks_channel_map: '//trim(message)
    stop 1
  end subroutine fail_channel_map

end module fks_channel_map
