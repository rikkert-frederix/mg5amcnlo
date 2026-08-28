module multiplicative_production_channels
  implicit none
  private

  public :: multiplicative_production_channel_partition

contains

  double precision function multiplicative_production_channel_partition( &
       categories, channel_count, active_channel)
    integer, intent(in) :: categories(:)
    integer, intent(in) :: channel_count, active_channel
    integer :: channel, matching_count, active_category

    if (channel_count < 1 .or. channel_count > size(categories)) then
      call fail_multiplicative_production_channels( &
           'the production integration-channel count is invalid')
    end if
    if (active_channel < 1 .or. active_channel > channel_count) then
      call fail_multiplicative_production_channels( &
           'the active production integration channel is invalid')
    end if

    ! MINT samples one Born phase-space channel and supplies the inverse
    ! sampling probability through VEGAS_WGT.  The FKS symmetry factors
    ! already partition different initial/final radiation categories, so the
    ! remaining channel partition must sum to one separately inside the
    ! active category.  A flat partition is exact: every retained Born map
    ! covers the complete production phase space and every density tuple is
    ! evaluated with the full matrix element.  Diagram-dependent weights are
    ! only a variance optimization and would unnecessarily couple the
    ! otherwise independent production and decay density operators.
    active_category = categories(active_channel)
    matching_count = 0
    do channel = 1, channel_count
      if (categories(channel) == active_category) &
           matching_count = matching_count + 1
    end do
    if (matching_count < 1) then
      call fail_multiplicative_production_channels( &
           'the active production category has no integration channel')
    end if
    multiplicative_production_channel_partition = &
         1d0/dble(matching_count)
  end function multiplicative_production_channel_partition


  subroutine fail_multiplicative_production_channels(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_production_channels: '//trim(message)
    stop 1
  end subroutine fail_multiplicative_production_channels
end module multiplicative_production_channels
