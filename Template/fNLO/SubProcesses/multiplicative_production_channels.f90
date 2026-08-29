module multiplicative_production_channels
  implicit none
  private

  public :: multiplicative_production_channel_partition

contains

  double precision function multiplicative_production_channel_partition( &
       categories, configurations, channel_count, active_channel, &
       diagram_weights, configuration_map, symmetry_factor)
    integer, intent(in) :: categories(:), configurations(:)
    integer, intent(in) :: channel_count, active_channel
    double precision, intent(in) :: diagram_weights(:), symmetry_factor
    integer, intent(in) :: configuration_map(0:)
    double precision :: denominator, weight
    integer :: active_category, active_configuration, configuration, graph

    if (channel_count < 1 .or. channel_count > size(categories) .or. &
        channel_count > size(configurations)) then
      call fail_multiplicative_production_channels( &
           'the production integration-channel count is invalid')
    end if
    if (active_channel < 1 .or. active_channel > channel_count) then
      call fail_multiplicative_production_channels( &
           'the active production integration channel is invalid')
    end if

    active_category = categories(active_channel)
    active_configuration = configurations(active_channel)
    if (active_category < 0 .or. active_category > 2) then
      call fail_multiplicative_production_channels( &
           'the active production FKS category is invalid')
    end if
    if (configuration_map(0) < 1 .or. &
        configuration_map(0) > ubound(configuration_map, 1)) then
      call fail_multiplicative_production_channels( &
           'the Born configuration map is empty or inconsistent')
    end if
    if (active_configuration < 1 .or. &
        active_configuration > configuration_map(0)) then
      call fail_multiplicative_production_channels( &
           'the active Born configuration is invalid')
    end if
    if (symmetry_factor < 0d0 .or. symmetry_factor /= symmetry_factor) then
      call fail_multiplicative_production_channels( &
           'the diagram symmetry factor is invalid')
    end if

    ! Each production map evaluates the complete density-matrix product.
    ! Partition that common integrand with the same positive single-diagram
    ! proxy used by the original fNLO multi-channel integrator.  The weights
    ! are computed solely from the production Born amplitudes, summed over
    ! helicities, so decay blocks remain independent.  CONFIGURATION_MAP can
    ! contain the same graph more than once; SYMMETRY_FACTOR accounts for
    ! equivalent maps omitted from the MINT channel list.
    denominator = 0d0
    do configuration = 1, configuration_map(0)
      graph = configuration_map(configuration)
      if (graph < 1 .or. graph > size(diagram_weights)) then
        call fail_multiplicative_production_channels( &
             'a Born configuration references an invalid diagram')
      end if
      weight = diagram_weights(graph)
      if (weight < 0d0 .or. weight /= weight) then
        call fail_multiplicative_production_channels( &
             'a production diagram weight is invalid')
      end if
      denominator = denominator + weight
    end do

    graph = configuration_map(active_configuration)
    if (denominator == 0d0) then
      multiplicative_production_channel_partition = 0d0
    else
      multiplicative_production_channel_partition = &
           symmetry_factor*diagram_weights(graph)/denominator
    end if
  end function multiplicative_production_channel_partition


  subroutine fail_multiplicative_production_channels(message)
    character(len=*), intent(in) :: message

    write (*, '(a)') &
         'ERROR in multiplicative_production_channels: '//trim(message)
    stop 1
  end subroutine fail_multiplicative_production_channels
end module multiplicative_production_channels
