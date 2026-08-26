module fks_random_module
  use ranmar_module, only: ntuple
  implicit none
  private

  public :: random_unit_interval

contains

  double precision function random_unit_interval(configuration)
    implicit none
    integer, intent(in) :: configuration
    double precision :: lower_bound, upper_bound

    lower_bound = 0d0
    upper_bound = 1d0
    call ntuple(random_unit_interval, lower_bound, upper_bound, &
                configuration)
  end function random_unit_interval

end module fks_random_module
