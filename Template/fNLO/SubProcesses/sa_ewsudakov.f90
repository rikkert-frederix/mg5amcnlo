module standalone_ewsudakov_driver
  use process_dimensions, only: nexternal
  implicit none
  private

  public :: run_standalone_ewsudakov

contains

  subroutine run_standalone_ewsudakov()
    implicit none
    double precision, allocatable :: momenta(:,:)
    double precision :: gstrong
    double precision :: result(3)
    double precision :: time_before
    double precision :: time_after
    integer :: i

    call init_process_dimensions_bridge()
    call init_born_dimensions_bridge()
    call init_fks_metadata_bridge()

    allocate(momenta(0:3, nexternal))
    momenta = 0d0

    do
      write(*, *) 'enter gstrong'
      read(*, *) gstrong
      call update_as_param()
      write(*, *) 'enter momenta'
      do i = 1, nexternal - 1
        read(*, *) momenta(0:3, i)
      end do

      call cpu_time(time_before)
      call ewsudakov(momenta, gstrong, result)
      call cpu_time(time_after)
      write(*, *) 'RES', result
      write(*, *) 'TIME', time_after - time_before
    end do
  end subroutine run_standalone_ewsudakov

end module standalone_ewsudakov_driver
