module open_output_files_module
  implicit none
  private

  public :: HwU_write_file, HwU_write_file_impl

  interface
    ! Explicit interface to the legacy external bridge used by callers.
    subroutine HwU_write_file()
      implicit none
    end subroutine HwU_write_file

    subroutine APPL_term()
      implicit none
    end subroutine APPL_term

    subroutine HwU_output(output_unit, xnorm)
      implicit none
      integer, intent(in) :: output_unit
      double precision, intent(in) :: xnorm
    end subroutine HwU_output
  end interface

contains

  subroutine HwU_write_file_impl(pineappl_enabled, num_observables, &
       observable_number)
    implicit none
    logical, intent(in) :: pineappl_enabled
    integer, intent(in) :: num_observables
    integer, intent(inout) :: observable_number
    double precision :: xnorm
    integer :: j

    if (pineappl_enabled) then
      do j = 1, num_observables
        observable_number = j
        call APPL_term()
      end do
    end if

    open(unit=99, file='MADatNLO.HwU', status='unknown')
    xnorm = 1d0
    call HwU_output(99, xnorm)
    close(99)
  end subroutine HwU_write_file_impl

end module open_output_files_module
