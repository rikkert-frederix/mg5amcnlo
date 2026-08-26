module open_output_files_module
  implicit none
  private

  public :: HwU_write_file

  interface
    subroutine HwU_output(output_unit, xnorm)
      implicit none
      integer, intent(in) :: output_unit
      double precision, intent(in) :: xnorm
    end subroutine HwU_output
  end interface

contains

  subroutine HwU_write_file()
    implicit none
    double precision :: xnorm

    open(unit=99, file='MADatNLO.HwU', status='unknown')
    xnorm = 1d0
    call HwU_output(99, xnorm)
    close(99)
  end subroutine HwU_write_file

end module open_output_files_module
