module hwu_dummy_module
  implicit none
  private

  public :: HwU_output, HwU_add_points, HwU_accum_iter, accum, addfil

contains

  subroutine HwU_output(idummy, dummy)
    implicit none
    integer, intent(in) :: idummy
    double precision, intent(in) :: dummy

    write (*, *) 'HwU_output should not be called', idummy
    stop 1
  end subroutine HwU_output

  subroutine HwU_add_points()
    implicit none
  end subroutine HwU_add_points

  subroutine HwU_accum_iter(ldummy, idummy, dummy)
    implicit none
    logical, intent(in) :: ldummy
    integer, intent(in) :: idummy
    double precision, intent(in) :: dummy(2)
  end subroutine HwU_accum_iter

  subroutine accum(idummy)
    implicit none
    integer, intent(in) :: idummy
  end subroutine accum

  subroutine addfil(string)
    implicit none
    character(len=*), intent(inout) :: string
  end subroutine addfil

end module hwu_dummy_module
