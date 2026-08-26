module hwu_dummy_module
  implicit none
  private

  public :: HwU_output, HwU_add_points, HwU_accum_iter

contains

  subroutine HwU_output(idummy)
    implicit none
    integer, intent(in) :: idummy

    write (*, *) 'HwU_output should not be called', idummy
    stop 1
  end subroutine HwU_output

  subroutine HwU_add_points()
    implicit none
  end subroutine HwU_add_points

  subroutine HwU_accum_iter()
    implicit none
  end subroutine HwU_accum_iter

end module hwu_dummy_module
