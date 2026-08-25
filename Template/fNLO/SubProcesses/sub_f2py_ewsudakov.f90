module ewsudakov_python_interface
  implicit none
  private

  public :: evaluate_ewsudakov

contains

  subroutine evaluate_ewsudakov(p_born, nexternal, gstrong, results)
    implicit none
    integer, intent(in) :: nexternal
    double precision, intent(in) :: p_born(0:3, nexternal)
    double precision, intent(in) :: gstrong
    double precision, intent(out) :: results(6)

    call ewsudakov_f77(p_born, gstrong, results)
  end subroutine evaluate_ewsudakov

end module ewsudakov_python_interface
