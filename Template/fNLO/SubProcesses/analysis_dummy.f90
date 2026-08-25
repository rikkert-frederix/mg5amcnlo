module analysis_dummy_module
  use process_dimensions, only: nexternal
  implicit none
  private

  public :: analysis_begin, analysis_end, analysis_fill

contains

  ! Dummy analysis routines used when no fixed-order analysis is selected.
  subroutine analysis_begin(nwgt, weights_info)
    implicit none
    integer, intent(in) :: nwgt
    character(len=*), intent(in) :: weights_info(*)
  end subroutine analysis_begin

  subroutine analysis_end(xnorm)
    implicit none
    double precision, intent(in) :: xnorm
  end subroutine analysis_end

  subroutine analysis_fill(p, istatus, ipdg, wgts, ibody)
    implicit none
    double precision, intent(in) :: p(0:4, nexternal)
    integer, intent(in) :: istatus(nexternal), ipdg(nexternal)
    double precision, intent(in) :: wgts(*)
    integer, intent(in) :: ibody
  end subroutine analysis_fill

end module analysis_dummy_module
