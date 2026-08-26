module analysis_dummy_module
  implicit none
  private

  public :: analysis_begin, analysis_end, analysis_fill

contains

  ! Dummy analysis routines used when no fixed-order analysis is selected.
  subroutine analysis_begin()
    implicit none
  end subroutine analysis_begin

  subroutine analysis_end()
    implicit none
  end subroutine analysis_end

  subroutine analysis_fill()
    implicit none
  end subroutine analysis_fill

end module analysis_dummy_module
