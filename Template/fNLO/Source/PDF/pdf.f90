module pdf_dispatch_module
  implicit none
  private

  integer, parameter :: dp = kind(1.0d0)
  real(dp) :: nnfx_workspace(-6:7)

  public :: pftopdg_impl

  interface
    subroutine nnevolvepdf(x, q, xpdf)
      double precision, intent(in) :: x
      double precision, intent(in) :: q
      double precision, intent(out) :: xpdf(-6:7)
    end subroutine nnevolvepdf
  end interface

contains

  subroutine pftopdg_impl(ih, x, q, pdf, pdlabel)
    implicit none
    integer, intent(in) :: ih
    real(dp), intent(in) :: x, q
    real(dp), intent(out) :: pdf(-7:7)
    character(len=*), intent(in) :: pdlabel

    call fdist_impl(ih, x, q, pdf, pdlabel)
  end subroutine pftopdg_impl


  subroutine fdist_impl(ih, x, xmu, fx, pdlabel)
    implicit none
    integer, intent(in) :: ih
    real(dp), intent(in) :: x, xmu
    real(dp), intent(out) :: fx(-7:7)
    character(len=*), intent(in) :: pdlabel

    fx = 0d0
    if (x >= 1d0) return

    if (abs(ih) /= 1) then
      write (*, *) 'The bundled NNPDF backend supports proton beams only'
      stop 1
    end if

    select case (trim(pdlabel))
    case ('nn23lo', 'nn23lo1', 'nn23nlo')
      call nnevolvepdf(x, xmu, nnfx_workspace)
      fx(-5:5) = nnfx_workspace(-5:5) / x
      fx(7) = nnfx_workspace(7) / x
    case default
      write (*, *) 'Unsupported bundled PDF label: ', trim(pdlabel)
      write (*, *) 'Use nn23lo, nn23lo1, nn23nlo, or LHAPDF 6'
      stop 1
    end select
  end subroutine fdist_impl

end module pdf_dispatch_module
