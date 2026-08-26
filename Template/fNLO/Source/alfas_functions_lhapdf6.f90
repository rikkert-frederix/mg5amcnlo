module alfas_functions_module
  implicit none
  private

  integer, parameter :: dp = kind(1d0)
  interface
    subroutine getnset(iset)
      integer, intent(out) :: iset
    end subroutine getnset

    double precision function alphasPDFm(iset, q)
      integer, intent(in) :: iset
      double precision, intent(in) :: q
    end function alphasPDFm
  end interface

  public :: alphas

contains

  double precision function alphas(q)
    implicit none
    real(dp), intent(in) :: q
    integer :: iset

    call getnset(iset)
    alphas = alphasPDFm(iset, q)
  end function alphas

end module alfas_functions_module
