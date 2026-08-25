module alfas_functions_module
  implicit none
  private

  integer, parameter :: dp = kind(1d0)
  real(dp), parameter :: pi = 3.14159265358979323846d0
  real(dp), parameter :: tmass = 174d0

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

    alphas = alphas_impl(q)
  end function alphas

  double precision function alphas_impl(q)
    implicit none
    real(dp), intent(in) :: q
    integer :: iset
    real :: tafter, tbefore

    call cpu_time(tbefore)
    call getnset(iset)
    alphas_impl = alphasPDFm(iset, q)
    call cpu_time(tafter)
  end function alphas_impl


  double precision function mfrun_impl(mf, scale, asmz, nloop)
    implicit none
    real(dp), intent(in) :: mf, scale, asmz
    integer, intent(in) :: nloop
    real(dp) :: a1, as, asmf, beta0, beta1, gamma0, gamma1, l2
    integer :: nf

    if (mf > tmass) then
      nf = 6
    else
      nf = 5
    end if

    beta0 = (11d0 - 2d0 / 3d0 * nf) / 4d0
    beta1 = (102d0 - 38d0 / 3d0 * nf) / 16d0
    gamma0 = 1d0
    gamma1 = (202d0 / 3d0 - 20d0 / 9d0 * nf) / 16d0
    a1 = -beta1 * gamma0 / beta0**2 + gamma1 / beta0
    as = alphas_impl(scale)
    asmf = alphas_impl(mf)
    l2 = (1d0 + a1 * as / pi) / (1d0 + a1 * asmf / pi)

    mfrun_impl = mf * (as / asmf)**(gamma0 / beta0)
    if (nloop == 2) mfrun_impl = mfrun_impl * l2
  end function mfrun_impl

end module alfas_functions_module
