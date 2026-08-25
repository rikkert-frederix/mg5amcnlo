module alfas_functions_module
  implicit none
  private

  integer, parameter :: dp = kind(1d0)
  real(dp), parameter :: pi = 3.14159265358979323846d0
  real(dp), parameter :: zmass = 91.188d0
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

  public :: alfa_impl
  public :: alfaw_impl
  public :: alphas
  public :: alphas_impl
  public :: mfrun_impl

contains

  double precision function alphas(q)
    implicit none
    real(dp), intent(in) :: q

    alphas = alphas_impl(q)
  end function alphas

  double precision function alfa_impl(alfa0, qsq)
    implicit none
    real(dp), intent(in) :: alfa0, qsq

    alfa_impl = alfa0 / &
         (1d0 - alfa0 * log(qsq / zmass**2) / 3d0 / pi)
  end function alfa_impl


  double precision function alfaw_impl(alfaw0, qsq, nh)
    implicit none
    real(dp), intent(in) :: alfaw0, qsq
    integer, intent(in) :: nh
    real(dp) :: dum
    integer :: nq

    if (qsq >= tmass**2) then
      nq = 6
    else
      nq = 5
    end if
    dum = (22d0 - 4d0 * nq - nh / 2d0) / (12d0 * pi)
    alfaw_impl = alfaw0 / &
         (1d0 + dum * alfaw0 * log(qsq / zmass**2))
  end function alfaw_impl


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
