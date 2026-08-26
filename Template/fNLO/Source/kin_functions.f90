module kin_functions_module
  implicit none
  private

  public :: delta_phi_impl
  public :: dot_impl
  public :: et_impl
  public :: eta_impl
  public :: pt_impl
  public :: rho_impl
  public :: switchmom_impl
  public :: theta_impl

contains

  double precision function delta_phi_impl(p1, p2)
    implicit none
    double precision, intent(in) :: p1(0:3), p2(0:3)
    double precision :: denom, temp

    denom = sqrt(p1(1)**2 + p1(2)**2) * &
      sqrt(p2(1)**2 + p2(2)**2)
    temp = max(-0.99999999d0, &
      (p1(1) * p2(1) + p1(2) * p2(2)) / denom)
    temp = min(0.99999999d0, temp)
    delta_phi_impl = acos(temp)
  end function delta_phi_impl


  double precision function et_impl(p)
    implicit none
    double precision, intent(in) :: p(0:3)
    double precision :: transverse_momentum

    transverse_momentum = sqrt(p(1)**2 + p(2)**2)
    if (transverse_momentum > 0d0) then
      et_impl = p(0) * transverse_momentum / &
        sqrt(transverse_momentum**2 + p(3)**2)
    else
      et_impl = 0d0
    end if
  end function et_impl


  double precision function pt_impl(p)
    implicit none
    double precision, intent(in) :: p(0:3)

    pt_impl = sqrt(p(1)**2 + p(2)**2)
  end function pt_impl


  subroutine switchmom_impl(p1, p, ic, jc, nexternal)
    implicit none
    integer, intent(in) :: nexternal
    double precision, intent(in) :: p1(0:3, nexternal)
    double precision, intent(out) :: p(0:3, nexternal)
    integer, intent(in) :: ic(nexternal)
    integer, intent(out) :: jc(nexternal)
    integer :: i, j

    do i = 1, nexternal
      do j = 0, 3
        p(j, ic(i)) = p1(j, i)
      end do
    end do
    do i = 1, nexternal
      jc(i) = 1
    end do
    jc(ic(1)) = -1
    jc(ic(2)) = -1
  end subroutine switchmom_impl


  double precision function dot_impl(p1, p2)
    implicit none
    double precision, intent(in) :: p1(0:3), p2(0:3)

    dot_impl = p1(0) * p2(0) - p1(1) * p2(1) - &
      p1(2) * p2(2) - p1(3) * p2(3)
    if (abs(dot_impl) < 1d-6) dot_impl = 0d0
  end function dot_impl


  double precision function threedot_impl(p1, p2)
    implicit none
    double precision, intent(in) :: p1(0:3), p2(0:3)

    threedot_impl = p1(1) * p2(1) + p1(2) * p2(2) + p1(3) * p2(3)
  end function threedot_impl


  double precision function rho_impl(p1)
    implicit none
    double precision, intent(in) :: p1(0:3)

    rho_impl = sqrt(threedot_impl(p1, p1))
  end function rho_impl


  double precision function theta_impl(p)
    implicit none
    double precision, intent(in) :: p(0:3)

    theta_impl = acos(p(3) / sqrt(p(1)**2 + p(2)**2 + p(3)**2))
  end function theta_impl


  double precision function eta_impl(p)
    implicit none
    double precision, intent(in) :: p(0:3)
    double precision, parameter :: pi = &
      3.14159265358979323846264338327950d0
    double precision :: polar_angle

    polar_angle = theta_impl(p)
    if (abs(polar_angle) < 1d-5) then
      eta_impl = 25d0
    else if (abs(polar_angle - pi) < 1d-5) then
      eta_impl = -25d0
    else if (polar_angle /= polar_angle) then
      eta_impl = -99d99
    else
      eta_impl = -log(tan(polar_angle / 2d0))
    end if
  end function eta_impl


end module kin_functions_module
