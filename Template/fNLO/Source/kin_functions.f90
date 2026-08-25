module kin_functions_module
  implicit none
  private

  public :: delta_phi_impl
  public :: dj
  public :: dj1_impl
  public :: dj1
  public :: dj2_impl
  public :: dj_impl
  public :: djb
  public :: djb_impl
  public :: dot_impl
  public :: et_impl
  public :: eta_impl
  public :: four_momentum_impl
  public :: four_momentum_set2_impl
  public :: phi_impl
  public :: pt_impl
  public :: ptdot_impl
  public :: r2_impl
  public :: rap
  public :: rap_impl
  public :: rap2_impl
  public :: rho_impl
  public :: sumdot_impl
  public :: switchhel_impl
  public :: switchmom_impl
  public :: theta_impl
  public :: threedot_impl

contains

  double precision function rap(p)
    use run_state, only: ebeam, xbk
    implicit none
    double precision, intent(in) :: p(0:3)

    rap = rap_impl(p, xbk, ebeam)
  end function rap


  double precision function dj(p1, p2)
    use run_state, only: jet_distance_parameter, lpp
    implicit none
    double precision, intent(in) :: p1(0:4), p2(0:4)

    dj = dj_impl(p1, p2, lpp, jet_distance_parameter)
  end function dj


  double precision function dj1(p1, p2)
    use run_state, only: lpp
    implicit none
    double precision, intent(in) :: p1(0:3), p2(0:3)

    dj1 = dj1_impl(p1, p2, lpp)
  end function dj1


  double precision function djb(p1)
    use run_state, only: lpp
    implicit none
    double precision, intent(in) :: p1(0:4)

    djb = djb_impl(p1, lpp)
  end function djb

  double precision function r2_impl(p1, p2)
    implicit none
    double precision, intent(in) :: p1(0:3), p2(0:3)

    r2_impl = delta_phi_impl(p1, p2)**2 + &
      (eta_impl(p1) - eta_impl(p2))**2
  end function r2_impl


  double precision function sumdot_impl(p1, p2, dsign)
    implicit none
    double precision, intent(in) :: p1(0:3), p2(0:3), dsign
    double precision :: ptot(0:3)

    ptot = p1 + dsign * p2
    sumdot_impl = dot_impl(ptot, ptot)
  end function sumdot_impl


  double precision function ptdot_impl(p1, p2)
    implicit none
    double precision, intent(in) :: p1(0:3), p2(0:3)

    ptdot_impl = (p1(1) + p2(1))**2 + (p1(2) + p2(2))**2
  end function ptdot_impl


  double precision function rap_impl(p, xbk, ebeam)
    implicit none
    double precision, intent(in) :: p(0:3), xbk(2), ebeam(2)
    double precision :: pm

    pm = p(0)
    rap_impl = 0.5d0 * log((pm + p(3)) / (pm - p(3))) + &
      0.5d0 * log(xbk(1) * ebeam(1) / (xbk(2) * ebeam(2)))
  end function rap_impl


  double precision function rap2_impl(p)
    implicit none
    double precision, intent(in) :: p(0:3)
    double precision :: pm

    pm = p(0)
    rap2_impl = 0.5d0 * log((pm + p(3)) / (pm - p(3)))
  end function rap2_impl


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


  double precision function dj_impl(p1, p2, lpp, d)
    implicit none
    double precision, intent(in) :: p1(0:4), p2(0:4), d
    integer, intent(in) :: lpp(2)
    double precision :: costh, eta1, eta2, p1a, p2a, pt1, pt2

    if ((lpp(1) /= 1) .and. (lpp(2) /= 1)) then
      p1a = sqrt(p1(1)**2 + p1(2)**2 + p1(3)**2)
      p2a = sqrt(p2(1)**2 + p2(2)**2 + p2(3)**2)
      if (p1a * p2a /= 0d0) then
        costh = (p1(1) * p2(1) + p1(2) * p2(2) + &
          p1(3) * p2(3)) / (p1a * p2a)
        dj_impl = 2d0 * min(p1(0)**2, p2(0)**2) * (1d0 - costh)
      else
        dj_impl = 0d0
      end if
    else
      pt1 = p1(1)**2 + p1(2)**2
      pt2 = p2(1)**2 + p2(2)**2
      if (pt1 == 0d0 .or. pt2 == 0d0) then
        dj_impl = 0d0
        return
      end if
      p1a = sqrt(pt1 + p1(3)**2)
      p2a = sqrt(pt2 + p2(3)**2)
      eta1 = 0.5d0 * log((p1a + p1(3)) / (p1a - p1(3)))
      eta2 = 0.5d0 * log((p2a + p2(3)) / (p2a - p2(3)))
      dj_impl = max(p1(4), p2(4)) + min(pt1, pt2) * 2d0 * &
        (cosh(eta1 - eta2) - &
        (p1(1) * p2(1) + p1(2) * p2(2)) / sqrt(pt1 * pt2)) / d**2
    end if
  end function dj_impl


  double precision function dj1_impl(p1, p2, lpp)
    implicit none
    double precision, intent(in) :: p1(0:3), p2(0:3)
    integer, intent(in) :: lpp(2)
    double precision :: costh, eta1, eta2, p1a, p2a, pt1, pt2, ptm1

    if ((lpp(1) == 0) .and. (lpp(2) == 0)) then
      p1a = sqrt(p1(1)**2 + p1(2)**2 + p1(3)**2)
      p2a = sqrt(p2(1)**2 + p2(2)**2 + p2(3)**2)
      if (p1a * p2a /= 0d0) then
        costh = (p1(1) * p2(1) + p1(2) * p2(2) + &
          p1(3) * p2(3)) / (p1a * p2a)
        dj1_impl = 2d0 * p1(0)**2 * (1d0 - costh)
      else
        dj1_impl = 0d0
      end if
    else
      pt1 = p1(1)**2 + p1(2)**2
      pt2 = p2(1)**2 + p2(2)**2
      p1a = sqrt(pt1 + p1(3)**2)
      p2a = sqrt(pt2 + p2(3)**2)
      eta1 = 0.5d0 * log((p1a + p1(3)) / (p1a - p1(3)))
      eta2 = 0.5d0 * log((p2a + p2(3)) / (p2a - p2(3)))
      ptm1 = max((p1(0) - p1(3)) * (p1(0) + p1(3)), 0d0)
      dj1_impl = 2d0 * ptm1 * (cosh(eta1 - eta2) - &
        (p1(1) * p2(1) + p1(2) * p2(2)) / sqrt(pt1 * pt2))
    end if
  end function dj1_impl


  double precision function djb_impl(p1, lpp)
    implicit none
    double precision, intent(in) :: p1(0:4)
    integer, intent(in) :: lpp(2)

    if ((lpp(1) == 0) .and. (lpp(2) == 0)) then
      djb_impl = max(p1(0), 0d0)**2
    else
      djb_impl = (p1(0) - p1(3)) * (p1(0) + p1(3))
    end if
  end function djb_impl


  double precision function dj2_impl(p1, p2)
    implicit none
    double precision, intent(in) :: p1(0:3), p2(0:3)

    dj2_impl = dot_impl(p1, p1) + 2d0 * dot_impl(p1, p2) + &
      dot_impl(p2, p2)
  end function dj2_impl


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


  subroutine switchhel_impl(hel, hel1, ic, nexternal)
    implicit none
    integer, intent(in) :: nexternal
    integer, intent(in) :: hel(nexternal), ic(nexternal)
    integer, intent(out) :: hel1(nexternal)
    integer :: i

    do i = 1, nexternal
      hel1(ic(i)) = hel(i)
    end do
  end subroutine switchhel_impl


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


  subroutine four_momentum_impl(theta, phi, rho, mass, p)
    implicit none
    double precision, intent(in) :: theta, phi, rho, mass
    double precision, intent(out) :: p(0:3)

    p(1) = rho * sin(theta) * cos(phi)
    p(2) = rho * sin(theta) * sin(phi)
    p(3) = rho * cos(theta)
    p(0) = sqrt(rho**2 + mass**2)
  end subroutine four_momentum_impl


  subroutine four_momentum_set2_impl(eta, phi, pt, mass, p)
    implicit none
    double precision, intent(in) :: eta, phi, pt, mass
    double precision, intent(out) :: p(0:3)

    p(1) = pt * cos(phi)
    p(2) = pt * sin(phi)
    p(3) = pt * sinh(eta)
    p(0) = sqrt(p(1)**2 + p(2)**2 + p(3)**2 + mass**2)
  end subroutine four_momentum_set2_impl


  double precision function phi_impl(p)
    implicit none
    double precision, intent(in) :: p(0:3)
    double precision, parameter :: pi = 3.141592654d0

    if (p(1) > 0d0) then
      phi_impl = atan(p(2) / p(1))
    else if (p(1) < 0d0) then
      phi_impl = atan(p(2) / p(1)) + pi
    else if (p(2) >= 0d0) then
      phi_impl = pi / 2d0
    else
      phi_impl = -pi / 2d0
    end if
    if (phi_impl < 0d0) phi_impl = phi_impl + 2d0 * pi
  end function phi_impl

end module kin_functions_module
