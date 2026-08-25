module gaussian_integration
  implicit none
  private
  public :: dgauss

contains

  double precision function dgauss(f, a, b, eps) result(h)
    implicit none

    interface
      double precision function f(x)
        implicit none
        double precision, intent(in) :: x
      end function f
    end interface

    double precision, intent(in) :: a, b, eps
    character(len=6), parameter :: name = 'DGAUSS'
    double precision, parameter :: z1 = 1d0
    double precision, parameter :: hf = z1 / 2d0
    double precision, parameter :: cst = 5d0 * z1 / 1000d0
    double precision, parameter :: x(12) = (/ &
         9.6028985649753623d-1, 7.9666647741362674d-1, &
         5.2553240991632899d-1, 1.8343464249564980d-1, &
         9.8940093499164993d-1, 9.4457502307323258d-1, &
         8.6563120238783174d-1, 7.5540440835500303d-1, &
         6.1787624440264375d-1, 4.5801677765722739d-1, &
         2.8160355077925891d-1, 9.5012509837637440d-2 /)
    double precision, parameter :: w(12) = (/ &
         1.0122853629037626d-1, 2.2238103445337447d-1, &
         3.1370664587788729d-1, 3.6268378337836198d-1, &
         2.7152459411754095d-2, 6.2253523938647893d-2, &
         9.5158511682492785d-2, 1.2462897125553387d-1, &
         1.4959598881657673d-1, 1.6915651939500254d-1, &
         1.8260341504492359d-1, 1.8945061045506850d-1 /)
    double precision :: aa, bb, c1, c2, const, refinement_test
    double precision :: s8, s16, u
    integer :: i

    h = 0d0
    if (b <= a .and. b >= a) return
    const = cst / abs(b - a)
    bb = a

10  aa = bb
    bb = b

20  c1 = hf * (bb + aa)
    c2 = hf * (bb - aa)
    s8 = 0d0
    do i = 1, 4
      u = c2 * x(i)
      s8 = s8 + w(i) * (f(c1 + u) + f(c1 - u))
    end do

    s16 = 0d0
    do i = 5, 12
      u = c2 * x(i)
      s16 = s16 + w(i) * (f(c1 + u) + f(c1 - u))
    end do
    s16 = c2 * s16

    if (abs(s16 - c2 * s8) <= eps * (1d0 + abs(s16))) then
      h = h + s16
      if (.not. (bb <= b .and. bb >= b)) go to 10
    else
      bb = c1
      refinement_test = 1d0 + const * abs(c2)
      if (.not. (refinement_test <= 1d0 .and. &
                 refinement_test >= 1d0)) go to 20
      h = 0d0
      write (*, *) name, 'ERROR: TOO HIGH ACCURACY REQUIRED'
    end if
  end function dgauss

end module gaussian_integration
