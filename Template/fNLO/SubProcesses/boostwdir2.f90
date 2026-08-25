module boostwdir2_module
  implicit none
  private

  public :: boostwdir2
  public :: boostwdir2_in_place

contains

  subroutine boostwdir2(chybst, shybst, chybstmo, xd, xin, xout)
    ! chybstmo = chybst - 1; computing it analytically improves the
    ! numerical accuracy.
    implicit none

    double precision, intent(in) :: chybst, shybst, chybstmo
    double precision, intent(in) :: xd(1:3), xin(0:3)
    double precision, intent(out) :: xout(0:3)
    double precision :: en, pz
    integer :: i

    if (abs(xd(1)**2 + xd(2)**2 + xd(3)**2 - 1.0d0) > 1.0d-6) then
      write (*, *) 'Error #1 in boostwdir2', xd
      stop
    end if

    en = xin(0)
    pz = xin(1)*xd(1) + xin(2)*xd(2) + xin(3)*xd(3)
    xout(0) = en*chybst - pz*shybst
    do i = 1, 3
      xout(i) = xin(i) + xd(i)*(pz*chybstmo - en*shybst)
    end do
  end subroutine boostwdir2


  subroutine boostwdir2_in_place(chybst, shybst, chybstmo, xd, momentum)
    implicit none

    double precision, intent(in) :: chybst, shybst, chybstmo
    double precision, intent(in) :: xd(1:3)
    double precision, intent(inout) :: momentum(0:3)
    double precision :: input_momentum(0:3)

    input_momentum = momentum
    call boostwdir2(chybst, shybst, chybstmo, xd, input_momentum, &
         momentum)
  end subroutine boostwdir2_in_place

end module boostwdir2_module
