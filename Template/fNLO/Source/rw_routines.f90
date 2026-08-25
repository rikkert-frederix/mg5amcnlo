module rw_routines
  implicit none
  private

  public :: case_trap2

contains

  subroutine case_trap2(name)
    ! Convert the ASCII uppercase characters in NAME to lowercase.
    character(len=*), intent(inout) :: name

    integer :: i
    integer :: k

    do i = 1, len(name)
      k = ichar(name(i:i))
      if (k >= 65 .and. k <= 90) then
        k = k + 32
        name(i:i) = char(k)
      end if
    end do
  end subroutine case_trap2


  subroutine to_upper(name)
    ! Preserve the legacy conversion semantics exactly: every character
    ! whose ASCII code is greater than 90 is shifted by 32 positions.
    character(len=*), intent(inout) :: name

    integer :: i
    integer :: k

    do i = 1, len(name)
      k = ichar(name(i:i))
      if (k > 90) then
        k = k - 32
        name(i:i) = char(k)
      end if
    end do
  end subroutine to_upper

end module rw_routines
