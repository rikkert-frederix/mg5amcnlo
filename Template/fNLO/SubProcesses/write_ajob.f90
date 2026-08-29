module ajob_file_module
  implicit none
  private

  integer :: file_counter = 0

  public :: open_bash_file
  public :: close_bash_file

contains



  subroutine open_bash_file(lun, fname, lname)
    implicit none
    integer, intent(in) :: lun
    character(len=30), intent(inout) :: fname
    integer, intent(inout) :: lname
    character(len=256) :: buffer

    file_counter = file_counter + 1
    if (file_counter < 10) then
      write (fname(5:5), '(i1)') file_counter
      lname = lname + 1
    else if (file_counter < 100) then
      write (fname(5:6), '(i2)') file_counter
      lname = lname + 2
    else if (file_counter < 1000) then
      write (fname(5:7), '(i3)') file_counter
      lname = lname + 3
    end if

    open(unit=lun, file=fname, status='unknown')
    open(unit=lun + 1, file='../ajob_template', status='old')
    do
      read(lun + 1, '(a)', err=100, end=100) buffer
      if (index(buffer, 'TAGTAGTAGTAGTAG') /= 0) exit
      write(lun, '(a)') buffer
    end do
    write(lun, '(a)', advance='no') 'for i in $channel '
    return

100 continue
    write(*, *) 'ajob_template or ajob_template_cluster ', &
                'does not have the correct format'
    stop
  end subroutine open_bash_file


  subroutine close_bash_file(lun)
    implicit none
    integer, intent(in) :: lun
    character(len=256) :: buffer

    write(lun, '(a)') '; do'
    do
      read(lun + 1, '(a)', err=100, end=100) buffer
      write(lun, '(a)') buffer
    end do

100 continue
    close(lun + 1)
    close(lun)
  end subroutine close_bash_file

end module ajob_file_module
