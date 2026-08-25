      subroutine open_bash_file(lun,fname,lname)
      use ajob_file_module, only:
     $     module_open_bash_file => open_bash_file
      implicit none
      integer, intent(in) :: lun
      character(len=30), intent(inout) :: fname
      integer, intent(inout) :: lname

      call module_open_bash_file(lun,fname,lname)
      return
      end


      subroutine close_bash_file(lun)
      use ajob_file_module, only:
     $     module_close_bash_file => close_bash_file
      implicit none
      integer, intent(in) :: lun

      call module_close_bash_file(lun)
      return
      end
