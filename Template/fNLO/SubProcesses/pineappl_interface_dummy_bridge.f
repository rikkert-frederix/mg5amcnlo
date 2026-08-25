c Legacy PineAPPL ABI for the module-owned dummy implementation.

      subroutine APPL_fill()
      use pineappl_dummy_module, only: module_APPL_fill => APPL_fill
      implicit none

      call module_APPL_fill()
      return
      end


      subroutine APPL_init()
      use pineappl_dummy_module, only: module_APPL_init => APPL_init
      implicit none

      call module_APPL_init()
      return
      end


      subroutine APPL_delete_itype()
      use pineappl_dummy_module, only:
     &     module_APPL_delete_itype => APPL_delete_itype
      implicit none

      call module_APPL_delete_itype()
      return
      end


      subroutine APPL_term()
      use pineappl_dummy_module, only: module_APPL_term => APPL_term
      implicit none

      call module_APPL_term()
      return
      end
