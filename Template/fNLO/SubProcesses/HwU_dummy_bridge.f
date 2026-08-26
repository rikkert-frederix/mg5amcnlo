c Legacy HwU ABI for the module-owned dummy implementation.

      subroutine HwU_output(idummy,dummy)
      use hwu_dummy_module, only: module_HwU_output => HwU_output
      implicit none
      integer idummy
      double precision dummy

      call module_HwU_output(idummy)
      return
      end


      subroutine HwU_add_points()
      use hwu_dummy_module, only:
     &     module_HwU_add_points => HwU_add_points
      implicit none

      call module_HwU_add_points()
      return
      end


      subroutine HwU_accum_iter(ldummy,idummy,dummy)
      use hwu_dummy_module, only:
     &     module_HwU_accum_iter => HwU_accum_iter
      implicit none
      logical ldummy
      integer idummy
      double precision dummy(2)

      call module_HwU_accum_iter()
      return
      end
