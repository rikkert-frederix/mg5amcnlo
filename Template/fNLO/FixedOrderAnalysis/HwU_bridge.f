      subroutine HwU_add_points
      use HwU_module,
     &     only: module_HwU_add_points => HwU_add_points
      implicit none
      call module_HwU_add_points
      end


      subroutine HwU_accum_iter(inclde,nPSpoints,values)
      use HwU_module,
     &     only: module_HwU_accum_iter => HwU_accum_iter
      implicit none
      logical inclde
      integer nPSpoints
      double precision values(2)
      call module_HwU_accum_iter(inclde,nPSpoints,values)
      end
      subroutine HwU_output(unit,xnorm)
      use HwU_module, only: module_HwU_output => HwU_output
      implicit none
      integer unit
      double precision xnorm
      call module_HwU_output(unit,xnorm)
      end
