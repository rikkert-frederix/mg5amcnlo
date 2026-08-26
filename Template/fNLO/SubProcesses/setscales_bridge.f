      subroutine set_model_ren_scale_bridge(mur,g_value)
      implicit none
      double precision mur,g_value
      include 'coupl.inc'

      mu_r=mur
      g=g_value
      call update_as_param()
      end


      subroutine set_model_qes_scale_bridge(qes_squared)
      use fnlo_process_common, only: qes2
      implicit none
      double precision qes_squared

      qes2=qes_squared
      end


      subroutine update_model_momenta_bridge(p,reset_momenta,
     &     copy_momenta)
      use fnlo_process_common, only: nexternal,max_particles,
     $     model_momenta
      implicit none
      double precision p(0:3,nexternal)
      logical reset_momenta,copy_momenta
      integer i,j

      if (reset_momenta) then
         do j=1,max_particles
            do i=0,3
               model_momenta(i,j)=0d0
            enddo
         enddo
      endif
      if (copy_momenta) then
         do j=1,nexternal
            do i=0,3
               model_momenta(i,j)=p(i,j)
            enddo
         enddo
      endif
      end
