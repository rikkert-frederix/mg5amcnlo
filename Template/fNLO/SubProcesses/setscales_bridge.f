      subroutine sync_setscales_bridge()
      use kinematic_runtime_state, only: sync_kinematic_state
      implicit none
      include 'nexternal.inc'
      logical is_a_j(nexternal),is_a_lp(nexternal),is_a_lm(nexternal)
      logical is_a_ph(nexternal)
      common /to_specisa/is_a_j,is_a_lp,is_a_lm,is_a_ph
      call sync_kinematic_state(is_a_j,is_a_lp,is_a_lm,is_a_ph)
      end


      subroutine set_model_ren_scale_bridge(mur,g_value)
      implicit none
      double precision mur,g_value
      include 'coupl.inc'

      mu_r=mur
      g=g_value
      call update_as_param()
      end


      subroutine set_model_qes_scale_bridge(qes_squared)
      implicit none
      double precision qes_squared
      include 'q_es.inc'

      qes2=qes_squared
      end


      subroutine update_model_momenta_bridge(p,reset_momenta,
     &     copy_momenta)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      double precision p(0:3,nexternal)
      logical reset_momenta,copy_momenta
      double precision model_momenta(0:3,max_particles)
      common /momenta_pp/model_momenta
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
