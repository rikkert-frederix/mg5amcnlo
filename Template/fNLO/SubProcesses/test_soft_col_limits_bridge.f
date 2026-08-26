      program test_soft_col_limits
c-----------------------------------------------------------------------
c     Keep the historical program target while the retained matrix-
c     element limit tests live in an include-free Fortran 90 module.
c-----------------------------------------------------------------------
      use test_soft_col_limits_module, only:
     $     run_test_soft_col_limits,finalize_test_soft_col_limits
      implicit none
      interface
         subroutine init_process_dimensions_bridge()
         end subroutine init_process_dimensions_bridge
         subroutine init_born_dimensions_bridge()
         end subroutine init_born_dimensions_bridge
         subroutine init_fks_metadata_bridge()
         end subroutine init_fks_metadata_bridge
         subroutine init_genps_fks_bridge()
         end subroutine init_genps_fks_bridge
      end interface

      call init_process_dimensions_bridge()
      call init_born_dimensions_bridge()
      call init_fks_metadata_bridge()
      call init_genps_fks_bridge()
      call run_test_soft_col_limits()
      call finalize_test_soft_col_limits()
      end


      subroutine init_test_limits_data_bridge()
c-----------------------------------------------------------------------
c     Materialize Python-generated Born configurations and model masses,
c     then copy them into allocatable module state.
c-----------------------------------------------------------------------
      use test_soft_col_limits_module, only:
     $     init_test_limits_data
      use fnlo_process_common, only: nexternal
      implicit none
      include 'coupl.inc'
      double precision zero
      parameter (zero=0d0)
      double precision pmass(nexternal)
      include 'born_conf.inc'
      include 'pmass.inc'

      call init_test_limits_data(mapconfig,pmass)
      return
      end


      subroutine test_limits_select_fks_bridge(configuration,
     $     i_fks_out,j_fks_out)
c-----------------------------------------------------------------------
c     Select the generated FKS and Les Houches records and expose only
c     the selected indices to the module implementation.
c-----------------------------------------------------------------------
      use fnlo_process_common, only: nfksprocess,i_fks,j_fks
      implicit none
      integer configuration,i_fks_out,j_fks_out
      interface
         subroutine fks_inc_chooser()
         end subroutine fks_inc_chooser
         subroutine leshouche_inc_chooser()
         end subroutine leshouche_inc_chooser
      end interface

      nfksprocess=configuration
      call fks_inc_chooser()
      call leshouche_inc_chooser()
      i_fks_out=i_fks
      j_fks_out=j_fks
      return
      end


      subroutine test_limits_set_nndim_bridge(nndim_value)
      use fnlo_process_common, only: nndim
      implicit none
      integer nndim_value

      nndim=nndim_value
      return
      end


      subroutine test_limits_set_controls_bridge(xi_fixed,y_fixed,
     $     calculated_born_in,soft_test_in,collinear_test_in)
      use fnlo_process_common, only: xi_i_fks_fix,y_ij_fks_fix,
     $     calculated_born,softtest,colltest
      implicit none
      double precision xi_fixed,y_fixed
      logical calculated_born_in,soft_test_in,collinear_test_in

      xi_i_fks_fix=xi_fixed
      y_ij_fks_fix=y_fixed
      calculated_born=calculated_born_in
      softtest=soft_test_in
      colltest=collinear_test_in
      return
      end


      subroutine test_limits_sync_state_bridge(p1_out,
     $     jac_out,p_born_out,xi_event_out,y_event_out,p_i_event_out,
     $     p_i_counter_out,xi_counter_out)
c-----------------------------------------------------------------------
c     Snapshot the state produced by GENPS_FKS.  Weight-only counterevent
c     arrays are deliberately omitted because this test never reads them.
c-----------------------------------------------------------------------
      use fnlo_process_common, only: nexternal,p1_cnt,jac_cnt,
     $     p_born,xi_i_fks_ev,y_ij_fks_ev,p_i_fks_ev,p_i_fks_cnt,
     $     xi_i_fks_cnt
      implicit none
      double precision p1_out(0:3,nexternal,0:2)
      double precision jac_out(0:2)
      double precision p_born_out(0:3,nexternal-1)
      double precision xi_event_out,y_event_out
      double precision p_i_event_out(0:3)
      double precision p_i_counter_out(0:3,0:2)
      double precision xi_counter_out(0:2)
      p1_out=p1_cnt
      jac_out=jac_cnt
      p_born_out=p_born
      xi_event_out=xi_i_fks_ev
      y_event_out=y_ij_fks_ev
      p_i_event_out=p_i_fks_ev
      p_i_counter_out=p_i_fks_cnt
      xi_counter_out=xi_i_fks_cnt
      return
      end


      subroutine test_limits_setcuts_bridge(etmin_out,etmax_out,
     $     mxxmin_out)
      use fnlo_process_common, only: nexternal,nincoming,etmin,etmax,
     $     mxxmin
      implicit none
      double precision etmin_out(nincoming+1:nexternal-1)
      double precision etmax_out(nincoming+1:nexternal-1)
      double precision mxxmin_out(nincoming+1:nexternal-1,
     $     nincoming+1:nexternal-1)
      interface
         subroutine setcuts()
         end subroutine setcuts
      end interface

      call setcuts()
      etmin_out=etmin
      etmax_out=etmax
      mxxmin_out=mxxmin
      return
      end


      subroutine test_limits_sreal_bridge(momentum,xi_fks,y_fks,
     $     weight,split_weights)
      use fks_singular_module, only: sreal
      use fnlo_process_common, only: nexternal,amp_split_size,
     $     amp_split
      implicit none
      double precision momentum(0:3,nexternal)
      double precision xi_fks,y_fks,weight
      double precision split_weights(amp_split_size)
      integer iamp
      call sreal(momentum,xi_fks,y_fks,weight)
      do iamp=1,amp_split_size
         split_weights(iamp)=amp_split(iamp)
      enddo
      return
      end
