      program driver
      use driver_mintfo_module, only: run_mintfo_driver
      use fnlo_process_common, only: nndim,flat_grid,
     $     i_momcmp_count,xratmax,ntot,nsun,nsps,nups,neps,
     $     n100,nddp,nqdp,nini,n10,n1
      implicit none

      call run_mintfo_driver(nndim,flat_grid,i_momcmp_count,
     $     xratmax,ntot,nsun,nsps,nups,neps,n100,nddp,nqdp,
     $     nini,n10,n1)
      end


      integer function drv_getpid()
      implicit none
      intrinsic getpid

      drv_getpid=getpid()
      return
      end


      double precision function sigint(xx,vegas_wgt,ifl,f)
      use driver_mintfo_module, only: sigint_impl
      use mint_module, only: ndimmax,nintegrals
      use fnlo_process_common, only: ini_fin_fks,nndim,nbody,
     $     event_momenta,p_born,virtual_over_born,calculated_born,abrv,
     $     wgt_me_born,wgt_me_real
      implicit none
      double precision xx(ndimmax),vegas_wgt,f(nintegrals)
      integer ifl

      sigint=sigint_impl(xx,vegas_wgt,ifl,f,ini_fin_fks,nndim,
     $     nbody,event_momenta,p_born,virtual_over_born,calculated_born,
     $     abrv,wgt_me_born,wgt_me_real)
      return
      end


      subroutine get_user_params(ncall,nitmax,irestart)
      use driver_mintfo_module, only: init_driver_generated_data,
     $     get_user_params_impl
      use fnlo_process_common, only: ini_fin_fks,isum_hel,
     $     multi_channel,abrv,nbody,
     $     mapconfig=>config_map,mc_hel
      use fnlo_runtime_common, only: random_offset_split
      implicit none
      integer ncall,nitmax,irestart

      call init_driver_generated_data(mapconfig)
      call get_user_params_impl(ncall,nitmax,irestart,ini_fin_fks,
     $     isum_hel,multi_channel,abrv,nbody,mc_hel,
     $     random_offset_split)
      return
      end
