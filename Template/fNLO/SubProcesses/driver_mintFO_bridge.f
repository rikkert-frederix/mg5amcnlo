      program driver
      use driver_mintfo_module, only: run_mintfo_driver
      implicit none
      integer nndim
      common /tosigint/nndim
      logical flat_grid
      common /to_readgrid/flat_grid
      integer i_momcmp_count
      double precision xratmax
      common /ccheckcnt/i_momcmp_count,xratmax
      character(len=4) abrv
      common /to_abrv/abrv
      integer ntot,nsun,nsps,nups,neps,n100,nddp,nqdp,nini,n10
      integer n1(0:9)
      common /ups_stats/ntot,nsun,nsps,nups,neps,n100,nddp,nqdp,
     $     nini,n10,n1
      logical useitmax
      common /cuseitmax/useitmax

      call run_mintfo_driver(nndim,flat_grid,i_momcmp_count,
     $     xratmax,abrv,ntot,nsun,nsps,nups,neps,n100,nddp,nqdp,
     $     nini,n10,n1,useitmax)
      end


      integer function drv_getpid()
      implicit none
      intrinsic getpid

      drv_getpid=getpid()
      return
      end


      double precision function sigint(xx,vegas_wgt,ifl,f)
      use driver_mintfo_module, only: sigint_impl
      use mint_module, only: ndimmax,nintegrals,maxchannels
      implicit none
      include 'nexternal.inc'
      double precision xx(ndimmax),vegas_wgt,f(nintegrals)
      integer ifl
      integer ini_fin_fks(maxchannels)
      common /fks_channels/ini_fin_fks
      integer nndim
      common /tosigint/nndim
      logical nbody
      common /cnbody/nbody
      double precision p1_cnt(0:3,nexternal,-2:2),wgt_cnt(-2:2)
      double precision pswgt_cnt(-2:2),jac_cnt(-2:2)
      common /counterevnts/p1_cnt,wgt_cnt,pswgt_cnt,jac_cnt
      double precision p_born(0:3,nexternal-1)
      common /pborn/p_born
      double precision virtual_over_born
      common /c_vob/virtual_over_born
      logical calculated_born
      common /ccalculatedborn/calculated_born
      character(len=4) abrv
      common /to_abrv/abrv
      double precision wgt_me_born,wgt_me_real
      common /c_wgt_me_tree/wgt_me_born,wgt_me_real
      integer fold,ifold_counter
      common /cfl/fold,ifold_counter
      logical use_evpr
      common /to_use_evpr/use_evpr

      sigint=sigint_impl(xx,vegas_wgt,ifl,f,ini_fin_fks,nndim,
     $     nbody,p1_cnt,p_born,virtual_over_born,calculated_born,
     $     abrv,wgt_me_born,wgt_me_real,fold,use_evpr)
      return
      end


      subroutine update_fks_dir(nfks)
      use driver_mintfo_module, only: update_fks_dir_impl
      implicit none
      integer nfks,nfksprocess
      common /c_nfksprocess/nfksprocess

      nfksprocess=nfks
      call update_fks_dir_impl()
      return
      end


      subroutine setup_ini_fin_fks_map(ini_fin_fks_map)
      use driver_mintfo_module, only: setup_ini_fin_fks_map_impl
      implicit none
      include 'nFKSconfigs.inc'
      integer ini_fin_fks_map(0:2,0:fks_configs)

      call setup_ini_fin_fks_map_impl(ini_fin_fks_map)
      return
      end


      subroutine get_born_nfksprocess(nfks_in,nfks_out)
      use driver_mintfo_module, only: get_born_nfksprocess_impl
      implicit none
      integer nfks_in,nfks_out

      call get_born_nfksprocess_impl(nfks_in,nfks_out)
      return
      end


      subroutine update_vegas_x(xx,x)
      use driver_mintfo_module, only: update_vegas_x_impl
      use mint_module, only: ndimmax
      implicit none
      double precision xx(ndimmax),x(99)
      integer nndim
      common /tosigint/nndim
      character(len=4) abrv
      common /to_abrv/abrv

      call update_vegas_x_impl(xx,x,nndim,abrv)
      return
      end


      subroutine get_user_params(ncall,nitmax,irestart)
      use driver_mintfo_module, only: init_driver_generated_data,
     $     get_user_params_impl
      use mint_module, only: maxchannels
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'nFKSconfigs.inc'
      include 'has_ewsudakov.inc'
      integer ncall,nitmax,irestart
      integer ini_fin_fks(maxchannels)
      common /fks_channels/ini_fin_fks
      integer isum_hel
      logical multi_channel
      common /to_matrix/isum_hel,multi_channel
      integer use_cut
      common /to_weight/use_cut
      integer lbw(0:nexternal)
      common /to_bw/lbw
      character(len=4) abrv
      common /to_abrv/abrv
      logical nbody
      common /cnbody/nbody
      double precision pmass(-nexternal:0,lmaxconfigs,0:fks_configs)
      double precision pwidth(-nexternal:0,lmaxconfigs,0:fks_configs)
      integer iforest(2,-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer sprop(-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer tprid(-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer mapconfig(0:lmaxconfigs,0:fks_configs)
      common /c_configurations/pmass,pwidth,iforest,sprop,tprid,
     $     mapconfig
      double precision volh
      integer mc_hel,ihel
      logical fillh
      common /mc_int2/volh,mc_hel,ihel,fillh
      integer random_offset_split
      common /c_random_offset_split/random_offset_split

      call init_driver_generated_data(mapconfig,has_ewsudakov)
      call get_user_params_impl(ncall,nitmax,irestart,ini_fin_fks,
     $     isum_hel,multi_channel,use_cut,lbw,abrv,nbody,mc_hel,
     $     random_offset_split)
      return
      end


      subroutine drv_soft_collinear()
      implicit none
      interface
      subroutine compute_soft_collinear_counter_term()
      end subroutine compute_soft_collinear_counter_term
      end interface

      call compute_soft_collinear_counter_term()
      return
      end
