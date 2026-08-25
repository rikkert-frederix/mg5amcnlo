c Generated dimensions and mutable COMMON boundary for the non-EW dummy.

      subroutine init_ewsudakov_defaults_bridge()
      use ewsudakov_dummy_module,
     &     only: initialize_ewsudakov_defaults
      implicit none
      integer sud_mod
      common /to_sud_mod/ sud_mod
      logical sud_filter_hel,sud_mc_hel,fav4,s_to_rij,cs_run
      logical rij_ge_mw
      common /to_filter_hel/ sud_filter_hel
      common /to_mc_hel/ sud_mc_hel
      common /to_FAV4/ fav4
      common /to_s_to_rij/ s_to_rij
      common /to_cs_run/ cs_run
      common /rij_ge_mw/ rij_ge_mw

      call initialize_ewsudakov_defaults(sud_mod,sud_filter_hel,
     &     sud_mc_hel,fav4,s_to_rij,cs_run,rij_ge_mw)
      return
      end

      subroutine get_lo2_orders(lo2_orders)
      use ewsudakov_dummy_module,
     &     only: get_lo2_orders_impl
      implicit none
      include 'orders.inc'
      integer lo2_orders(nsplitorders)

      call get_lo2_orders_impl(nsplitorders,born_orders,qcd_pos,
     &     qed_pos,lo2_orders)
      return
      end


      subroutine sdk_get_invariants(p,iflist,invariants)
      use ewsudakov_dummy_module,
     &     only: sdk_get_invariants_impl
      implicit none
      include 'nexternal.inc'
      include 'coupl.inc'
      double precision p(0:3,nexternal-1)
      integer iflist(nexternal-1)
      double precision invariants(nexternal-1,nexternal-1)
      logical rij_ge_mw
      common /rij_ge_mw/ rij_ge_mw

      call sdk_get_invariants_impl(nexternal-1,p,iflist,mdl_mw,
     &     rij_ge_mw,invariants)
      return
      end


      subroutine ewsudakov_f77(p_born_in,gstr_in,results)
      use ewsudakov_dummy_module,
     &     only: ewsudakov_f77_impl
      implicit none
      include 'nexternal.inc'
      double precision p_born_in(0:3,nexternal-1)
      double precision gstr_in,results(6)
      double precision p_born(0:3,nexternal-1)
      common /pborn/ p_born

      integer sud_mod,nfksprocess
      common /to_sud_mod/ sud_mod
      common /c_nfksprocess/ nfksprocess
      logical sud_mc_hel
      common /to_mc_hel/ sud_mc_hel
      logical s_to_rij,rij_ge_mw
      common /to_s_to_rij/ s_to_rij
      common /rij_ge_mw/ rij_ge_mw

      call ewsudakov_f77_impl(nexternal,p_born_in,gstr_in,results,
     &     p_born,sud_mod,nfksprocess,sud_mc_hel,s_to_rij,rij_ge_mw)
      return
      end


      subroutine ewsud_dummy_set_coupling(gstr)
      implicit none
      include 'coupl.inc'
      double precision gstr

      g=gstr
      call update_as_param()
      return
      end


      subroutine ewsud_dummy_evaluate(nexternal,
     &     p_born,wgt_born,ewsud_lsc,ewsud_ssc,ewsud_xxc,ewsud_par)
      implicit none
      include 'orders.inc'
      integer nexternal
      double precision p_born(0:3,nexternal-1),wgt_born,born
      double complex ewsud_lsc,ewsud_ssc,ewsud_xxc,ewsud_par
      double complex amp_split_ewsud_lsc(amp_split_size)
      common /to_amp_ewsud_lsc/ amp_split_ewsud_lsc
      double complex amp_split_ewsud_ssc(amp_split_size)
      common /to_amp_ewsud_ssc/ amp_split_ewsud_ssc
      double complex amp_split_ewsud_xxc(amp_split_size)
      common /to_amp_ewsud_xxc/ amp_split_ewsud_xxc
      double complex amp_split_ewsud_par(amp_split_size)
      common /to_amp_ewsud_par/ amp_split_ewsud_par

      call sborn(p_born,born)
      wgt_born=amp_split(1)
      call sudakov_wrapper(p_born)
      ewsud_lsc=amp_split_ewsud_lsc(1)
      ewsud_ssc=amp_split_ewsud_ssc(1)
      ewsud_xxc=amp_split_ewsud_xxc(1)
      ewsud_par=amp_split_ewsud_par(1)
      return
      end
