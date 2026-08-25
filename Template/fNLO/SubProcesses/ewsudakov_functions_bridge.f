c Generated data and the historical external ABI live here.

      subroutine init_ewsudakov_defaults_bridge()
      use ewsudakov_functions_module,
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


      subroutine init_ewsudakov_state_bridge()
      use ewsudakov_functions_module,
     &     only: initialize_ewsudakov_state
      use run_state, only: scale
      implicit none
      include 'coupl.inc'
      include '../../Source/MODEL/input.inc'
      include 'q_es.inc'
      include 'ewsudakov_haslo.inc'
      integer sud_mod,deb_settozero
      common /to_sud_mod/ sud_mod
      common /to_deb_settozero/ deb_settozero
      logical fav4,s_to_rij,rij_ge_mw,printinewsdkf
      common /to_FAV4/ fav4
      common /to_s_to_rij/ s_to_rij
      common /rij_ge_mw/ rij_ge_mw
      common /to_printinewsdkf/ printinewsdkf

      call init_ewsudakov_defaults_bridge()
      call init_process_dimensions_bridge()
      call initialize_ewsudakov_state(g,gal,mdl_mh,mdl_mt,mdl_mw,
     &     mdl_mz,mdl_ee,aewm1,mdl_gf,qes2,scale,sud_mod,fav4,
     &     s_to_rij,rij_ge_mw,printinewsdkf,deb_settozero,
     &     has_lo1,has_lo2)
      return
      end


      logical function has_lo2()
      use ewsudakov_functions_module,
     &     only: has_lo2_core => has_lo2_impl
      implicit none
      call init_ewsudakov_state_bridge()
      has_lo2=has_lo2_core()
      return
      end


      subroutine get_lo2_orders(lo2_orders)
      use ewsudakov_functions_module,
     &     only: get_lo2_orders_core => get_lo2_orders
      implicit none
      include 'orders.inc'
      integer lo2_orders(nsplitorders)
      call init_ewsudakov_state_bridge()
      call get_lo2_orders_core(lo2_orders)
      return
      end


      double complex function get_imlog(s)
      use ewsudakov_functions_module,
     &     only: get_imlog_core => get_imlog
      implicit none
      double precision s
      call init_ewsudakov_state_bridge()
      get_imlog=get_imlog_core(s)
      return
      end


      double complex function get_lsc_diag(pdglist,hels,iflist,
     &     invariants)
      use ewsudakov_functions_module,
     &     only: get_lsc_diag_core => get_lsc_diag
      implicit none
      include 'nexternal.inc'
      integer pdglist(nexternal-1),hels(nexternal-1)
      integer iflist(nexternal-1)
      double precision invariants(nexternal-1,nexternal-1)
      call init_ewsudakov_state_bridge()
      get_lsc_diag=get_lsc_diag_core(pdglist,hels,iflist,invariants)
      return
      end


      double complex function get_lsc_nondiag(pdglist,hels,iflist,
     &     invariants,ileg,pdg_old,pdg_new)
      use ewsudakov_functions_module,
     &     only: get_lsc_nondiag_core => get_lsc_nondiag
      implicit none
      include 'nexternal.inc'
      integer pdglist(nexternal-1),hels(nexternal-1)
      integer iflist(nexternal-1),ileg,pdg_old,pdg_new
      double precision invariants(nexternal-1,nexternal-1)
      call init_ewsudakov_state_bridge()
      get_lsc_nondiag=get_lsc_nondiag_core(pdglist,hels,iflist,
     &     invariants,ileg,pdg_old,pdg_new)
      return
      end


      double complex function get_ssc_c(ileg1,ileg2,pdglist,pdgp1,
     &     pdgp2,hels,iflist,invariants)
      use ewsudakov_functions_module,
     &     only: get_ssc_c_core => get_ssc_c
      implicit none
      include 'nexternal.inc'
      integer ileg1,ileg2,pdgp1,pdgp2
      integer pdglist(nexternal-1),hels(nexternal-1)
      integer iflist(nexternal-1)
      double precision invariants(nexternal-1,nexternal-1)
      call init_ewsudakov_state_bridge()
      get_ssc_c=get_ssc_c_core(ileg1,ileg2,pdglist,pdgp1,pdgp2,
     &     hels,iflist,invariants)
      return
      end


      double complex function get_ssc_n_diag(pdglist,hels,iflist,
     &     invariants)
      use ewsudakov_functions_module,
     &     only: get_ssc_n_diag_core => get_ssc_n_diag
      implicit none
      include 'nexternal.inc'
      integer pdglist(nexternal-1),hels(nexternal-1)
      integer iflist(nexternal-1)
      double precision invariants(nexternal-1,nexternal-1)
      call init_ewsudakov_state_bridge()
      get_ssc_n_diag=get_ssc_n_diag_core(pdglist,hels,iflist,
     &     invariants)
      return
      end


      double complex function get_ssc_n_nondiag_1(pdglist,hels,
     &     iflist,invariants,ileg,pdg_old,pdg_new)
      use ewsudakov_functions_module,
     &     only: get_ssc_n_nondiag_1_core => get_ssc_n_nondiag_1
      implicit none
      include 'nexternal.inc'
      integer pdglist(nexternal-1),hels(nexternal-1)
      integer iflist(nexternal-1),ileg,pdg_old,pdg_new
      double precision invariants(nexternal-1,nexternal-1)
      call init_ewsudakov_state_bridge()
      get_ssc_n_nondiag_1=get_ssc_n_nondiag_1_core(pdglist,hels,
     &     iflist,invariants,ileg,pdg_old,pdg_new)
      return
      end


      double complex function get_ssc_n_nondiag_2(pdglist,hels,
     &     iflist,invariants,ileg1,pdg_old1,pdg_new1,ileg2,pdg_old2,
     &     pdg_new2)
      use ewsudakov_functions_module,
     &     only: get_ssc_n_nondiag_2_core => get_ssc_n_nondiag_2
      implicit none
      include 'nexternal.inc'
      integer pdglist(nexternal-1),hels(nexternal-1)
      integer iflist(nexternal-1),ileg1,pdg_old1,pdg_new1
      integer ileg2,pdg_old2,pdg_new2
      double precision invariants(nexternal-1,nexternal-1)
      call init_ewsudakov_state_bridge()
      get_ssc_n_nondiag_2=get_ssc_n_nondiag_2_core(pdglist,hels,
     &     iflist,invariants,ileg1,pdg_old1,pdg_new1,ileg2,pdg_old2,
     &     pdg_new2)
      return
      end


      double complex function get_xxc_diag(pdglist,hels,iflist,
     &     invariants)
      use ewsudakov_functions_module,
     &     only: get_xxc_diag_core => get_xxc_diag
      implicit none
      include 'nexternal.inc'
      integer pdglist(nexternal-1),hels(nexternal-1)
      integer iflist(nexternal-1)
      double precision invariants(nexternal-1,nexternal-1)
      call init_ewsudakov_state_bridge()
      get_xxc_diag=get_xxc_diag_core(pdglist,hels,iflist,invariants)
      return
      end


      double complex function get_xxc_nondiag(pdglist,hels,iflist,
     &     invariants,ileg,pdg_old,pdg_new)
      use ewsudakov_functions_module,
     &     only: get_xxc_nondiag_core => get_xxc_nondiag
      implicit none
      include 'nexternal.inc'
      integer pdglist(nexternal-1),hels(nexternal-1)
      integer iflist(nexternal-1),ileg,pdg_old,pdg_new
      double precision invariants(nexternal-1,nexternal-1)
      call init_ewsudakov_state_bridge()
      get_xxc_nondiag=get_xxc_nondiag_core(pdglist,hels,iflist,
     &     invariants,ileg,pdg_old,pdg_new)
      return
      end


      double complex function bigL(s)
      use ewsudakov_functions_module, only: bigl_core => bigl
      implicit none
      double precision s
      call init_ewsudakov_state_bridge()
      bigL=bigl_core(s)
      return
      end


      double complex function smallL(s)
      use ewsudakov_functions_module, only: smalll_core => smalll
      implicit none
      double precision s
      call init_ewsudakov_state_bridge()
      smallL=smalll_core(s)
      return
      end


      double complex function log_a_over_b_sing(a,b)
      use ewsudakov_functions_module,
     &     only: log_a_over_b_sing_core => log_a_over_b_sing
      implicit none
      double precision a,b
      call init_ewsudakov_state_bridge()
      log_a_over_b_sing=log_a_over_b_sing_core(a,b)
      return
      end


      double complex function smallL_a_over_b_sing(a,b)
      use ewsudakov_functions_module,
     &     only: smalll_ab_core => smalll_a_over_b_sing
      implicit none
      double precision a,b
      call init_ewsudakov_state_bridge()
      smallL_a_over_b_sing=smalll_ab_core(a,b)
      return
      end


      double complex function bigL_a_over_b_sing(a,b)
      use ewsudakov_functions_module,
     &     only: bigl_ab_core => bigl_a_over_b_sing
      implicit none
      double precision a,b
      call init_ewsudakov_state_bridge()
      bigL_a_over_b_sing=bigl_ab_core(a,b)
      return
      end


      double complex function bigLem(s,m2k)
      use ewsudakov_functions_module, only: biglem_core => biglem
      implicit none
      double precision s,m2k
      call init_ewsudakov_state_bridge()
      bigLem=biglem_core(s,m2k)
      return
      end


      double complex function smallLem(m2k)
      use ewsudakov_functions_module, only: smalllem_core => smalllem
      implicit none
      double precision m2k
      call init_ewsudakov_state_bridge()
      smallLem=smalllem_core(m2k)
      return
      end


      double complex function dZemAA_logs()
      use ewsudakov_functions_module,
     &     only: dzemaa_logs_core => dzemaa_logs
      implicit none
      call init_ewsudakov_state_bridge()
      dZemAA_logs=dzemaa_logs_core()
      return
      end


      double complex function sdk_chargesq(pdg,hel,ifsign)
      use ewsudakov_functions_module,
     &     only: sdk_chargesq_core => sdk_chargesq
      implicit none
      integer pdg,hel,ifsign
      call init_ewsudakov_state_bridge()
      sdk_chargesq=sdk_chargesq_core(pdg,hel,ifsign)
      return
      end


      double complex function sdk_charge(pdg,hel,ifsign)
      use ewsudakov_functions_module,
     &     only: sdk_charge_core => sdk_charge
      implicit none
      integer pdg,hel,ifsign
      call init_ewsudakov_state_bridge()
      sdk_charge=sdk_charge_core(pdg,hel,ifsign)
      return
      end


      double complex function sdk_tpm(pdg,hel,ifsign,pdgp)
      use ewsudakov_functions_module, only: sdk_tpm_core => sdk_tpm
      implicit none
      integer pdg,hel,ifsign,pdgp
      call init_ewsudakov_state_bridge()
      sdk_tpm=sdk_tpm_core(pdg,hel,ifsign,pdgp)
      return
      end


      double complex function sdk_t3_diag(pdg,hel,ifsign)
      use ewsudakov_functions_module,
     &     only: sdk_t3_diag_core => sdk_t3_diag
      implicit none
      integer pdg,hel,ifsign
      call init_ewsudakov_state_bridge()
      sdk_t3_diag=sdk_t3_diag_core(pdg,hel,ifsign)
      return
      end


      double complex function sdk_yo2_diag(pdg,hel,ifsign)
      use ewsudakov_functions_module,
     &     only: sdk_yo2_diag_core => sdk_yo2_diag
      implicit none
      integer pdg,hel,ifsign
      call init_ewsudakov_state_bridge()
      sdk_yo2_diag=sdk_yo2_diag_core(pdg,hel,ifsign)
      return
      end


      double complex function sdk_iz_diag(pdg,hel,ifsign)
      use ewsudakov_functions_module,
     &     only: sdk_iz_diag_core => sdk_iz_diag
      implicit none
      integer pdg,hel,ifsign
      call init_ewsudakov_state_bridge()
      sdk_iz_diag=sdk_iz_diag_core(pdg,hel,ifsign)
      return
      end


      double complex function sdk_iz_nondiag(pdg,hel,ifsign)
      use ewsudakov_functions_module,
     &     only: sdk_iz_nondiag_core => sdk_iz_nondiag
      implicit none
      integer pdg,hel,ifsign
      call init_ewsudakov_state_bridge()
      sdk_iz_nondiag=sdk_iz_nondiag_core(pdg,hel,ifsign)
      return
      end


      double complex function sdk_ia_diag(pdg,hel,ifsign)
      use ewsudakov_functions_module,
     &     only: sdk_ia_diag_core => sdk_ia_diag
      implicit none
      integer pdg,hel,ifsign
      call init_ewsudakov_state_bridge()
      sdk_ia_diag=sdk_ia_diag_core(pdg,hel,ifsign)
      return
      end


      double complex function sdk_iz2_diag(pdg,hel,ifsign)
      use ewsudakov_functions_module,
     &     only: sdk_iz2_diag_core => sdk_iz2_diag
      implicit none
      integer pdg,hel,ifsign
      call init_ewsudakov_state_bridge()
      sdk_iz2_diag=sdk_iz2_diag_core(pdg,hel,ifsign)
      return
      end


      double complex function sdk_cew_diag(pdg,hel,ifsign)
      use ewsudakov_functions_module,
     &     only: sdk_cew_diag_core => sdk_cew_diag
      implicit none
      integer pdg,hel,ifsign
      call init_ewsudakov_state_bridge()
      sdk_cew_diag=sdk_cew_diag_core(pdg,hel,ifsign)
      return
      end


      double complex function sdk_cew_nondiag()
      use ewsudakov_functions_module,
     &     only: sdk_cew_nondiag_core => sdk_cew_nondiag
      implicit none
      call init_ewsudakov_state_bridge()
      sdk_cew_nondiag=sdk_cew_nondiag_core()
      return
      end


      double complex function sdk_betaew_diag(pdg)
      use ewsudakov_functions_module,
     &     only: sdk_betaew_diag_core => sdk_betaew_diag
      implicit none
      integer pdg
      call init_ewsudakov_state_bridge()
      sdk_betaew_diag=sdk_betaew_diag_core(pdg)
      return
      end


      double complex function sdk_betaew_nondiag()
      use ewsudakov_functions_module,
     &     only: sdk_betaew_nondiag_core => sdk_betaew_nondiag
      implicit none
      call init_ewsudakov_state_bridge()
      sdk_betaew_nondiag=sdk_betaew_nondiag_core()
      return
      end


      subroutine sdk_get_invariants(p,iflist,invariants)
      use ewsudakov_functions_module,
     &     only: sdk_get_invariants_core => sdk_get_invariants
      implicit none
      include 'nexternal.inc'
      double precision p(0:3,nexternal-1)
      integer iflist(nexternal-1)
      double precision invariants(nexternal-1,nexternal-1)
      call init_ewsudakov_state_bridge()
      call sdk_get_invariants_core(p,iflist,invariants)
      return
      end


      subroutine sdk_test_functions()
      use ewsudakov_functions_module,
     &     only: sdk_test_functions_core => sdk_test_functions
      implicit none
      call init_ewsudakov_state_bridge()
      call sdk_test_functions_core()
      return
      end


      double precision function get_isopart_mass_from_id(pdg)
      use ewsudakov_functions_module,
     &     only: get_isopart_mass_core => get_isopart_mass_from_id
      implicit none
      integer pdg
      call init_ewsudakov_state_bridge()
      get_isopart_mass_from_id=get_isopart_mass_core(pdg)
      return
      end


      subroutine get_par_ren_alphamz(invariants)
      use ewsudakov_functions_module,
     &     only: get_par_ren_alphamz_core => get_par_ren_alphamz
      implicit none
      include 'nexternal.inc'
      include 'orders.inc'
      double precision invariants(nexternal-1,nexternal-1)
      integer imaxpara
      parameter (imaxpara=6)
      double complex amp_split_ewsud_der(amp_split_size,imaxpara)
      double complex amp_split_ewsud_der2(amp_split_size,imaxpara)
      double complex amp_split_ewsud(amp_split_size)
      common /to_amp_split_ewsud_der/ amp_split_ewsud_der
      common /to_amp_split_ewsud_der2/ amp_split_ewsud_der2
      common /to_amp_split_ewsud/ amp_split_ewsud
      call init_ewsudakov_state_bridge()
      call get_par_ren_alphamz_core(invariants,amp_split_ewsud_der,
     &     amp_split_ewsud_der2,amp_split_ewsud)
      return
      end


      subroutine get_par_ren_gmu(invariants)
      use ewsudakov_functions_module,
     &     only: get_par_ren_gmu_core => get_par_ren_gmu
      implicit none
      include 'nexternal.inc'
      include 'orders.inc'
      double precision invariants(nexternal-1,nexternal-1)
      integer imaxpara
      parameter (imaxpara=6)
      double complex amp_split_ewsud_der(amp_split_size,imaxpara)
      double complex amp_split_ewsud_der2(amp_split_size,imaxpara)
      double complex amp_split_ewsud(amp_split_size)
      common /to_amp_split_ewsud_der/ amp_split_ewsud_der
      common /to_amp_split_ewsud_der2/ amp_split_ewsud_der2
      common /to_amp_split_ewsud/ amp_split_ewsud
      call init_ewsudakov_state_bridge()
      call get_par_ren_gmu_core(invariants,amp_split_ewsud_der,
     &     amp_split_ewsud_der2,amp_split_ewsud)
      return
      end


      double complex function get_qcd_lo2(pdglist,hels,iflist,
     &     invariants,iamp)
      use ewsudakov_functions_module,
     &     only: get_qcd_lo2_core => get_qcd_lo2
      implicit none
      include 'nexternal.inc'
      integer pdglist(nexternal-1),hels(nexternal-1)
      integer iflist(nexternal-1),iamp
      double precision invariants(nexternal-1,nexternal-1)
      call init_ewsudakov_state_bridge()
      get_qcd_lo2=get_qcd_lo2_core(pdglist,hels,iflist,
     &     invariants,iamp)
      return
      end
