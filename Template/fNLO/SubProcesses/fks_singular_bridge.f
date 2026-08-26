c Python-generated declarations stay on this fixed-form boundary.  The
c free-form implementation aliases the same live COMMON objects as the
c generated matrix elements, chooser, genps, and fixed-order driver.

      subroutine init_fks_singular_bridge()
      use fks_singular_module, only:
     $     initialize_fks_model_state,initialize_fks_phase_state,
     $     initialize_fks_amplitude_state,initialize_fks_config_state,
     $     initialize_fks_pineappl_state,
     $     initialize_fks_generated_state,
     $     validate_fks_singular_state
      implicit none
      include 'nexternal.inc'
      include 'genps.inc'
      include 'nFKSconfigs.inc'
      include 'coupl.inc'
      include 'q_es.inc'
      include 'pineappl_common.inc'
      include 'leshouche_decl.inc'

      double precision zero
      parameter (zero=0d0)
      double precision pmass(nexternal)
      logical initialized
      save initialized
      data initialized /.false./

      double precision p_born(0:3,nexternal-1)
      double precision p_born_coll(0:3,nexternal-1)
      double precision p_born_norad(0:3,nexternal-1)
      double precision p_ev(0:3,nexternal)
      target p_born,p_born_coll,p_born_norad,p_ev
      common /pborn/p_born
      common /pborn_coll/p_born_coll
      common /pborn_norad/p_born_norad
      common /pev/p_ev

      double precision p1_cnt(0:3,nexternal,0:2)
      double precision wgt_cnt(0:2),pswgt_cnt(0:2)
      double precision jac_cnt(0:2)
      target p1_cnt,wgt_cnt,pswgt_cnt,jac_cnt
      common /counterevnts/p1_cnt,wgt_cnt,pswgt_cnt,jac_cnt

      integer idup_common(nexternal,maxproc)
      integer mothup_common(2,nexternal,maxproc)
      integer icolup_common(2,nexternal,maxflow)
      integer niprocs_common
      target idup_common,mothup_common,icolup_common,niprocs_common
      common /c_leshouche_inc/idup_common,mothup_common,
     $     icolup_common,niprocs_common
      target idup_d
      common /c_leshouche_idup_d/idup_d

      double precision amp_split_virt(amp_split_size)
      double precision amp_split_born_for_virt(amp_split_size)
      double precision amp_split_avv(amp_split_size)
      target amp_split_virt,amp_split_born_for_virt,amp_split_avv
      common /to_amp_split_virt/amp_split_virt,
     $     amp_split_born_for_virt,amp_split_avv

      double precision amp_split_wgtnstmp(amp_split_size)
      double precision amp_split_wgtwnstmpmuf(amp_split_size)
      double precision amp_split_wgtwnstmpmur(amp_split_size)
      target amp_split_wgtnstmp,amp_split_wgtwnstmpmuf
      target amp_split_wgtwnstmpmur
      common /to_amp_split_bsv/amp_split_wgtnstmp,
     $     amp_split_wgtwnstmpmuf,amp_split_wgtwnstmpmur

      double precision amp_split_wgtdegrem_xi(amp_split_size)
      double precision amp_split_wgtdegrem_lxi(amp_split_size)
      double precision amp_split_wgtdegrem_muf(amp_split_size)
      target amp_split_wgtdegrem_xi,amp_split_wgtdegrem_lxi
      target amp_split_wgtdegrem_muf
      common /to_amp_split_deg/amp_split_wgtdegrem_xi,
     $     amp_split_wgtdegrem_lxi,amp_split_wgtdegrem_muf

      double precision amp_split_wgtpsch_p(amp_split_size)
      double precision amp_split_wgtpsch_l(amp_split_size)
      double precision amp_split_wgtpsch_d(amp_split_size)
      target amp_split_wgtpsch_p,amp_split_wgtpsch_l
      target amp_split_wgtpsch_d
      common /to_amp_split_dis/amp_split_wgtpsch_p,
     $     amp_split_wgtpsch_l,amp_split_wgtpsch_d

      double precision amp_split_soft(amp_split_size)
      double precision amp_split_finite_ml(amp_split_size)
      double precision amp_split_poles_fks(amp_split_size,2)
      target amp_split_soft,amp_split_finite_ml,amp_split_poles_fks
      common /to_amp_split_soft/amp_split_soft
      common /to_amp_split_finite/amp_split_finite_ml
      common /to_amp_split_poles_fks/amp_split_poles_fks

      integer fks_j_from_i(nexternal,0:nexternal)
      integer particle_type(nexternal),pdg_type(nexternal)
      logical split_type(nsplitorders)
      logical split_type_used(nsplitorders)
      double complex ans_cnt(2,nsplitorders)
      target fks_j_from_i,particle_type,pdg_type
      target split_type,split_type_used,ans_cnt
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type
      common /c_split_type/split_type
      common /c_born_cnt/ans_cnt
      common /to_split_type_used/split_type_used

      double precision amp2(ngraphs),jamp2(0:ncolor)
      target amp2,jamp2
      common /to_amps/amp2,jamp2

      double precision config_mass(-nexternal:0,lmaxconfigs,
     $     0:fks_configs)
      double precision config_width(-nexternal:0,lmaxconfigs,
     $     0:fks_configs)
      integer config_forest(2,-max_branch:-1,lmaxconfigs,
     $     0:fks_configs)
      integer config_sprop(-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer config_tprid(-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer config_map(0:lmaxconfigs,0:fks_configs)
      target config_mass,config_width,config_forest
      target config_sprop,config_tprid,config_map
      common /c_configurations/config_mass,config_width,
     $     config_forest,config_sprop,config_tprid,config_map

      integer real_forest(2,-max_branch:-1,lmaxconfigs)
      integer real_sprop(-max_branch:-1,lmaxconfigs)
      integer real_tprid(-max_branch:-1,lmaxconfigs)
      integer real_map(0:lmaxconfigs)
      double precision real_mass(-max_branch:nexternal,lmaxconfigs)
      double precision real_width(-max_branch:-1,lmaxconfigs)
      integer real_prow(-max_branch:-1,lmaxconfigs)
      target real_forest,real_sprop,real_tprid,real_map
      target real_mass,real_width,real_prow
      common /c_configs_inc/real_forest,real_sprop,real_tprid,real_map
      common /c_props_inc/real_mass,real_width,real_prow

      double precision subproc_pd(0:maxproc)
      integer subproc_iproc
      target subproc_pd,subproc_iproc
      common /subproc/subproc_pd,subproc_iproc

      integer flavour_map(fks_configs),iproc_save(fks_configs)
      integer eto(maxproc,fks_configs),etoi(maxproc,fks_configs)
      integer maxproc_found
      target flavour_map,iproc_save,eto,etoi,maxproc_found
      common /c_flavour_map/flavour_map
      common /cproc_combination/iproc_save,eto,etoi,maxproc_found

      logical is_aorg(nexternal)
      target is_aorg
      common /c_is_aorg/is_aorg

      target g,qes2
      target amp_split,amp_split_cnt
      target appl_amp_split_size,appl_qcdpower,appl_qedpower
      target appl_nproc,appl_x1,appl_x2,appl_muf2,appl_mur2
      target appl_qes2,appl_w0,appl_wr,appl_wf,appl_wb
      target appl_flavmap,appl_event_weight,appl_vegaswgt

      if (initialized) return
      call init_process_dimensions_bridge()
      call init_fks_metadata_bridge()
      include 'pmass.inc'

      call initialize_fks_model_state(g,qes2,nf,pmass)
      call initialize_fks_phase_state(p_born,p_born_coll,
     $     p_born_norad,p_ev,p1_cnt,wgt_cnt,pswgt_cnt,jac_cnt,
     $     idup_common,mothup_common,icolup_common,niprocs_common,
     $     is_aorg,amp2,jamp2,subproc_pd,subproc_iproc,flavour_map,
     $     iproc_save,eto,etoi,maxproc_found)
      call initialize_fks_amplitude_state(amp_split,amp_split_cnt,
     $     amp_split_virt,amp_split_born_for_virt,amp_split_avv,
     $     amp_split_wgtnstmp,amp_split_wgtwnstmpmuf,
     $     amp_split_wgtwnstmpmur,amp_split_wgtdegrem_xi,
     $     amp_split_wgtdegrem_lxi,amp_split_wgtdegrem_muf,
     $     amp_split_wgtpsch_p,amp_split_wgtpsch_l,
     $     amp_split_wgtpsch_d,amp_split_soft,amp_split_finite_ml,
     $     amp_split_poles_fks,fks_j_from_i,particle_type,pdg_type,
     $     split_type,ans_cnt,split_type_used,idup_d)
      call initialize_fks_config_state(config_mass,config_width,
     $     config_forest,config_sprop,config_tprid,config_map,
     $     real_forest,real_sprop,real_tprid,real_map,real_mass,
     $     real_width,real_prow)
      call initialize_fks_pineappl_state(appl_amp_split_size,
     $     appl_qcdpower,appl_qedpower,appl_nproc,appl_x1,appl_x2,
     $     appl_muf2,appl_mur2,appl_qes2,appl_w0,appl_wr,
     $     appl_wf,appl_wb,appl_flavmap,appl_event_weight,
     $     appl_vegaswgt)
      call init_fks_singular_born_config_bridge()
      call validate_fks_singular_state()
      initialized=.true.
      return
      end

      subroutine init_fks_singular_born_config_bridge()
      use fks_singular_module, only: initialize_fks_generated_state
      implicit none
      include 'nexternal.inc'
      include 'maxconfigs.inc'
      include 'born_conf.inc'
      double precision born_mass(-nexternal:0,lmaxconfigs)
      double precision born_width(-nexternal:0,lmaxconfigs)

      born_mass=0d0
      born_width=0d0
c born_props.inc writes PMASS and PWIDTH by those exact names.
      call fill_fks_born_props_bridge(born_mass,born_width)
      call initialize_fks_generated_state(max_branchb_used,
     $     lmaxconfigsb_used,iforest,sprop,tprid,mapconfig,
     $     born_mass,born_width)
      return
      end


      subroutine fill_fks_born_props_bridge(pmass,pwidth)
      implicit none
      include 'nexternal.inc'
      include 'maxconfigs.inc'
      include 'coupl.inc'
      double precision zero
      parameter (zero=0d0)
      double precision pmass(-nexternal:0,lmaxconfigs)
      double precision pwidth(-nexternal:0,lmaxconfigs)
      include 'born_props.inc'
      return
      end

c Generated Born helicity code still calls this random-number entry point.
      double precision function ran2()
      use fks_singular_module, only: implementation => ran2
      implicit none
      call init_fks_singular_bridge()
      ran2=implementation()
      end
