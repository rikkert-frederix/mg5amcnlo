c Python-generated declarations stay on this fixed-form boundary. The
c free-form implementation aliases the same live fnlo_process_common
c objects as generated matrix elements, chooser, genps, and driver.

      subroutine init_fks_singular_bridge()
      use fks_singular_module, only:
     $     initialize_fks_model_state,initialize_fks_phase_state,
     $     initialize_fks_amplitude_state,initialize_fks_config_state,
     $     initialize_fks_pineappl_state,
     $     initialize_fks_generated_state,
     $     validate_fks_singular_state
      use fnlo_process_common, only: nexternal,qes2,p_born,
     $     p_born_coll,p_born_norad,event_momenta,event_jacobian,
     $     idup_common=>idup,mothup_common=>mothup,
     $     icolup_common=>icolup,niprocs_common=>niprocs,idup_d,
     $     amp_split,amp_split_cnt,amp_split_virt,
     $     amp_split_born_for_virt,amp_split_avv,
     $     amp_split_wgtnstmp,amp_split_wgtwnstmpmuf,
     $     amp_split_wgtwnstmpmur,amp_split_wgtdegrem_xi,
     $     amp_split_wgtdegrem_lxi,amp_split_wgtdegrem_muf,
     $     amp_split_wgtpsch_p,amp_split_wgtpsch_l,
     $     amp_split_wgtpsch_d,amp_split_soft,amp_split_finite_ml,
     $     amp_split_poles_fks,fks_j_from_i,particle_type,pdg_type,
     $     split_type,ans_cnt,split_type_used,amp2,jamp2,
     $     config_mass,config_width,config_forest,config_sprop,
     $     config_tprid,config_map,real_forest,real_sprop,
     $     real_tprid,real_map,real_mass,real_width,real_prow,
     $     subproc_pd,subproc_iproc,flavour_map,iproc_save,eto,etoi,
     $     maxproc_found,is_aorg,appl_amp_split_size,appl_qcdpower,
     $     appl_qedpower,appl_nproc,appl_x1,appl_x2,appl_muf2,
     $     appl_mur2,appl_qes2,appl_w0,appl_wr,appl_wf,appl_wb,
     $     appl_flavmap,appl_event_weight,appl_vegaswgt
      implicit none
      include 'coupl.inc'

      double precision zero
      parameter (zero=0d0)
      double precision pmass(nexternal)
      logical initialized
      save initialized
      data initialized /.false./
      target g

      if (initialized) return
      call init_process_dimensions_bridge()
      call init_fks_metadata_bridge()
      include 'pmass.inc'

      call initialize_fks_model_state(g,qes2,nf,pmass)
      call initialize_fks_phase_state(p_born,p_born_coll,
     $     p_born_norad,event_momenta,event_jacobian,
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
