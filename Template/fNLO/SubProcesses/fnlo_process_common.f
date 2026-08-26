c Process-sized shared storage for the fixed-order template.
c
c Hand-written fNLO code accesses these objects through USE association.
c Named COMMON layout is retained where generated routines still require
c that external ABI; maintained code uses module association.

      module fnlo_process_common
      use extra_weights, only: maxscales,maxPDFs,maxPDFsets,
     $     maxdynscales
      implicit none

c Generated compile-time dimensions and split-order storage.
      include 'nexternal.inc'
      include 'genps.inc'
      include 'nFKSconfigs.inc'
      include 'orders.inc'
      include 'leshouche_decl.inc'
      include 'pineappl_maxproc.inc'

      integer fnlo_maxchannels
      parameter (fnlo_maxchannels=20)

c Event slots used by real-emission and FKS subtraction kinematics.
      integer soft_counterevent,collinear_counterevent
      integer soft_collinear_counterevent,real_event
      integer first_counterevent,last_counterevent
      parameter (soft_counterevent=0,collinear_counterevent=1)
      parameter (soft_collinear_counterevent=2,real_event=3)
      parameter (first_counterevent=soft_counterevent)
      parameter (last_counterevent=soft_collinear_counterevent)

c Selected subprocess and generated process records.
      integer nfksprocess
      common /c_nfksprocess/nfksprocess

      integer fks_j_from_i(nexternal,0:nexternal)
      integer particle_type(nexternal),pdg_type(nexternal)
      common /c_fks_inc/fks_j_from_i,particle_type,pdg_type

      integer i_fks,j_fks
      common /fks_indices/i_fks,j_fks

      logical need_color_links
      common /c_need_links/need_color_links

      integer iextra_cnt,isplitorder_born,isplitorder_cnt
      common /c_extra_cnt/iextra_cnt,isplitorder_born,
     $     isplitorder_cnt

      logical split_type(nsplitorders)
      logical split_type_used(nsplitorders)
      common /c_split_type/split_type
      common /to_split_type_used/split_type_used

      logical is_aorg(nexternal)
      common /c_is_aorg/is_aorg

      double complex ans_cnt(2,nsplitorders)
      common /c_born_cnt/ans_cnt

      integer idup(nexternal,maxproc)
      integer mothup(2,nexternal,maxproc)
      integer icolup(2,nexternal,maxflow),niprocs
      common /c_leshouche_inc/idup,mothup,icolup,niprocs
      common /c_leshouche_idup_d/idup_d

      double precision subproc_pd(0:maxproc)
      integer subproc_iproc
      common /subproc/subproc_pd,subproc_iproc

      integer flavour_map(fks_configs)
      integer iproc_save(fks_configs)
      integer eto(maxproc,fks_configs),etoi(maxproc,fks_configs)
      integer maxproc_found
      common /c_flavour_map/flavour_map
      common /cproc_combination/iproc_save,eto,etoi,maxproc_found

c Generated phase-space configurations.
      double precision config_mass(-nexternal:0,lmaxconfigs,
     $     0:fks_configs)
      double precision config_width(-nexternal:0,lmaxconfigs,
     $     0:fks_configs)
      integer config_forest(2,-max_branch:-1,lmaxconfigs,
     $     0:fks_configs)
      integer config_sprop(-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer config_tprid(-max_branch:-1,lmaxconfigs,0:fks_configs)
      integer config_map(0:lmaxconfigs,0:fks_configs)
      common /c_configurations/config_mass,config_width,
     $     config_forest,config_sprop,config_tprid,config_map

      integer real_forest(2,-max_branch:-1,lmaxconfigs)
      integer real_sprop(-max_branch:-1,lmaxconfigs)
      integer real_tprid(-max_branch:-1,lmaxconfigs)
      integer real_map(0:lmaxconfigs)
      double precision real_mass(-max_branch:nexternal,lmaxconfigs)
      double precision real_width(-max_branch:-1,lmaxconfigs)
      integer real_prow(-max_branch:-1,lmaxconfigs)
      common /c_configs_inc/real_forest,real_sprop,real_tprid,
     $     real_map
      common /c_props_inc/real_mass,real_width,real_prow

      integer config_tree(2,-max_branch:-1),config_index
      common /to_itree/config_tree,config_index

      integer born_tree(2,-max_branch:-1)
      integer born_ns,born_nt,born_onebody,born_nbranch
      logical born_one_body
      common /born_trees/born_tree,born_ns,born_nt,born_onebody,
     $     born_nbranch,born_one_body

      integer iconfig0,this_config
      common /ciconfig0/iconfig0
      common /to_mconfigs/this_config

      double precision cbw_mass(-1:1,-nexternal:-1)
      double precision cbw_width(-1:1,-nexternal:-1)
      integer cbw_level_max,cbw(-nexternal:-1)
      integer cbw_level(-nexternal:-1)
      common /c_conflictingbw/cbw_mass,cbw_width,cbw_level_max,
     $     cbw,cbw_level

      double precision particle_masses(nexternal)
      double precision schannel_masses(-nexternal:nexternal)
      common /to_mass/particle_masses
      common /to_phase_space_s_channel/schannel_masses

c Born and event momenta.
      double precision p_born(0:3,nexternal-1)
      double precision p_born_l(0:3,nexternal-1)
      double precision p_born_coll(0:3,nexternal-1)
      double precision p_born_norad(0:3,nexternal-1)
      common /pborn/p_born
      common /pborn_l/p_born_l
      common /pborn_coll/p_born_coll
      common /pborn_norad/p_born_norad

c The real event and its three FKS counterevents share one layout.
      double precision event_momenta(0:3,nexternal,
     $     soft_counterevent:real_event)
      double precision event_jacobian(soft_counterevent:real_event)
      double precision event_xi(soft_counterevent:real_event)
      double precision event_y(soft_counterevent:real_event)
      double precision event_xi_hat(soft_counterevent:real_event)
      double precision event_fks_momentum(0:3,
     $     soft_counterevent:real_event)
      double precision event_xi_max(soft_counterevent:real_event)
      double precision event_xi_norm(soft_counterevent:real_event)
      double precision event_bjorken_x(2,
     $     soft_counterevent:real_event)
      double precision event_sqrt_shat(soft_counterevent:real_event)
      double precision event_shat(soft_counterevent:real_event)
      double precision ybst_til_tolab(soft_counterevent:real_event)
      double precision ybst_til_tocm(soft_counterevent:real_event)

      double precision veckn_ev,veckbarn_ev,xp0jfks
      common /cgenps_fks/veckn_ev,veckbarn_ev,xp0jfks

      complex(kind=kind(0d0)) xij_aor
      common /cxij_aor/xij_aor

      double precision xi_i_fks_fix,y_ij_fks_fix
      logical softtest,colltest
      common /cxiyfix/xi_i_fks_fix,y_ij_fks_fix
      common /sctests/softtest,colltest

      double precision tau_born_lower_bound
      double precision tau_lower_bound_resonance,tau_lower_bound
      common /ctau_lower_bound/tau_born_lower_bound,
     $     tau_lower_bound_resonance,tau_lower_bound

c Matrix-element and subtraction state.
      double precision amp2(ngraphs),jamp2(0:ncolor)
      common /to_amps/amp2,jamp2

      double precision amp_split_virt(amp_split_size)
      double precision amp_split_born_for_virt(amp_split_size)
      double precision amp_split_avv(amp_split_size)
      common /to_amp_split_virt/amp_split_virt,
     $     amp_split_born_for_virt,amp_split_avv

      double precision amp_split_wgtnstmp(amp_split_size)
      double precision amp_split_wgtwnstmpmuf(amp_split_size)
      double precision amp_split_wgtwnstmpmur(amp_split_size)
      common /to_amp_split_bsv/amp_split_wgtnstmp,
     $     amp_split_wgtwnstmpmuf,amp_split_wgtwnstmpmur

      double precision amp_split_wgtdegrem_xi(amp_split_size)
      double precision amp_split_wgtdegrem_lxi(amp_split_size)
      double precision amp_split_wgtdegrem_muf(amp_split_size)
      common /to_amp_split_deg/amp_split_wgtdegrem_xi,
     $     amp_split_wgtdegrem_lxi,amp_split_wgtdegrem_muf

      double precision amp_split_wgtpsch_p(amp_split_size)
      double precision amp_split_wgtpsch_l(amp_split_size)
      double precision amp_split_wgtpsch_d(amp_split_size)
      common /to_amp_split_dis/amp_split_wgtpsch_p,
     $     amp_split_wgtpsch_l,amp_split_wgtpsch_d

      double precision amp_split_soft(amp_split_size)
      double precision amp_split_finite_ml(amp_split_size)
      double precision amp_split_poles_fks(amp_split_size,2)
      common /to_amp_split_soft/amp_split_soft
      common /to_amp_split_finite/amp_split_finite_ml
      common /to_amp_split_poles_fks/amp_split_poles_fks

      double precision wgt_me_born,wgt_me_real
      logical calculated_born
      common /c_wgt_me_tree/wgt_me_born,wgt_me_real
      common /ccalculatedborn/calculated_born

      double precision f_b,f_nb
      double precision f_r,f_s,f_c,f_dc,f_sc,f_dsc(4)
      double precision f_pdfsch_d,f_pdfsch_p,f_pdfsch_l
      common /factor_nbody/f_b,f_nb
      common /factor_n1body/f_r,f_s,f_c,f_dc,f_sc,f_dsc
      common /factor_pdfsch/f_pdfsch_d,f_pdfsch_p,f_pdfsch_l

      double precision fkssymmetryfactor,fkssymmetryfactorborn
      double precision fkssymmetryfactordeg
      integer ngluons
      common /numberofparticles/fkssymmetryfactor,
     $     fkssymmetryfactorborn,fkssymmetryfactordeg,ngluons

      double precision delta_used,xicut_used
      double precision xiscut_used,xibsvcut_used
      common /cdelta_used/delta_used
      common /cxicut_used/xicut_used
      common /cxiscut_used/xiscut_used,xibsvcut_used

      double precision diagramsymmetryfactor,iden_comp
      common /dsymfactor/diagramsymmetryfactor
      common /c_iden_comp/iden_comp

      double precision c(0:1),gamma(0:1),gammap(0:1),beta0
      common /fks_colors/c,gamma,gammap
      common /cbeta0/beta0

      integer i_type,j_type,m_type
      common /cparticle_types/i_type,j_type,m_type

      logical nbody,nocntevents,use_evpr
      common /cnbody/nbody
      common /cnocntevents/nocntevents
      common /to_use_evpr/use_evpr

      logical fixed_order
      integer nndim,ini_fin_fks(fnlo_maxchannels)
      common /c_fixed_order/fixed_order
      common /tosigint/nndim
      common /fks_channels/ini_fin_fks

      integer isum_hel,use_cut,lbw(0:nexternal)
      logical multi_channel
      common /to_matrix/isum_hel,multi_channel
      common /to_weight/use_cut
      common /to_bw/lbw

      character*4 abrv
      integer fold,ifold_counter
      double precision virtual_over_born
      common /to_abrv/abrv
      common /cfl/fold,ifold_counter
      common /c_vob/virtual_over_born

      integer i_momcmp_count
      double precision xratmax

      double precision rndec(10)
      common /crndec/rndec

c Driver, loop-provider, and integration state.
      logical flat_grid,useitmax
      common /to_readgrid/flat_grid
      common /cuseitmax/useitmax

      integer ntot,nsun,nsps,nups,neps,n100,nddp,nqdp,nini,n10
      integer n1(0:9)
      common /ups_stats/ntot,nsun,nsps,nups,neps,n100,nddp,nqdp,
     $     nini,n10,n1

      double precision volh
      integer mc_hel,ihel
      logical fillh
      common /mc_int2/volh,mc_hel,ihel,fillh

      logical force_polecheck,polecheck_passed,cs_run,updateloop
      integer ret_code_common
      common /to_polecheck/force_polecheck,polecheck_passed
      common /to_ret_code/ret_code_common
      common /to_cs_run/cs_run
      common /to_updateloop/updateloop

      double precision qes2
      common /coupl_es/qes2

      double precision model_momenta(0:3,max_particles)
      common /momenta_pp/model_momenta

      logical is_a_j(nexternal),is_a_lp(nexternal)
      logical is_a_lm(nexternal),is_a_ph(nexternal)
      common /to_specisa/is_a_j,is_a_lp,is_a_lm,is_a_ph

      double precision etmin(nincoming+1:nexternal-1)
      double precision etmax(nincoming+1:nexternal-1)
      double precision mxxmin(nincoming+1:nexternal-1,
     $     nincoming+1:nexternal-1)
      common /to_cuts/etmin,etmax,mxxmin

      character*80 mur_id_str,muf1_id_str,muf2_id_str,qes_id_str
      character*80 temp_scale_id
      common /cscales_id_string/mur_id_str,muf1_id_str,muf2_id_str,
     $     qes_id_str
      common /ctemp_scale_id/temp_scale_id

c Analysis and histogram state.
      integer orders_tag_plot,amp_pos_plot
      common /corderstagplot/orders_tag_plot
      common /campposplot/amp_pos_plot

      double precision xsecscale_acc(maxscales,maxscales,
     $     maxdynscales)
      double precision xsecpdfr_acc(0:maxPDFs,maxPDFsets)
      common /scale_pdf_print/xsecscale_acc,xsecpdfr_acc

      integer iappl
      common /for_applgrid/iappl

c PineAPPL interface storage.
      integer appl_amp_split_size
      integer appl_qcdpower(amp_split_size)
      integer appl_qedpower(amp_split_size)
      common /appl_common_fixed/appl_amp_split_size,appl_qcdpower,
     $     appl_qedpower

      integer appl_nlumi,appl_nproc(mxpdflumi)
      integer appl_lumimap(2,max_nproc,mxpdflumi)
      common /appl_common_lumi/appl_lumimap,appl_nproc,appl_nlumi

      double precision appl_x1(4),appl_x2(4)
      double precision appl_muf2(4),appl_mur2(4),appl_qes2(4)
      double precision appl_w0(4,amp_split_size)
      double precision appl_wr(4,amp_split_size)
      double precision appl_wf(4,amp_split_size)
      double precision appl_wb(4,amp_split_size)
      integer appl_flavmap(4)
      common /appl_common_weights/appl_x1,appl_x2,appl_muf2,
     $     appl_mur2,appl_qes2,appl_w0,appl_wr,appl_wf,appl_wb,
     $     appl_flavmap

      double precision appl_q2min,appl_q2max,appl_xmin,appl_xmax
      integer appl_nq2,appl_q2order,appl_nx,appl_xorder
      common /appl_common_grid/appl_q2min,appl_q2max,appl_xmin,
     $     appl_xmax,appl_nq2,appl_q2order,appl_nx,appl_xorder

      double precision appl_www_histo,appl_norm_histo
      double precision appl_obs_histo,appl_obs_min,appl_obs_max
      double precision appl_obs_bins(0:100)
      integer appl_obs_nbins,appl_itype,appl_amp_pos,appl_obs_num
      common /appl_common_histokin/appl_www_histo,appl_norm_histo,
     $     appl_obs_histo,appl_obs_min,appl_obs_max,appl_obs_bins,
     $     appl_obs_nbins,appl_itype,appl_amp_pos,appl_obs_num

      double precision appl_event_weight,appl_vegaswgt
      common /appl_common_reco/appl_event_weight,appl_vegaswgt

      integer nh_obs,ih_obs(50)
      common /appl_ident/nh_obs,ih_obs

c Pointer-binding bridges require TARGET on their actual arguments.
      target amp_split,amp_split_cnt,idup_d
      target qes2,p_born,p_born_coll,p_born_norad
      target event_momenta,event_jacobian
      target idup,mothup,icolup,niprocs,is_aorg,amp2,jamp2
      target fks_j_from_i,particle_type,pdg_type,split_type,ans_cnt
      target split_type_used,subproc_pd,subproc_iproc
      target flavour_map,iproc_save,eto,etoi,maxproc_found
      target config_mass,config_width,config_forest,config_sprop
      target config_tprid,config_map,real_forest,real_sprop
      target real_tprid,real_map,real_mass,real_width,real_prow
      target config_tree,config_index,born_tree,born_ns,born_nt
      target born_onebody,born_nbranch,born_one_body,p_born_l
      target cbw_mass,cbw_width,cbw_level_max,cbw,cbw_level
      target particle_masses,schannel_masses
      target appl_amp_split_size,appl_qcdpower,appl_qedpower
      target appl_nproc,appl_x1,appl_x2,appl_muf2,appl_mur2
      target appl_qes2,appl_w0,appl_wr,appl_wf,appl_wb
      target appl_flavmap,appl_event_weight,appl_vegaswgt

      end module fnlo_process_common
