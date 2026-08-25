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

      double precision p1_cnt(0:3,nexternal,-2:2)
      double precision wgt_cnt(-2:2),pswgt_cnt(-2:2)
      double precision jac_cnt(-2:2)
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

      double precision amp_split_6to5f(amp_split_size)
      double precision amp_split_6to5f_muf(amp_split_size)
      double precision amp_split_6to5f_mur(amp_split_size)
      target amp_split_6to5f,amp_split_6to5f_muf,
     $     amp_split_6to5f_mur
      common /to_amp_split_6to5f/amp_split_6to5f,
     $     amp_split_6to5f_muf,amp_split_6to5f_mur

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

      target g,mdl_mt,qes2
      target amp_split,amp_split_cnt
      target appl_amp_split_size,appl_qcdpower,appl_qedpower
      target appl_nproc,appl_x1,appl_x2,appl_muf2,appl_mur2
      target appl_qes2,appl_w0,appl_wr,appl_wf,appl_wb
      target appl_flavmap,appl_event_weight,appl_vegaswgt

      if (initialized) return
      call init_process_dimensions_bridge()
      call init_fks_metadata_bridge()
      include 'pmass.inc'

      call initialize_fks_model_state(g,qes2,mdl_mt,nf,pmass)
      call initialize_fks_phase_state(p_born,p_born_coll,
     $     p_born_norad,p_ev,p1_cnt,wgt_cnt,pswgt_cnt,jac_cnt,
     $     idup_common,mothup_common,icolup_common,niprocs_common,
     $     is_aorg,amp2,jamp2,subproc_pd,subproc_iproc,flavour_map,
     $     iproc_save,eto,etoi,maxproc_found)
      call initialize_fks_amplitude_state(amp_split,amp_split_cnt,
     $     amp_split_6to5f,amp_split_6to5f_muf,
     $     amp_split_6to5f_mur,
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
      include 'maxparticles.inc'
      include 'maxconfigs.inc'
      include 'coupl.inc'
      integer i,j
      double precision zero
      parameter (zero=0d0)
      include 'born_conf.inc'
      double precision born_mass(-nexternal:0,lmaxconfigs)
      double precision born_width(-nexternal:0,lmaxconfigs)
      integer born_pow(-nexternal:0,lmaxconfigs)

      born_mass=0d0
      born_width=0d0
      born_pow=0
c born_props.inc writes PMASS, PWIDTH and POW by those exact names.
      call fill_fks_born_props_bridge(born_mass,born_width,born_pow)
      call initialize_fks_generated_state(max_branchb_used,
     $     lmaxconfigsb_used,iforest,sprop,tprid,mapconfig,
     $     born_mass,born_width)
      return
      end


      subroutine fill_fks_born_props_bridge(pmass,pwidth,pow)
      implicit none
      include 'nexternal.inc'
      include 'maxconfigs.inc'
      include 'coupl.inc'
      double precision zero
      parameter (zero=0d0)
      double precision pmass(-nexternal:0,lmaxconfigs)
      double precision pwidth(-nexternal:0,lmaxconfigs)
      integer pow(-nexternal:0,lmaxconfigs)
      include 'born_props.inc'
      return
      end


c This historical entry point exceeds the Fortran 95 identifier limit.  It
c stays at the compiler-extension boundary while the module implementation
c uses the standard-conforming shorter name.
      subroutine compute_soft_collinear_counter_term()
      use fks_singular_module, only:
     $     implementation => compute_soft_collinear_ct_impl
      implicit none
      call init_fks_singular_bridge()
      call implementation()
      return
      end
c Compatibility entry points for generated fixed-form callers.  The
c implementation lives in fks_singular_module; generated declarations and
c historical COMMON storage are initialized by init_fks_singular_bridge.

      subroutine compute_born()
      use fks_singular_module, only: implementation => compute_born
      implicit none
      call init_fks_singular_bridge()
      call implementation()
      end subroutine compute_born

      subroutine compute_6to5flav_cnt()
      use fks_singular_module, only: implementation => compute_6to5flav_cnt
      implicit none
      call init_fks_singular_bridge()
      call implementation()
      end subroutine compute_6to5flav_cnt

      subroutine compute_nbody_noborn()
      use fks_singular_module, only: implementation => compute_nbody_noborn
      implicit none
      call init_fks_singular_bridge()
      call implementation()
      end subroutine compute_nbody_noborn

      subroutine compute_real_emission(p)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => compute_real_emission
      implicit none
      double precision :: p(0:3,nexternal)
      call init_fks_singular_bridge()
      call implementation(p)
      end subroutine compute_real_emission

      subroutine compute_soft_counter_term()
      use fks_singular_module, only: implementation => compute_soft_counter_term
      implicit none
      call init_fks_singular_bridge()
      call implementation()
      end subroutine compute_soft_counter_term

      subroutine compute_collinear_counter_term()
      use fks_singular_module, only: implementation => compute_collinear_counter_term
      implicit none
      call init_fks_singular_bridge()
      call implementation()
      end subroutine compute_collinear_counter_term

      logical function pdg_equal(pdg1,pdg2)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => pdg_equal
      implicit none
      integer :: pdg1(nexternal),pdg2(nexternal)
      call init_fks_singular_bridge()
      pdg_equal=implementation(pdg1,pdg2)
      end function pdg_equal

      logical function colour_con_equal(n,icol1,icol2)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => colour_con_equal
      implicit none
      integer :: n,icol1(2,nexternal),icol2(2,nexternal)
      call init_fks_singular_bridge()
      colour_con_equal=implementation(n,icol1,icol2)
      end function colour_con_equal

      logical function momenta_equal(p1,p2)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => momenta_equal
      implicit none
      double precision :: p1(0:3,nexternal),p2(0:3,nexternal)
      call init_fks_singular_bridge()
      momenta_equal=implementation(p1,p2)
      end function momenta_equal

      logical function momenta_equal_uborn(p1,p2,jfks1,ifks1,jfks2,ifks2)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => momenta_equal_uborn
      implicit none
      integer :: jfks1,ifks1,jfks2,ifks2
      double precision :: p1(0:3,nexternal),p2(0:3,nexternal)
      call init_fks_singular_bridge()
      momenta_equal_uborn=implementation(p1,p2,jfks1,ifks1,jfks2,ifks2)
      end function momenta_equal_uborn

      subroutine compute_prefactors_nbody(vegas_wgt)
      use fks_singular_module, only: implementation => compute_prefactors_nbody
      implicit none
      double precision :: vegas_wgt
      call init_fks_singular_bridge()
      call implementation(vegas_wgt)
      end subroutine compute_prefactors_nbody

      subroutine include_multichannel_enhance(imode)
      use fks_singular_module, only: implementation => include_multichannel_enhance
      implicit none
      integer :: imode
      call init_fks_singular_bridge()
      call implementation(imode)
      end subroutine include_multichannel_enhance

      subroutine compute_prefactors_n1body(vegas_wgt,jac_ev)
      use fks_singular_module, only: implementation => compute_prefactors_n1body
      implicit none
      double precision :: vegas_wgt,jac_ev
      call init_fks_singular_bridge()
      call implementation(vegas_wgt,jac_ev)
      end subroutine compute_prefactors_n1body

      subroutine add_wgt(type,orders,wgt1,wgt2,wgt3)
      use process_dimensions, only: nsplitorders
      use fks_singular_module, only: implementation => add_wgt
      implicit none
      integer :: type,orders(nsplitorders)
      double precision :: wgt1,wgt2,wgt3
      call init_fks_singular_bridge()
      call implementation(type,orders,wgt1,wgt2,wgt3)
      end subroutine add_wgt

      subroutine include_PDF_and_alphas()
      use fks_singular_module, only: implementation => include_PDF_and_alphas
      implicit none
      call init_fks_singular_bridge()
      call implementation()
      end subroutine include_PDF_and_alphas

      subroutine separate_flavour_config(ict)
      use fks_singular_module, only: implementation => separate_flavour_config
      implicit none
      integer :: ict
      call init_fks_singular_bridge()
      call implementation(ict)
      end subroutine separate_flavour_config

      subroutine set_pdg_codes(iproc,pd,iFKS,ict)
      use process_dimensions, only: maxproc
      use fks_singular_module, only: implementation => set_pdg_codes
      implicit none
      integer :: iproc,iFKS,ict
      double precision :: pd(0:maxproc)
      call init_fks_singular_bridge()
      call implementation(iproc,pd,iFKS,ict)
      end subroutine set_pdg_codes

      subroutine reweight_scale()
      use fks_singular_module, only: implementation => reweight_scale
      implicit none
      call init_fks_singular_bridge()
      call implementation()
      end subroutine reweight_scale

      subroutine reweight_pdf()
      use fks_singular_module, only: implementation => reweight_pdf
      implicit none
      call init_fks_singular_bridge()
      call implementation()
      end subroutine reweight_pdf

      subroutine fill_pineappl_weights(vegas_wgt)
      use fks_singular_module, only: implementation => fill_pineappl_weights
      implicit none
      double precision :: vegas_wgt
      call init_fks_singular_bridge()
      call implementation(vegas_wgt)
      end subroutine fill_pineappl_weights

      subroutine get_wgt_nbody(sig)
      use fks_singular_module, only: implementation => get_wgt_nbody
      implicit none
      double precision :: sig
      call init_fks_singular_bridge()
      call implementation(sig)
      end subroutine get_wgt_nbody

      subroutine get_wgt_no_nbody(sig)
      use fks_singular_module, only: implementation => get_wgt_no_nbody
      implicit none
      double precision :: sig
      call init_fks_singular_bridge()
      call implementation(sig)
      end subroutine get_wgt_no_nbody

      subroutine fill_plots()
      use fks_singular_module, only: implementation => fill_plots
      implicit none
      call init_fks_singular_bridge()
      call implementation()
      end subroutine fill_plots

      subroutine fill_mint_function(f)
      use mint_module, only: nintegrals
      use fks_singular_module, only: implementation => fill_mint_function
      implicit none
      double precision :: f(nintegrals)
      call init_fks_singular_bridge()
      call implementation(f)
      end subroutine fill_mint_function

      subroutine rotate_invar(pin,pout,cth,sth,cphi,sphi)
      use fks_singular_module, only: implementation => rotate_invar
      implicit none
      double precision :: pin(0:3),pout(0:3),cth,sth,cphi,sphi
      call init_fks_singular_bridge()
      call implementation(pin,pout,cth,sth,cphi,sphi)
      end subroutine rotate_invar

      subroutine trp_rotate_invar(pin,pout,cth,sth,cphi,sphi)
      use fks_singular_module, only: implementation => trp_rotate_invar
      implicit none
      double precision :: pin(0:3),pout(0:3),cth,sth,cphi,sphi
      call init_fks_singular_bridge()
      call implementation(pin,pout,cth,sth,cphi,sphi)
      end subroutine trp_rotate_invar

      subroutine getaziangles(p,cphi,sphi)
      use fks_singular_module, only: implementation => getaziangles
      implicit none
      double precision :: p(0:3),cphi,sphi
      call init_fks_singular_bridge()
      call implementation(p,cphi,sphi)
      end subroutine getaziangles

      subroutine phspncheck_born(ecm,xmass,xmom,pass)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => phspncheck_born
      implicit none
      double precision :: ecm,xmass(nexternal-1),xmom(0:3,nexternal-1)
      logical :: pass
      call init_fks_singular_bridge()
      call implementation(ecm,xmass,xmom,pass)
      end subroutine phspncheck_born

      subroutine phspncheck_nocms(npart,ecm,xmass,xmom,pass)
      use process_dimensions, only: nexternal,max_branch,max_particles
      use fks_singular_module, only: implementation => phspncheck_nocms
      implicit none
      integer :: npart
      double precision :: ecm,xmass(-max_branch:max_particles)
      double precision :: xmom(0:3,nexternal)
      logical :: pass
      call init_fks_singular_bridge()
      call implementation(npart,ecm,xmass,xmom,pass)
      end subroutine phspncheck_nocms

      double precision function xlen4(v)
      use fks_singular_module, only: implementation => xlen4
      implicit none
      double precision :: v(0:3)
      call init_fks_singular_bridge()
      xlen4=implementation(v)
      end function xlen4

      subroutine sreal(pp,xi_i_fks,y_ij_fks,wgt)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => sreal
      implicit none
      double precision :: pp(0:3,nexternal),xi_i_fks,y_ij_fks,wgt
      call init_fks_singular_bridge()
      call implementation(pp,xi_i_fks,y_ij_fks,wgt)
      end subroutine sreal

      subroutine sborncol_fsr(p,xi_i_fks,y_ij_fks,wgt)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => sborncol_fsr
      implicit none
      double precision :: p(0:3,nexternal),xi_i_fks,y_ij_fks,wgt
      call init_fks_singular_bridge()
      call implementation(p,xi_i_fks,y_ij_fks,wgt)
      end subroutine sborncol_fsr

      subroutine sborncol_isr(p,xi_i_fks,y_ij_fks,wgt)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => sborncol_isr
      implicit none
      double precision :: p(0:3,nexternal),xi_i_fks,y_ij_fks,wgt
      call init_fks_singular_bridge()
      call implementation(p,xi_i_fks,y_ij_fks,wgt)
      end subroutine sborncol_isr

      subroutine xkplus(PDFscheme,col1,col2,x,xkk)
      use fks_singular_module, only: implementation => xkplus
      implicit none
      integer :: PDFscheme,col1,col2
      double precision :: x,xkk(2)
      call init_fks_singular_bridge()
      call implementation(PDFscheme,col1,col2,x,xkk)
      end subroutine xkplus

      subroutine xklog(PDFscheme,col1,col2,x,xkk)
      use fks_singular_module, only: implementation => xklog
      implicit none
      integer :: PDFscheme,col1,col2
      double precision :: x,xkk(2)
      call init_fks_singular_bridge()
      call implementation(PDFscheme,col1,col2,x,xkk)
      end subroutine xklog

      subroutine xkdelta(PDFscheme,col1,col2,xkk)
      use fks_singular_module, only: implementation => xkdelta
      implicit none
      integer :: PDFscheme,col1,col2
      double precision :: xkk(2)
      call init_fks_singular_bridge()
      call implementation(PDFscheme,col1,col2,xkk)
      end subroutine xkdelta

      subroutine AP_reduced(col1,col2,t,z,ap)
      use fks_singular_module, only: implementation => AP_reduced
      implicit none
      integer :: col1,col2
      double precision :: t,z,ap(2)
      call init_fks_singular_bridge()
      call implementation(col1,col2,t,z,ap)
      end subroutine AP_reduced

      subroutine AP_reduced_prime(col1,col2,t,z,apprime)
      use fks_singular_module, only: implementation => AP_reduced_prime
      implicit none
      integer :: col1,col2
      double precision :: t,z,apprime(2)
      call init_fks_singular_bridge()
      call implementation(col1,col2,t,z,apprime)
      end subroutine AP_reduced_prime

      subroutine Qterms_reduced_timelike(col1,col2,t,z,Qterms)
      use fks_singular_module, only: implementation => Qterms_reduced_timelike
      implicit none
      integer :: col1,col2
      double precision :: t,z,Qterms(2)
      call init_fks_singular_bridge()
      call implementation(col1,col2,t,z,Qterms)
      end subroutine Qterms_reduced_timelike

      subroutine Qterms_reduced_spacelike(col1,col2,t,z,Qterms)
      use fks_singular_module, only: implementation => Qterms_reduced_spacelike
      implicit none
      integer :: col1,col2
      double precision :: t,z,Qterms(2)
      call init_fks_singular_bridge()
      call implementation(col1,col2,t,z,Qterms)
      end subroutine Qterms_reduced_spacelike

      subroutine AP_reduced_SUSY(col1,col2,t,z,ap)
      use fks_singular_module, only: implementation => AP_reduced_SUSY
      implicit none
      integer :: col1,col2
      double precision :: t,z,ap(2)
      call init_fks_singular_bridge()
      call implementation(col1,col2,t,z,ap)
      end subroutine AP_reduced_SUSY

      subroutine AP_reduced_massive(col1,col2,t,z,q2,m2,ap)
      use fks_singular_module, only: implementation => AP_reduced_massive
      implicit none
      integer :: col1,col2
      double precision :: t,z,q2,m2,ap(2)
      call init_fks_singular_bridge()
      call implementation(col1,col2,t,z,q2,m2,ap)
      end subroutine AP_reduced_massive

      subroutine sbornsoft(pp,xi_i_fks,y_ij_fks,wgt)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => sbornsoft
      implicit none
      double precision :: pp(0:3,nexternal),xi_i_fks,y_ij_fks,wgt
      call init_fks_singular_bridge()
      call implementation(pp,xi_i_fks,y_ij_fks,wgt)
      end subroutine sbornsoft

      subroutine eikonal_reduced(pp,m,n,i_fks,j_fks,xi_i_fks,y_ij_fks,eik)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => eikonal_reduced
      implicit none
      double precision :: pp(0:3,nexternal),xi_i_fks,y_ij_fks,eik
      integer :: m,n,i_fks,j_fks
      call init_fks_singular_bridge()
      call implementation(pp,m,n,i_fks,j_fks,xi_i_fks,y_ij_fks,eik)
      end subroutine eikonal_reduced

      subroutine sreal_deg(p,xi_i_fks,y_ij_fks,collrem_xi,collrem_lxi)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => sreal_deg
      implicit none
      double precision :: p(0:3,nexternal),xi_i_fks,y_ij_fks
      double precision :: collrem_xi,collrem_lxi
      call init_fks_singular_bridge()
      call implementation(p,xi_i_fks,y_ij_fks,collrem_xi,collrem_lxi)
      end subroutine sreal_deg

      subroutine set_cms_stuff(icountevts)
      use fks_singular_module, only: implementation => set_cms_stuff
      implicit none
      integer :: icountevts
      call init_fks_singular_bridge()
      call implementation(icountevts)
      end subroutine set_cms_stuff

      subroutine xmom_compare(i_fks,j_fks,jac,jac_cnt,p,p1_cnt,pass)
      use process_dimensions, only: nexternal,max_branch,max_particles
      use fks_singular_module, only: implementation => xmom_compare
      implicit none
      integer :: i_fks,j_fks
      double precision :: jac,jac_cnt(-2:2)
      double precision :: p(0:3,-max_branch:max_particles)
      double precision :: p1_cnt(0:3,nexternal,-2:2)
      logical :: pass
      call init_fks_singular_bridge()
      call implementation(i_fks,j_fks,jac,jac_cnt,p,p1_cnt,pass)
      end subroutine xmom_compare

      subroutine xmcompare(verbose,pass0,inum,iden,i_fks,j_fks,p,p1_cnt)
      use process_dimensions, only: nexternal,max_branch,max_particles
      use fks_singular_module, only: implementation => xmcompare
      implicit none
      logical :: verbose,pass0
      integer :: inum,iden,i_fks,j_fks
      double precision :: p(0:3,-max_branch:max_particles)
      double precision :: p1_cnt(0:3,nexternal,-2:2)
      call init_fks_singular_bridge()
      call implementation(verbose,pass0,inum,iden,i_fks,j_fks,p,p1_cnt)
      end subroutine xmcompare

      subroutine xprintout(iunit,xv,xlim)
      use fks_singular_module, only: implementation => xprintout
      implicit none
      integer :: iunit
      double precision :: xv,xlim
      call init_fks_singular_bridge()
      call implementation(iunit,xv,xlim)
      end subroutine xprintout

      subroutine checkres(xsecvc,xseclvc,wgt,wgtl,xp,lxp,iflag,imax,iev,
     $nexternal,i_fks,j_fks,iret)
      use fks_singular_module, only: implementation => checkres
      implicit none
      double precision :: xsecvc(15),xseclvc,wgt(15),wgtl
      double precision :: xp(15,0:3,21),lxp(0:3,21)
      integer :: iflag,imax,iev,nexternal,i_fks,j_fks,iret
      call init_fks_singular_bridge()
      call implementation(xsecvc,xseclvc,wgt,wgtl,xp,lxp,iflag,imax,iev,
     $nexternal,i_fks,j_fks,iret)
      end subroutine checkres

      subroutine checkres2(xsecvc,xseclvc,wgt,wgtl,xp,lxp,iflag,imax,iev,
     $i_fks,j_fks,iret)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => checkres2
      implicit none
      double precision :: xsecvc(15),xseclvc(15),wgt(15),wgtl(15)
      double precision :: xp(15,0:3,nexternal+1),lxp(0:3,nexternal+1)
      integer :: iflag,imax,iev,i_fks,j_fks,iret
      call init_fks_singular_bridge()
      call implementation(xsecvc,xseclvc,wgt,wgtl,xp,lxp,iflag,imax,iev,
     $i_fks,j_fks,iret)
      end subroutine checkres2

      subroutine bornsoftvirtual(p,bsv_wgt,virt_wgt,born_wgt)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => bornsoftvirtual
      implicit none
      double precision :: p(0:3,nexternal),bsv_wgt,virt_wgt,born_wgt
      call init_fks_singular_bridge()
      call implementation(p,bsv_wgt,virt_wgt,born_wgt)
      end subroutine bornsoftvirtual

      subroutine eikonal_Ireg(p,m,n,xicut_used,eikIreg)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => eikonal_Ireg
      implicit none
      double precision :: p(0:3,nexternal),xicut_used,eikIreg
      integer :: m,n
      call init_fks_singular_bridge()
      call implementation(p,m,n,xicut_used,eikIreg)
      end subroutine eikonal_Ireg

      double precision function xj1a(x,y,tHVv,tHVvl)
      use fks_singular_module, only: implementation => xj1a
      implicit none
      double precision :: x,y,tHVv,tHVvl
      call init_fks_singular_bridge()
      xj1a=implementation(x,y,tHVv,tHVvl)
      end function xj1a

      double precision function DDILOG(x)
      use fks_singular_module, only: implementation => DDILOG
      implicit none
      double precision :: x
      call init_fks_singular_bridge()
      DDILOG=implementation(x)
      end function DDILOG

      subroutine getpoles(p,xmu2,double,single,fksprefact)
      use process_dimensions, only: nexternal
      use fks_singular_module, only: implementation => getpoles
      implicit none
      double precision :: p(0:3,nexternal),xmu2,double,single
      logical :: fksprefact
      call init_fks_singular_bridge()
      call implementation(p,xmu2,double,single,fksprefact)
      end subroutine getpoles

      subroutine setfksfactor()
      use fks_singular_module, only: implementation => setfksfactor
      implicit none
      call init_fks_singular_bridge()
      call implementation()
      end subroutine setfksfactor

      subroutine set_mu_central(ic,dd,c_mu2_r,c_mu2_f)
      use fks_singular_module, only: implementation => set_mu_central
      implicit none
      integer :: ic,dd
      double precision :: c_mu2_r,c_mu2_f
      call init_fks_singular_bridge()
      call implementation(ic,dd,c_mu2_r,c_mu2_f)
      end subroutine set_mu_central

      double precision function ran2()
      use fks_singular_module, only: implementation => ran2
      implicit none
      call init_fks_singular_bridge()
      ran2=implementation()
      end function ran2

      subroutine fill_configurations_common()
      use fks_singular_module, only: implementation => fill_configurations_common
      implicit none
      call init_fks_singular_bridge()
      call implementation()
      end subroutine fill_configurations_common

      subroutine fill_configurations_born(iforest_in,sprop_in,tprid_in,
     $mapconfig_in,pmass_in,pwidth_in)
      use process_dimensions, only: nexternal,max_branch,lmaxconfigs
      use fks_singular_module, only: implementation => fill_configurations_born
      implicit none
      integer :: iforest_in(2,-max_branch:-1,lmaxconfigs)
      integer :: sprop_in(-max_branch:-1,lmaxconfigs)
      integer :: tprid_in(-max_branch:-1,lmaxconfigs)
      integer :: mapconfig_in(0:lmaxconfigs)
      double precision :: pmass_in(-nexternal:0,lmaxconfigs)
      double precision :: pwidth_in(-nexternal:0,lmaxconfigs)
      call init_fks_singular_bridge()
      call implementation(iforest_in,sprop_in,tprid_in,mapconfig_in,
     $pmass_in,pwidth_in)
      end subroutine fill_configurations_born

      subroutine fill_configurations_real(iforest_in,sprop_in,tprid_in,
     $mapconfig_in,pmass_in,pwidth_in)
      use process_dimensions, only: nexternal,max_branch,lmaxconfigs
      use fks_singular_module, only: implementation => fill_configurations_real
      implicit none
      integer :: iforest_in(2,-max_branch:-1,lmaxconfigs)
      integer :: sprop_in(-max_branch:-1,lmaxconfigs)
      integer :: tprid_in(-max_branch:-1,lmaxconfigs)
      integer :: mapconfig_in(0:lmaxconfigs)
      double precision :: pmass_in(-nexternal:0,lmaxconfigs)
      double precision :: pwidth_in(-nexternal:0,lmaxconfigs)
      call init_fks_singular_bridge()
      call implementation(iforest_in,sprop_in,tprid_in,mapconfig_in,
     $pmass_in,pwidth_in)
      end subroutine fill_configurations_real
