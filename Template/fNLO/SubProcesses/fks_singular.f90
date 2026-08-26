module fks_singular_module
  use process_dimensions, only: nexternal, nincoming, max_particles, &
                                max_branch, lmaxconfigs, maxproc, ngraphs, ncolor, maxflow, fks_configs, &
                                nsplitorders, qcd_pos, qed_pos, amp_split_size, amp_split_size_born, &
                                validate_process_dimensions
  use run_state
  use timing_state
  use mint_module, only: maxchannels
  use fks_metadata, only: validate_fks_metadata, fks_i_d, fks_j_d, &
                          extra_cnt_d, isplitorder_born_d, isplitorder_cnt_d, split_type_d
  use setscales_module, only: set_alphas, set_ren_scale, set_fac_scale
  use split_orders, only: get_orders_tag, amp_split_pos_to_orders, &
                          lo_qcd_to_amp_pos, nlo_qcd_to_amp_pos
  use kin_functions_module, only: dot => dot_impl, rho => rho_impl
  use fks_sij_module, only: initialize_fks_sij_module, &
                            set_fks_sij_partition_state, fks_sij_impl
  use FKSParams, only: use_poly_virtual
  use chooser_functions_module, only: get_mother_colour_impl, set_pdg_impl
  use madfks_plot_module, only: initplot_impl, outfun_impl
  use fnlo_process_common, only: nfksprocess, i_fks, j_fks, &
                                 ybst_til_tolab, ybst_til_tocm, &
                                 sqrtshat, shat, xi_i_fks_ev, &
                                 y_ij_fks_ev, p_i_fks_ev, p_i_fks_cnt, &
                                 xi_i_fks_cnt, xi_i_hat_ev, &
                                 xiimax_ev, xiimax_cnt, xinorm_ev, &
                                 xinorm_cnt, f_b, f_nb, f_r, f_s, f_c, &
                                 f_dc, f_sc, f_dsc, f_pdfsch_d, &
                                 f_pdfsch_p, f_pdfsch_l, xiscut_used, &
                                 xibsvcut_used, delta_used, xicut_used, &
                                 fkssymmetryfactor, &
                                 fkssymmetryfactorborn, &
                                 fkssymmetryfactordeg, ngluons, &
                                 rndec, this_config, &
                                 diagramsymmetryfactor, &
                                 calculatedborn => calculated_born, &
                                 use_evpr, nocntevents, wgt_me_born, &
                                 wgt_me_real, &
                                 fold, ifold_counter, fixed_order, &
                                 orders_tag_plot, amp_pos_plot, &
                                 virtual_over_born, softtest, colltest, &
                                 need_color_links, xij_aor, &
                                 i_type, j_type, m_type, &
                                 iextra_cnt, isplitorder_born, &
                                 isplitorder_cnt, iden_comp, sqrtshat_ev, &
                                 shat_ev, sqrtshat_cnt, shat_cnt, &
                                 ycm_ev, ycm_cnt, xbjrk_ev, &
                                 xbjrk_cnt, i_momcmp_count, xratmax, c, &
                                 gamma, gammap, beta0, abrv, &
                                 multi_channel, nbody
  implicit none
  private

  double precision, parameter :: fks_a = 1.5d0, fks_b = 1.5d0
  double precision, parameter :: a_h_damp = 1d0, one_h_damp = 1d-2
  logical, parameter :: useenergy = .true., usebeta = .true.
  double precision, parameter :: deltao = 1d0, deltai = 1d0, xicut = 0.5d0
  double precision, parameter :: deltas = 1d0, xiscut = 0.5d0, xibsvcut = 1d0
  double precision, parameter :: deltaminy = 0.95d0, skewy = 10d0, alphay = 2d0

  double precision :: nf = 0d0

  double precision, pointer :: g => null(), qes2 => null()
  double precision, pointer :: amp_split(:) => null()
  complex(kind=kind(0d0)), pointer :: amp_split_cnt(:, :, :) => null()

  double precision, pointer :: p_born(:, :) => null()
  double precision, pointer :: p_born_coll(:, :) => null()
  double precision, pointer :: p_born_norad(:, :) => null()
  double precision, pointer :: p_ev(:, :) => null()
  double precision, pointer :: p1_cnt(:, :, :) => null()
  double precision, pointer :: wgt_cnt(:) => null(), pswgt_cnt(:) => null()
  double precision, pointer :: jac_cnt(:) => null()
  integer, pointer :: idup(:, :) => null(), mothup(:, :, :) => null()
  integer, pointer :: icolup(:, :, :) => null(), niprocs => null()
  integer, pointer :: idup_d(:, :, :) => null()
  integer, pointer :: fks_j_from_i(:, :) => null()
  integer, pointer :: particle_type(:) => null(), pdg_type(:) => null()
  logical, pointer :: split_type(:) => null()
  logical, pointer :: is_aorg(:) => null()
  complex(kind=kind(0d0)), pointer :: ans_cnt(:, :) => null()

  double precision, pointer :: amp_split_virt(:) => null()
  double precision, pointer :: amp_split_born_for_virt(:) => null()
  double precision, pointer :: amp_split_avv(:) => null()
  double precision, pointer :: amp_split_wgtnstmp(:) => null()
  double precision, pointer :: amp_split_wgtwnstmpmuf(:) => null()
  double precision, pointer :: amp_split_wgtwnstmpmur(:) => null()
  double precision, pointer :: amp_split_wgtdegrem_xi(:) => null()
  double precision, pointer :: amp_split_wgtdegrem_lxi(:) => null()
  double precision, pointer :: amp_split_wgtdegrem_muf(:) => null()
  double precision, pointer :: amp_split_wgtpsch_p(:) => null()
  double precision, pointer :: amp_split_wgtpsch_l(:) => null()
  double precision, pointer :: amp_split_wgtpsch_d(:) => null()
  double precision, pointer :: amp_split_soft(:) => null()
  double precision, pointer :: amp_split_finite_ml(:) => null()
  double precision, pointer :: amp_split_poles_fks(:, :) => null()
  logical, pointer :: split_type_used(:) => null()

  double precision, pointer :: amp2(:) => null(), jamp2(:) => null()

  double precision, pointer :: config_mass(:, :, :) => null()
  double precision, pointer :: config_width(:, :, :) => null()
  integer, pointer :: config_forest(:, :, :, :) => null()
  integer, pointer :: config_sprop(:, :, :) => null()
  integer, pointer :: config_tprid(:, :, :) => null()
  integer, pointer :: config_map(:, :) => null()
  integer, pointer :: real_forest(:, :, :) => null()
  integer, pointer :: real_sprop(:, :) => null(), real_tprid(:, :) => null()
  integer, pointer :: real_map(:) => null(), real_prow(:, :) => null()
  double precision, pointer :: real_mass(:, :) => null()
  double precision, pointer :: real_width(:, :) => null()

  double precision, pointer :: subproc_pd(:) => null()
  integer, pointer :: subproc_iproc => null()
  integer, pointer :: flavour_map(:) => null(), iproc_save(:) => null()
  integer, pointer :: eto(:, :) => null(), etoi(:, :) => null()
  integer, pointer :: maxproc_found => null()

  integer, pointer :: appl_amp_split_size => null()
  integer, pointer :: appl_qcdpower(:) => null(), appl_qedpower(:) => null()
  integer, pointer :: appl_nproc(:) => null()
  double precision, pointer :: appl_x1(:) => null(), appl_x2(:) => null()
  double precision, pointer :: appl_muf2(:) => null(), appl_mur2(:) => null()
  double precision, pointer :: appl_qes2(:) => null()
  double precision, pointer :: appl_w0(:, :) => null(), appl_wr(:, :) => null()
  double precision, pointer :: appl_wf(:, :) => null(), appl_wb(:, :) => null()
  integer, pointer :: appl_flavmap(:) => null()
  double precision, pointer :: appl_event_weight => null()
  double precision, pointer :: appl_vegaswgt => null()

  double precision, allocatable, save :: external_masses(:)
  integer, allocatable, save :: born_forest(:, :, :), born_sprop(:, :)
  integer, allocatable, save :: born_tprid(:, :), born_map(:)
  double precision, allocatable, save :: born_mass(:, :), born_width(:, :)
  integer, save :: born_max_branch_used = 0, born_lmaxconfigs_used = 0

  double precision, allocatable, save :: amp_split_virt_save(:)
  logical, allocatable, save :: firsttime_nfksprocess(:)
  double precision, allocatable, save :: diagramsymmetryfactor_save(:)
  integer, allocatable, save :: fac_i_fks(:), fac_j_fks(:)
  integer, allocatable, save :: i_type_fks(:), j_type_fks(:), m_type_fks(:)
  integer, allocatable, save :: ngluons_fks(:)
  integer, allocatable, save :: iden_real_fks(:), iden_born_fks(:)
  logical, save :: setfks_firsttime = .true.
  logical, save :: fks_singular_state_initialized = .false.

  public :: compute_born
  public :: compute_nbody_noborn, compute_real_emission
  public :: compute_soft_counter_term, compute_collinear_counter_term
  public :: compute_soft_collinear_ct_impl
  public :: compute_prefactors_nbody, include_multichannel_enhance
  public :: compute_prefactors_n1body, include_pdf_and_alphas
  public :: reweight_scale
  public :: reweight_pdf, fill_pineappl_weights, get_wgt_nbody
  public :: get_wgt_no_nbody, fill_plots, fill_mint_function
  public :: rotate_invar, phspncheck_born, phspncheck_nocms
  public :: sreal, set_cms_stuff, xmom_compare, xprintout, checkres2
  public :: getpoles, setfksfactor, ran2, fill_configurations_common
  public :: initialize_fks_model_state, initialize_fks_phase_state
  public :: initialize_fks_amplitude_state, initialize_fks_config_state
  public :: initialize_fks_pineappl_state, initialize_fks_generated_state
  public :: validate_fks_singular_state

contains

  double precision function evaluate_fks_sij(p, ii_fks, jj_fks, &
                                             xi_i_fks, y_ij_fks)
    implicit none
    double precision, intent(in) :: p(0:3, nexternal)
    double precision, intent(in) :: xi_i_fks, y_ij_fks
    integer, intent(in) :: ii_fks, jj_fks

    call require_fks_singular_state()
    call initialize_fks_sij_module(nexternal, nincoming, fks_a, fks_b, &
                                   a_h_damp, one_h_damp, useenergy, usebeta)
    call set_fks_sij_partition_state(fks_j_from_i, particle_type, is_aorg, &
                                     i_fks, j_fks, ybst_til_tocm, sqrtshat, shat, &
                                     p_i_fks_cnt, external_masses)
    evaluate_fks_sij = fks_sij_impl(p, ii_fks, jj_fks, xi_i_fks, y_ij_fks)
  end function evaluate_fks_sij

  subroutine initialize_fks_model_state(g_in, qes2_in, nf_in, &
                                        external_masses_in)
    implicit none
    double precision, target, intent(inout) :: g_in, qes2_in
    double precision, intent(in) :: nf_in, external_masses_in(:)
    call validate_process_dimensions()
    if (size(external_masses_in) /= nexternal) then
      call fail_fks_singular_state('invalid generated model-data shape')
    end if
    g => g_in
    qes2 => qes2_in
    nf = nf_in
    if (.not. allocated(external_masses)) allocate (external_masses(nexternal))
    external_masses = external_masses_in
    call allocate_fks_singular_caches()
  end subroutine initialize_fks_model_state

  subroutine initialize_fks_phase_state(p_born_in, p_born_coll_in, &
                                        p_born_norad_in, p_ev_in, p1_cnt_in, wgt_cnt_in, pswgt_cnt_in, &
                                        jac_cnt_in, idup_in, mothup_in, icolup_in, niprocs_in, &
                                        is_aorg_in, amp2_in, jamp2_in, subproc_pd_in, subproc_iproc_in, flavour_map_in, &
                                        iproc_save_in, eto_in, etoi_in, maxproc_found_in)
    implicit none
    double precision, target, intent(inout) :: p_born_in(0:, 1:)
    double precision, target, intent(inout) :: p_born_coll_in(0:, 1:)
    double precision, target, intent(inout) :: p_born_norad_in(0:, 1:)
    double precision, target, intent(inout) :: p_ev_in(0:, 1:)
    double precision, target, intent(inout) :: p1_cnt_in(0:, 1:, 0:)
    double precision, target, intent(inout) :: wgt_cnt_in(0:)
    double precision, target, intent(inout) :: pswgt_cnt_in(0:)
    double precision, target, intent(inout) :: jac_cnt_in(0:)
    integer, target, intent(inout) :: idup_in(1:, 1:), mothup_in(1:, 1:, 1:)
    integer, target, intent(inout) :: icolup_in(1:, 1:, 1:), niprocs_in
    logical, target, intent(inout) :: is_aorg_in(1:)
    double precision, target, intent(inout) :: amp2_in(1:), jamp2_in(0:)
    double precision, target, intent(inout) :: subproc_pd_in(0:)
    integer, target, intent(inout) :: subproc_iproc_in
    integer, target, intent(inout) :: flavour_map_in(1:), iproc_save_in(1:)
    integer, target, intent(inout) :: eto_in(1:, 1:), etoi_in(1:, 1:)
    integer, target, intent(inout) :: maxproc_found_in
    call validate_process_dimensions()
    if (size(p1_cnt_in, 1) /= 4 .or. &
        size(p1_cnt_in, 2) /= nexternal .or. &
        size(p1_cnt_in, 3) /= 3 .or. &
        size(wgt_cnt_in) /= 3 .or. &
        size(pswgt_cnt_in) /= 3 .or. &
        size(jac_cnt_in) /= 3) then
      call fail_fks_singular_state('invalid counterevent storage shape')
    end if
    p_born => p_born_in
    p_born_coll => p_born_coll_in
    p_born_norad => p_born_norad_in
    p_ev => p_ev_in
    p1_cnt => p1_cnt_in
    wgt_cnt => wgt_cnt_in
    pswgt_cnt => pswgt_cnt_in
    jac_cnt => jac_cnt_in
    idup => idup_in
    mothup => mothup_in
    icolup => icolup_in
    niprocs => niprocs_in
    is_aorg => is_aorg_in
    amp2 => amp2_in
    jamp2 => jamp2_in
    subproc_pd => subproc_pd_in
    subproc_iproc => subproc_iproc_in
    flavour_map => flavour_map_in
    iproc_save => iproc_save_in
    eto => eto_in
    etoi => etoi_in
    maxproc_found => maxproc_found_in
  end subroutine initialize_fks_phase_state

  subroutine initialize_fks_amplitude_state(amp_split_in, amp_split_cnt_in, &
                                            amp_virt_in, amp_born_virt_in, amp_avv_in, &
                                            amp_bsv_in, amp_bsv_muf_in, amp_bsv_mur_in, amp_deg_xi_in, amp_deg_lxi_in, &
                                            amp_deg_muf_in, amp_dis_p_in, amp_dis_l_in, amp_dis_d_in, amp_soft_in, &
                                            amp_finite_in, amp_poles_in, fks_j_from_i_in, particle_type_in, &
                                            pdg_type_in, split_type_in, ans_cnt_in, split_type_used_in, idup_d_in)
    implicit none
    double precision, target, intent(inout) :: amp_split_in(1:)
    complex(kind=kind(0d0)), target, intent(inout) :: amp_split_cnt_in(1:, 1:, 1:)
    double precision, target, intent(inout) :: amp_virt_in(1:), amp_born_virt_in(1:), amp_avv_in(1:)
    double precision, target, intent(inout) :: amp_bsv_in(1:), amp_bsv_muf_in(1:), amp_bsv_mur_in(1:)
    double precision, target, intent(inout) :: amp_deg_xi_in(1:), amp_deg_lxi_in(1:), amp_deg_muf_in(1:)
    double precision, target, intent(inout) :: amp_dis_p_in(1:), amp_dis_l_in(1:), amp_dis_d_in(1:)
    double precision, target, intent(inout) :: amp_soft_in(1:), amp_finite_in(1:), amp_poles_in(1:, 1:)
    integer, target, intent(inout) :: fks_j_from_i_in(1:, 0:), particle_type_in(1:), pdg_type_in(1:)
    logical, target, intent(inout) :: split_type_in(1:), split_type_used_in(1:)
    complex(kind=kind(0d0)), target, intent(inout) :: ans_cnt_in(1:, 1:)
    integer, target, intent(inout) :: idup_d_in(1:, 1:, 1:)
    amp_split => amp_split_in
    amp_split_cnt => amp_split_cnt_in
    amp_split_virt => amp_virt_in
    amp_split_born_for_virt => amp_born_virt_in
    amp_split_avv => amp_avv_in
    amp_split_wgtnstmp => amp_bsv_in
    amp_split_wgtwnstmpmuf => amp_bsv_muf_in
    amp_split_wgtwnstmpmur => amp_bsv_mur_in
    amp_split_wgtdegrem_xi => amp_deg_xi_in
    amp_split_wgtdegrem_lxi => amp_deg_lxi_in
    amp_split_wgtdegrem_muf => amp_deg_muf_in
    amp_split_wgtpsch_p => amp_dis_p_in
    amp_split_wgtpsch_l => amp_dis_l_in
    amp_split_wgtpsch_d => amp_dis_d_in
    amp_split_soft => amp_soft_in
    amp_split_finite_ml => amp_finite_in
    amp_split_poles_fks => amp_poles_in
    fks_j_from_i => fks_j_from_i_in
    particle_type => particle_type_in
    pdg_type => pdg_type_in
    split_type => split_type_in
    ans_cnt => ans_cnt_in
    split_type_used => split_type_used_in
    idup_d => idup_d_in
  end subroutine initialize_fks_amplitude_state

  subroutine initialize_fks_config_state(config_mass_in, config_width_in, &
                                         config_forest_in, config_sprop_in, config_tprid_in, config_map_in, &
                                         real_forest_in, real_sprop_in, real_tprid_in, real_map_in, real_mass_in, &
                                         real_width_in, real_prow_in)
    implicit none
    double precision, target, intent(inout) :: config_mass_in(-nexternal:, 1:, 0:), config_width_in(-nexternal:, 1:, 0:)
    integer, target, intent(inout) :: config_forest_in(1:, -max_branch:, 1:, 0:)
    integer, target, intent(inout) :: config_sprop_in(-max_branch:, 1:, 0:), config_tprid_in(-max_branch:, 1:, 0:)
    integer, target, intent(inout) :: config_map_in(0:, 0:)
    integer, target, intent(inout) :: real_forest_in(1:, -max_branch:, 1:)
    integer, target, intent(inout) :: real_sprop_in(-max_branch:, 1:), real_tprid_in(-max_branch:, 1:), real_map_in(0:)
    double precision, target, intent(inout) :: real_mass_in(-max_branch:, 1:), real_width_in(-max_branch:, 1:)
    integer, target, intent(inout) :: real_prow_in(-max_branch:, 1:)
    config_mass => config_mass_in
    config_width => config_width_in
    config_forest => config_forest_in
    config_sprop => config_sprop_in
    config_tprid => config_tprid_in
    config_map => config_map_in
    real_forest => real_forest_in
    real_sprop => real_sprop_in
    real_tprid => real_tprid_in
    real_map => real_map_in
    real_mass => real_mass_in
    real_width => real_width_in
    real_prow => real_prow_in
  end subroutine initialize_fks_config_state

  subroutine initialize_fks_pineappl_state(appl_amp_size_in, appl_qcd_in, &
                                           appl_qed_in, appl_nproc_in, appl_x1_in, appl_x2_in, appl_muf2_in, &
                                           appl_mur2_in, appl_qes2_in, appl_w0_in, appl_wr_in, appl_wf_in, appl_wb_in, &
                                           appl_flavmap_in, appl_event_weight_in, appl_vegaswgt_in)
    implicit none
    integer, target, intent(inout) :: appl_amp_size_in, appl_qcd_in(1:), appl_qed_in(1:), appl_nproc_in(1:)
    double precision, target, intent(inout) :: appl_x1_in(1:), appl_x2_in(1:), appl_muf2_in(1:), appl_mur2_in(1:), appl_qes2_in(1:)
    double precision, target, intent(inout) :: appl_w0_in(1:, 1:), appl_wr_in(1:, 1:), appl_wf_in(1:, 1:), appl_wb_in(1:, 1:)
    integer, target, intent(inout) :: appl_flavmap_in(1:)
    double precision, target, intent(inout) :: appl_event_weight_in, appl_vegaswgt_in
    appl_amp_split_size => appl_amp_size_in
    appl_qcdpower => appl_qcd_in
    appl_qedpower => appl_qed_in
    appl_nproc => appl_nproc_in
    appl_x1 => appl_x1_in
    appl_x2 => appl_x2_in
    appl_muf2 => appl_muf2_in
    appl_mur2 => appl_mur2_in
    appl_qes2 => appl_qes2_in
    appl_w0 => appl_w0_in
    appl_wr => appl_wr_in
    appl_wf => appl_wf_in
    appl_wb => appl_wb_in
    appl_flavmap => appl_flavmap_in
    appl_event_weight => appl_event_weight_in
    appl_vegaswgt => appl_vegaswgt_in
  end subroutine initialize_fks_pineappl_state

  subroutine initialize_fks_generated_state(max_branch_used, lmax_used, &
                                            forest_in, sprop_in, tprid_in, map_in, mass_in, width_in)
    implicit none
    integer, intent(in) :: max_branch_used, lmax_used
    integer, intent(in) :: forest_in(1:, -max_branch_used:, 1:)
    integer, intent(in) :: sprop_in(-max_branch_used:, 1:), tprid_in(-max_branch_used:, 1:), map_in(0:)
    double precision, intent(in) :: mass_in(-nexternal:, 1:), width_in(-nexternal:, 1:)
    if (max_branch_used < 1 .or. max_branch_used > max_branch .or. &
        lmax_used < 1 .or. lmax_used > lmaxconfigs) then
      call fail_fks_singular_state('invalid generated Born configuration bounds')
    end if
    if (.not. allocated(born_forest)) then
      allocate (born_forest(2, -max_branch_used:-1, lmax_used))
      allocate (born_sprop(-max_branch_used:-1, lmax_used))
      allocate (born_tprid(-max_branch_used:-1, lmax_used))
      allocate (born_map(0:lmax_used))
      allocate (born_mass(-nexternal:0, lmax_used))
      allocate (born_width(-nexternal:0, lmax_used))
    else if (born_max_branch_used /= max_branch_used .or. &
             born_lmaxconfigs_used /= lmax_used) then
      call fail_fks_singular_state('generated Born configuration shape changed')
    end if
    born_max_branch_used = max_branch_used
    born_lmaxconfigs_used = lmax_used
    born_forest = forest_in
    born_sprop = sprop_in
    born_tprid = tprid_in
    born_map = map_in
    born_mass = mass_in
    born_width = width_in
  end subroutine initialize_fks_generated_state

  subroutine allocate_fks_singular_caches()
    implicit none
    if (fks_singular_state_initialized) return
    if (.not. allocated(amp_split_virt_save)) allocate (amp_split_virt_save(amp_split_size))
    if (.not. allocated(firsttime_nfksprocess)) allocate (firsttime_nfksprocess(fks_configs))
    if (.not. allocated(diagramsymmetryfactor_save)) allocate (diagramsymmetryfactor_save(maxchannels))
    if (.not. allocated(fac_i_fks)) then
      allocate (fac_i_fks(fks_configs), fac_j_fks(fks_configs))
      allocate (i_type_fks(fks_configs), j_type_fks(fks_configs), m_type_fks(fks_configs))
      allocate (ngluons_fks(fks_configs))
      allocate (iden_real_fks(fks_configs), iden_born_fks(fks_configs))
    end if
    amp_split_virt_save = 0d0
    firsttime_nfksprocess = .true.
    diagramsymmetryfactor_save = 0d0
    fac_i_fks = 0
    fac_j_fks = 0
    i_type_fks = 0
    j_type_fks = 0
    m_type_fks = 0
    ngluons_fks = 0
    iden_real_fks = 0
    iden_born_fks = 0
    fks_singular_state_initialized = .true.
  end subroutine allocate_fks_singular_caches

  subroutine validate_fks_singular_state()
    implicit none
    call validate_process_dimensions()
    call validate_fks_metadata()
    if (.not. fks_singular_state_initialized) call fail_fks_singular_state('state is not initialized')
    if (.not. associated(g) .or. .not. associated(amp_split) .or. &
        .not. associated(amp_split_cnt)) then
      call fail_fks_singular_state('model/amplitude state is not bound')
    end if
    if (.not. associated(p_born) .or. .not. associated(p1_cnt) .or. &
        .not. associated(config_mass) .or. .not. associated(fks_j_from_i) .or. &
        .not. associated(is_aorg)) then
      call fail_fks_singular_state('phase-space state is not bound')
    end if
    if (.not. associated(appl_w0) .or. .not. associated(appl_nproc)) then
      call fail_fks_singular_state('PineAPPL state is not bound')
    end if
    if (.not. allocated(born_forest) .or. .not. allocated(external_masses)) then
      call fail_fks_singular_state('generated data are not initialized')
    end if
  end subroutine validate_fks_singular_state

  subroutine require_fks_singular_state()
    implicit none
    call validate_fks_singular_state()
  end subroutine require_fks_singular_state

  subroutine fail_fks_singular_state(message)
    implicit none
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in fks_singular_module: '//trim(message)
    stop 1
  end subroutine fail_fks_singular_state

  subroutine compute_born
! This subroutine computes the Born matrix elements and adds its value
! to the list of weights using the add_wgt subroutine
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp

    double precision wgt_c
    double precision wgt1
    call cpu_time(tBefore)
    if (f_b .eq. 0d0) return
    if (xi_i_hat_ev*xiimax_cnt(0) .gt. xiBSVcut_used) return
    call sborn(p_born, wgt_c)
    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      wgt1 = amp_split(iamp)*f_b/g**(qcd_power)
      call add_wgt(2, orders, wgt1, 0d0, 0d0)
    end do

    call cpu_time(tAfter)
    tBorn = tBorn + (tAfter - tBefore)
    return
  end subroutine compute_born

  subroutine compute_nbody_noborn
! This subroutine computes the soft-virtual matrix elements and adds its
! value to the list of weights using the add_wgt subroutine
    use extra_weights
    use mint_module
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp

    double precision wgt1, wgt2, wgt3, bsv_wgt, virt_wgt, born_wgt, g22, wgt4
    call cpu_time(tBefore)
    if (f_nb .eq. 0d0) return
    if (xi_i_hat_ev*xiimax_cnt(0) .gt. xiBSVcut_used) return
    call bornsoftvirtual(p1_cnt(:, :, 0), bsv_wgt, virt_wgt, born_wgt)
    do iamp = 1, amp_split_size
      if (amp_split_wgtnstmp(iamp) .eq. 0d0 .and. &
          amp_split_wgtwnstmpmur(iamp) .eq. 0d0 .and. &
          amp_split_wgtwnstmpmuf(iamp) .eq. 0d0 .and. &
          amp_split_avv(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = amp_split_wgtnstmp(iamp)*f_nb/g22
      wgt2 = amp_split_wgtwnstmpmur(iamp)*f_nb/g22
      wgt3 = amp_split_wgtwnstmpmuf(iamp)*f_nb/g22
      wgt4 = amp_split_avv(iamp)*f_nb/g22
      call add_wgt(3, orders, wgt1, wgt2, wgt3)
      call add_wgt(15, orders, wgt4, 0d0, 0d0)
    end do
! Special for the soft-virtual needed for the virt-tricks. The
! *_wgt_mint variable should be directly passed to the mint-integrator
! and not be part of the plots nor computation of the cross section.
    virt_wgt_mint(0) = virt_wgt_mint(0) + virt_wgt*f_nb
    born_wgt_mint(0) = born_wgt_mint(0) + born_wgt*f_b
    do iamp = 1, amp_split_size
      if (amp_split_virt(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      wgt1 = amp_split_virt(iamp)*f_nb
      virt_wgt_mint(iamp) = virt_wgt_mint(iamp) + wgt1
      born_wgt_mint(iamp) = born_wgt_mint(iamp) + amp_split_born_for_virt(iamp)*f_nb
      wgt1 = wgt1/g**(QCD_power)
      call add_wgt(14, orders, wgt1, 0d0, 0d0)
    end do

    call cpu_time(tAfter)
    tIS = tIS + (tAfter - tBefore)
    return
  end subroutine compute_nbody_noborn

  subroutine compute_real_emission(p)
! This subroutine computes the real-emission matrix elements and adds
! its value to the list of weights using the add_wgt subroutine
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp
    double precision s_ev, p(0:3, nexternal), wgt1, fx_ev
    call cpu_time(tBefore)
    if (f_r .eq. 0d0) return
    s_ev = evaluate_fks_sij(p, i_fks, j_fks, xi_i_fks_ev, y_ij_fks_ev)
    if (s_ev .le. 0.d0) return
    call sreal(p, xi_i_fks_ev, y_ij_fks_ev, fx_ev)
    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      wgt1 = amp_split(iamp)*s_ev*f_r/g**(qcd_power)
      call add_wgt(1, orders, wgt1, 0d0, 0d0)
    end do
    call cpu_time(tAfter)
    tReal = tReal + (tAfter - tBefore)
    return
  end subroutine compute_real_emission

  subroutine compute_soft_counter_term
! This subroutine computes the soft counter term and adds its value to
! the list of weights using the add_wgt subroutine
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp
    double precision wgt1, s_s, fx_s, zero, g22
    parameter(zero=0d0)
    call cpu_time(tBefore)
    if (f_s .eq. 0d0) return
    if (xi_i_hat_ev*xiimax_cnt(0) .gt. xiScut_used) return
    s_s = evaluate_fks_sij(p1_cnt(:, :, 0), i_fks, j_fks, zero, y_ij_fks_ev)
    if (s_s .le. 0d0) return
    call sreal(p1_cnt(:, :, 0), 0d0, y_ij_fks_ev, fx_s)

    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = 0d0
      if (xi_i_fks_ev .le. xiScut_used) then
        wgt1 = -amp_split(iamp)*s_s*f_s/g22
      end if
      if (wgt1 .ne. 0d0) call add_wgt(4, orders, wgt1, 0d0, 0d0)
    end do

    call cpu_time(tAfter)
    tCount = tCount + (tAfter - tBefore)
    return
  end subroutine compute_soft_counter_term

  subroutine compute_collinear_counter_term
! This subroutine computes the collinear counter term and adds its value
! to the list of weights using the add_wgt subroutine
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp
! amp_split for the PDF scheme
    double precision one, s_c, fx_c, deg_xi_c, deg_lxi_c, wgt1, wgt3, g22
    parameter(one=1d0)
    double precision pmass(nexternal)
    call cpu_time(tBefore)
    pmass = external_masses
    if (f_c .eq. 0d0 .and. f_dc .eq. 0d0) return
    if (y_ij_fks_ev .le. 1d0 - deltaS .or. pmass(j_fks) .ne. 0.d0) return
    s_c = evaluate_fks_sij(p1_cnt(:, :, 1), i_fks, j_fks, xi_i_fks_cnt(1), one)
    if (s_c .le. 0d0) return
! sreal_deg should be called **BEFORE** sreal
! in order not to overwrtie the amp_split array
    call sreal_deg(p1_cnt(:, :, 1), xi_i_fks_cnt(1), deg_xi_c, deg_lxi_c)
    call sreal(p1_cnt(:, :, 1), xi_i_fks_cnt(1), one, fx_c)

    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0 .and. &
          amp_split_wgtdegrem_xi(iamp) .eq. 0d0 .and. &
          amp_split_wgtdegrem_lxi(iamp) .eq. 0d0 .and. &
          amp_split_wgtpsch_p(iamp) .eq. 0d0 .and. &
          amp_split_wgtpsch_l(iamp) .eq. 0d0 .and. &
          amp_split_wgtpsch_d(iamp) .eq. 0d0) cycle

      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = -amp_split(iamp)*s_c*f_c/g22
      wgt1 = wgt1 + (amp_split_wgtdegrem_xi(iamp) + &
                     amp_split_wgtpsch_p(iamp) + &
                     (amp_split_wgtdegrem_lxi(iamp) + &
                      amp_split_wgtpsch_l(iamp))*log(xi_i_fks_cnt(1)))*f_dc/g22
      wgt3 = amp_split_wgtdegrem_muF(iamp)*f_dc/g22
      if (wgt1 .ne. 0d0 .or. wgt3 .ne. 0d0) call add_wgt(5, orders, wgt1, 0d0, wgt3)
    end do

    call cpu_time(tAfter)
    tCount = tCount + (tAfter - tBefore)
    return
  end subroutine compute_collinear_counter_term

  subroutine compute_soft_collinear_ct_impl
! This subroutine computes the soft-collinear counter term and adds its
! value to the list of weights using the add_wgt subroutine
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer iamp
! amp_split for the PDF scheme
    double precision zero, one, s_sc, fx_sc, wgt1, wgt3, deg_xi_sc, deg_lxi_sc, g22
    parameter(zero=0d0, one=1d0)
! PDF scheme prefactors
    double precision pmass(nexternal)
    pmass = external_masses
    call cpu_time(tBefore)
    if (f_sc .eq. 0d0 .and. f_dsc(1) .eq. 0d0 .and. f_dsc(2) .eq. 0d0 .and. f_dsc(3) .eq. 0d0 .and. f_dsc(4) .eq. 0d0) return
    if (xi_i_hat_ev*xiimax_cnt(1) .ge. xiScut_used .or. y_ij_fks_ev .le. 1d0 - deltaS .or. pmass(j_fks) .ne. 0.d0) return
    s_sc = evaluate_fks_sij(p1_cnt(:, :, 2), i_fks, j_fks, zero, one)
    if (s_sc .le. 0d0) return
! sreal_deg should be called **BEFORE** sreal
! in order not to overwrtie the amp_split array
    call sreal_deg(p1_cnt(:, :, 2), zero, deg_xi_sc, deg_lxi_sc)
    call sreal(p1_cnt(:, :, 2), zero, one, fx_sc)

    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0 .and. &
          amp_split_wgtdegrem_xi(iamp) .eq. 0d0 .and. &
          amp_split_wgtdegrem_lxi(iamp) .eq. 0d0 .and. &
          amp_split_wgtpsch_p(iamp) .eq. 0d0 .and. &
          amp_split_wgtpsch_l(iamp) .eq. 0d0 .and. &
          amp_split_wgtpsch_d(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = 0d0
      wgt3 = 0d0
      if (xi_i_fks_cnt(1) .lt. xiScut_used) then
        wgt1 = amp_split(iamp)*s_sc*f_sc/g22
        wgt1 = wgt1 + ( &
               -(amp_split_wgtdegrem_xi(iamp) + amp_split_wgtpsch_p(iamp) + &
                 (amp_split_wgtdegrem_lxi(iamp) + amp_split_wgtpsch_l(iamp))* &
                 log(xi_i_fks_cnt(1)))*f_dsc(1) &
               - (amp_split_wgtdegrem_xi(iamp)*f_dsc(2) + &
                  amp_split_wgtdegrem_lxi(iamp)*f_dsc(3)) &
               + amp_split_wgtpsch_d(iamp)*f_pdfsch_d &
               + amp_split_wgtpsch_p(iamp)*f_pdfsch_p &
               + amp_split_wgtpsch_l(iamp)*f_pdfsch_l)/g22
        wgt3 = -amp_split_wgtdegrem_muF(iamp)*f_dsc(4)/g22
      end if
      if (wgt1 .ne. 0d0 .or. wgt3 .ne. 0d0) call add_wgt(6, orders, wgt1, 0d0, wgt3)
    end do

    call cpu_time(tAfter)
    tCount = tCount + (tAfter - tBefore)
    return
  end subroutine compute_soft_collinear_ct_impl

  logical function pdg_equal(pdg1, pdg2)
! Returns .true. if the lists of PDG codes --'pdg1' and 'pdg2'-- are
! equal.
    implicit none
    integer i, pdg1(nexternal), pdg2(nexternal)
    pdg_equal = .true.
    do i = 1, nexternal
      if (pdg1(i) .ne. pdg2(i)) then
        pdg_equal = .false.
        return
      end if
    end do
  end function pdg_equal

  logical function momenta_equal(p1, p2)
! Returns .true. only if the momenta p1 and p2 are equal. To save time,
! it only checks the 0th and 3rd components (energy and z-direction).
    implicit none
    integer i, j
    double precision p1(0:3, nexternal), p2(0:3, nexternal), vtiny
    parameter(vtiny=1d-8)
    momenta_equal = .true.
    do i = 1, nexternal
      do j = 0, 3, 3
        if (p1(j, i) .eq. 0d0 .or. p2(j, i) .eq. 0d0) then
          if (abs(p1(j, i) - p2(j, i)) .gt. vtiny) then
            momenta_equal = .false.
            return
          end if
        else
          if (abs((p1(j, i) - p2(j, i))/max(abs(p1(j, i)), abs(p2(j, i)))) .gt. vtiny) then
            momenta_equal = .false.
            return
          end if
        end if
      end do
    end do
  end function momenta_equal

  logical function momenta_equal_uborn(p1, p2, jfks1, ifks1, jfks2, ifks2)
! Returns .true. only if the momenta p1 and p2 are equal, but with the
! momenta of i_fks and j_fks summed. To save time, it only checks the
! 0th and 3rd components (energy and z-direction).
    implicit none
    integer i, j, jfks1, ifks1, jfks2, ifks2
    double precision p1(0:3, nexternal), p2(0:3, nexternal), pb1(0:3, nexternal), pb2(0:3, nexternal)
! Fill the underlying Born momenta pb1 and pb2
    do i = 1, nexternal
      do j = 0, 3, 3 ! skip x and y components, since they are not used in
! the 'momenta_equal' function
        if (i .lt. ifks1) then
          pb1(j, i) = p1(j, i)
        elseif (i .eq. ifks1) then
! Sum the i_fks to the j_fks momenta (i_fks is always greater than
! j_fks, so this is fine: it will NOT be overwritten later in the
! do-loop)
          pb1(j, jfks1) = pb1(j, jfks1) + p1(j, i)
          pb1(j, nexternal) = 0d0 ! fill the final one with zero's
        else
          pb1(j, i - 1) = p1(j, i)   ! skip the i_fks momenta
        end if
        if (i .lt. ifks2) then
          pb2(j, i) = p2(j, i)
        elseif (i .eq. ifks2) then
          pb2(j, jfks2) = pb2(j, jfks2) + p2(j, i) ! sum i_fks to j_fks momenta
          pb2(j, nexternal) = 0d0 ! fill the final one with zero's
        else
          pb2(j, i - 1) = p2(j, i)   ! skip the i_fks momenta
        end if
      end do
    end do
! Check if they are equal
    momenta_equal_uborn = momenta_equal(pb1, pb2)
  end function momenta_equal_uborn

  subroutine compute_prefactors_nbody(vegas_wgt)
! Compute all the relevant prefactors for the Born and the soft-virtual,
! i.e. all the nbody contributions. Also initialises the plots and
! bpower.
    use extra_weights
    use mint_module
    implicit none
    real :: tBefore, tAfter
    double precision pi, vegas_wgt
    integer i, j
    logical firsttime
    data firsttime/.true./
    parameter(pi=3.1415926535897932385d0)
    logical needrndec
    parameter(needrndec=.false.)
    integer orders(nsplitorders)
!
    call cpu_time(tBefore)
! Random numbers to be used in the plotting routine: these numbers will
! not change between events, counter events and n-body contributions.
    if (needrndec) then
      do i = 1, 10
        rndec(i) = ran2()
      end do
    end if
    if (firsttime) then
      if (pineappl) then
! PineAPPL stuff
        appl_amp_split_size = amp_split_size
        do j = 1, amp_split_size
          call amp_split_pos_to_orders(j, orders)
          appl_qcdpower(j) = orders(qcd_pos)
          appl_qedpower(j) = orders(qed_pos)
        end do
      end if
! Initialize hiostograms for fixed order runs
      if (fixed_order) call initplot_impl()
      firsttime = .false.
    end if
    call set_cms_stuff(0)
! f_* multiplication factors for Born and nbody
    f_b = jac_cnt(0)*xinorm_ev/(min(xiimax_ev, xiBSVcut_used)*shat/(16*pi**2))*fkssymmetryfactorBorn*vegas_wgt
    f_nb = f_b
    call cpu_time(tAfter)
    tf_nb = tf_nb + (tAfter - tBefore)
    return
  end subroutine compute_prefactors_nbody

  subroutine include_multichannel_enhance(imode)
    implicit none
    real :: tBefore, tAfter
    double precision xnoborn_cnt, xtot, wgt_c, enhance, enhance_real, pas(0:3, nexternal)
    data xnoborn_cnt/0d0/
    integer inoborn_cnt, i, imode
    data inoborn_cnt/0/
    double precision p_born_used(0:3, nexternal - 1)

    call cpu_time(tBefore)

! Compute the multi-channel enhancement factor 'enhance'.
    enhance = 1.d0
    if (p_born(0, 1) .gt. 0d0) then
      call sborn(p_born, wgt_c)
    elseif (p_born(0, 1) .lt. 0d0) then
      enhance = 0d0
    end if
    if (enhance .eq. 0d0) then
      xnoborn_cnt = xnoborn_cnt + 1.d0
      if (log10(xnoborn_cnt) .gt. inoborn_cnt) then
        write (*, *) 'WARNING: no Born momenta more than 10**', inoborn_cnt, 'times'
        inoborn_cnt = inoborn_cnt + 1
      end if
    else
      xtot = 0d0
      if (config_map(0, 0) .eq. 0) then
        write (*, *) 'Fatal error in compute_prefactor_nbody:'//' no Born diagrams ', config_map, '. Check bornfromreal.inc'
        write (*, *) 'Is fks_singular compiled correctly?'
        stop 1
      end if
      do i = 1, config_map(0, 0)
        xtot = xtot + amp2(config_map(i, 0))
      end do
      if (xtot .ne. 0d0) then
        enhance = amp2(config_map(this_config, 0))/xtot
        enhance = enhance*diagramsymmetryfactor
      else
        enhance = 0d0
      end if
    end if

! When not doing event projection, use the Born point computed by that
! mapping for the real multi-channel enhancement.
    enhance_real = 1.d0
    if (.not. use_evpr .and. imode .eq. 2) then
      p_born_used(:, :) = p_born_norad(:, :)
      if (p_born_used(0, 1) .gt. 0d0) then
        calculatedBorn = .false.
        pas(0:3, nexternal) = 0d0
        pas(0:3, 1:nexternal - 1) = p_born_used(0:3, 1:nexternal - 1)
        call set_alphas(pas)
        call sborn(p_born_used, wgt_c)
        call set_alphas(p_ev)
        calculatedBorn = .false.
      elseif (p_born_used(0, 1) .lt. 0d0) then
        if (enhance .ne. 0d0) then
          enhance_real = enhance
        else
          enhance_real = 0d0
        end if
      end if
! Compute the multi-channel enhancement factor 'enhance_real'.
      if (enhance_real .eq. 0d0) then
        xnoborn_cnt = xnoborn_cnt + 1.d0
        if (log10(xnoborn_cnt) .gt. inoborn_cnt) then
          write (*, *) 'WARNING: no Born momenta more than 10**', inoborn_cnt, 'times'
          inoborn_cnt = inoborn_cnt + 1
        end if
      else
        xtot = 0d0
        if (config_map(0, 0) .eq. 0) then
          write (*, *) 'Fatal error in compute_prefactor_n1body,'//' no Born diagrams ', config_map, '. Check bornfromreal.inc'
          write (*, *) 'Is fks_singular compiled correctly?'
          stop 1
        end if
        do i = 1, config_map(0, 0)
          xtot = xtot + amp2(config_map(i, 0))
        end do
        if (xtot .ne. 0d0) then
          enhance_real = amp2(config_map(this_config, 0))/xtot
          enhance_real = enhance_real*diagramsymmetryfactor
        else
          enhance_real = 0d0
        end if
      end if
    else
      enhance_real = enhance
    end if

    if (imode .eq. 1) then
      f_b = f_b*enhance
      f_nb = f_nb*enhance
    elseif (imode .eq. 2) then
      f_r = f_r*enhance_real
    elseif (imode .eq. 3) then
      f_s = f_s*enhance
      f_c = f_c*enhance
      f_dc = f_dc*enhance
      f_sc = f_sc*enhance
      f_dsc(1) = f_dsc(1)*enhance
      f_dsc(2) = f_dsc(2)*enhance
      f_dsc(3) = f_dsc(3)*enhance
      f_dsc(4) = f_dsc(4)*enhance
      f_pdfsch_d = f_pdfsch_d*enhance
      f_pdfsch_p = f_pdfsch_p*enhance
      f_pdfsch_l = f_pdfsch_l*enhance
    end if
    call cpu_time(tAfter)
    tf_nb = tf_nb + (tAfter - tBefore)

    return
  end subroutine include_multichannel_enhance

  subroutine compute_prefactors_n1body(vegas_wgt, jac_ev)
! Compute all relevant prefactors for the real emission and counter
! terms.
    implicit none
    real :: tBefore, tAfter
    double precision vegas_wgt, prefact, prefact_cnt_ssc, prefact_deg
    double precision prefact_c, prefact_coll, jac_ev, pi
    double precision prefact_cnt_ssc_c, prefact_coll_c
    double precision prefact_deg_slxi, prefact_deg_sxi
    integer i
    parameter(pi=3.1415926535897932385d0)
! prefactors for the PDF scheme
    double precision prefact_pdfsch_d, prefact_pdfsch_p, prefact_pdfsch_l
    double precision pmass(nexternal)
    pmass = external_masses
    call cpu_time(tBefore)

! f_* multiplication factors for real-emission, soft counter, ... etc.
    prefact = xinorm_ev/xi_i_fks_ev/(1 - y_ij_fks_ev)
    f_r = prefact*jac_ev*fkssymmetryfactor*vegas_wgt
    if (.not. nocntevents) then
      prefact_cnt_ssc = xinorm_ev/min(xiimax_ev, xiScut_used)*log(xicut_used/min(xiimax_ev, xiScut_used))/(1 - y_ij_fks_ev)
      f_s = (prefact + prefact_cnt_ssc)*jac_cnt(0)*fkssymmetryfactor*vegas_wgt
      if (pmass(j_fks) .eq. 0d0) then
! For the soft-collinear, these should be itwo. But they are always
! equal to ione, so no need to define separate factors.
        prefact_c = xinorm_cnt(1)/xi_i_fks_cnt(1)/(1 - y_ij_fks_ev)
        prefact_coll = xinorm_cnt(1)/xi_i_fks_cnt(1)*log(delta_used/deltaS)/deltaS
        f_c = (prefact_c + prefact_coll)*jac_cnt(1)*fkssymmetryfactor*vegas_wgt
        call set_cms_stuff(1)
        prefact_deg = xinorm_cnt(1)/xi_i_fks_cnt(1)/deltaS
        prefact_cnt_ssc_c = xinorm_cnt(1)/min(xiimax_cnt(1), xiScut_used) &
                            *log(xicut_used/min(xiimax_cnt(1), xiScut_used)) &
                            /(1 - y_ij_fks_ev)
        prefact_coll_c = xinorm_cnt(1)/min(xiimax_cnt(1), xiScut_used) &
                         *log(xicut_used/min(xiimax_cnt(1), xiScut_used)) &
                         *log(delta_used/deltaS)/deltaS
        f_dc = jac_cnt(1)*prefact_deg/(shat/(32*pi**2))*fkssymmetryfactorDeg*vegas_wgt
        f_sc = (prefact_c + prefact_coll + prefact_cnt_ssc_c + prefact_coll_c)*jac_cnt(2)*fkssymmetryfactorDeg*vegas_wgt
        call set_cms_stuff(2)
        prefact_deg_sxi = xinorm_cnt(1)/min(xiimax_cnt(1), xiScut_used)*log(xicut_used/min(xiimax_cnt(1), xiScut_used))*1/deltaS
        prefact_deg_slxi = xinorm_cnt(1)/min(xiimax_cnt(1), xiScut_used) &
                           *(log(xicut_used)**2 &
                             - log(min(xiimax_cnt(1), xiScut_used))**2) &
                           /(2.d0*deltaS)
        f_dsc(1) = prefact_deg*jac_cnt(2)/(shat/(32*pi**2))*fkssymmetryfactorDeg*vegas_wgt
        f_dsc(2) = prefact_deg_sxi*jac_cnt(2)/(shat/(32*pi**2))*fkssymmetryfactorDeg*vegas_wgt
        f_dsc(3) = prefact_deg_slxi*jac_cnt(2)/(shat/(32*pi**2))*fkssymmetryfactorDeg*vegas_wgt
        f_dsc(4) = (prefact_deg + prefact_deg_sxi)*jac_cnt(2)/(shat/(32*pi**2))*fkssymmetryfactorDeg*vegas_wgt
! prefactor for the PDF scheme
        prefact_pdfsch_d = xinorm_cnt(1)/xiScut_used/deltaS
        f_pdfsch_d = prefact_pdfsch_d*jac_cnt(2)/(shat/(32*pi**2))*fkssymmetryfactorDeg*vegas_wgt
        prefact_pdfsch_p = xinorm_cnt(1)*dlog(xiScut_used)/xiScut_used/deltaS
        f_pdfsch_p = prefact_pdfsch_p*jac_cnt(2)/(shat/(32*pi**2))*fkssymmetryfactorDeg*vegas_wgt
        prefact_pdfsch_l = xinorm_cnt(1)*dlog(xiScut_used)**2/2d0/xiScut_used/deltaS
        f_pdfsch_l = prefact_pdfsch_l*jac_cnt(2)/(shat/(32*pi**2))*fkssymmetryfactorDeg*vegas_wgt
      else
        f_c = 0d0
        f_dc = 0d0
        f_sc = 0d0
        do i = 1, 4
          f_dsc(i) = 0d0
        end do
      end if
    else
      f_s = 0d0
      f_c = 0d0
      f_dc = 0d0
      f_sc = 0d0
      do i = 1, 4
        f_dsc(i) = 0d0
      end do
    end if
    call cpu_time(tAfter)
    tf_all = tf_all + (tAfter - tBefore)
    return
  end subroutine compute_prefactors_n1body

  subroutine add_wgt(type, orders, wgt1, wgt2, wgt3)
! Adds a contribution to the list in weight_lines. 'type' sets the type
! of the contribution and wgt1..wgt3 are the coefficients multiplying
! the logs. The arguments are:
!     type=1 : real-emission
!     type=2 : Born
!     type=3 : integrated counter terms
!     type=4 : soft counter-term
!     type=5 : collinear counter-term
!     type=6 : soft-collinear counter-term
!     type=14: virtual corrections
!     type=15: virt-trick: average born contribution
!     wgt1 : weight of the contribution not multiplying a scale log
!     wgt2 : coefficient of the weight multiplying the log(/mu_R^2/Q^2/)
!     wgt3 : coefficient of the weight multiplying the log(/mu_F^2/Q^2/)
!
!
! The argument orders specifies what are the squared coupling orders
! factorizing the particular set of weights added here. The position of
! the QCD order there can be obtained via qcd_pos from orders.inc
! This is solely used for now in order to apply a potential user-defined filer.
!
! This subroutine increments the 'icontr' counter: each new call to this
! function makes sure that it's considered a new contribution. For each
! contribution, we save the
!     The type: itype(icontr)
!     The weights: wgt(1,icontr),wgt(2,icontr) and wgt(3,icontr) for
!         wgt1, wgt2 and wgt3, respectively.
!     The Bjorken x's: bjx(1,icontr), bjx(2,icontr)
!     The Ellis-Sexton scale squared used to compute the weight:
!        scales2(1,icontr)
!     The renormalisation scale squared used to compute the weight:
!        scales2(2,icontr)
!     The factorisation scale squared used to compute the weight:
!       scales2(3,icontr)
!     The value of the strong coupling: g_strong(icontr)
!     The FKS configuration: nFKS(icontr)
!     The boost to go from the momenta in the C.o.M. frame to the
!         laboratory frame: y_bst(icontr)
!     The power of the strong coupling (g_strong) for the current
!       weight: QCDpower(icontr)
!     The momenta: momenta(j,i,icontr). For the Born contribution, the
!        counter-term momenta are used. This is okay for any IR-safe
!        observables.
!     The PDG codes: pdg(i,icontr). Always the ones with length
!        'nexternal' are used, because the momenta are also the
!        'nexternal' ones. This is okay for IR-safe observables.
!     The PDG codes of the underlying Born process:
!        pdg_uborn(i,icontr). The PDGs of j_fks and i_fks are combined
!        to get the PDG code of the mother. The extra parton is given a
!        PDG=21 (gluon) code.
!     If the contribution belongs to an H-event or S-event:
!        H_event(icontr)
!     The weight of the born or real-emission matrix element
!        corresponding to this contribution: wgt_ME_tree. This weight does
!        include the 'ngluon' correction factor for the Born.
!
! Not set in this subroutine, but included in the weight_lines module
! are the
!     wgts(iwgt,icontr) : weights including scale/PDFs/logs. These are
!        normalised so that they can be used directly to compute cross
!        sections and fill plots. 'iwgt' goes from 1 to the maximum
!        number of weights obtained from scale and PDF reweighting, with
!        the iwgt=1 element being the central value.
!     plot_id(icontr) : =20 for Born, 11 for real-emission and 12 for
!        anything else.
!     plot_wgts(iwgt,icontr) : same as wgts(), but only non-zero for
!        unique contributions and non-unique are added to the unique
!        ones. 'Unique' here is defined that they would be identical in
!        an analysis routine (i.e. same momenta and PDG codes)
!     niproc(icontr) : number of combined subprocesses in parton_lum_*.f
!     parton_iproc(iproc,icontr) : value of the PDF for the iproc
!        contribution
!     parton_pdg(nexternal,iproc,icontr) : value of the PDG codes for
!     the iproc contribution
!     ipr(icontr): for separate_flavour_configs: the iproc of current
!        contribution
    use weight_lines
    use extra_weights
    use FKSParams
    implicit none
    integer type, i, j
    logical foundIt, foundOrders
    double precision wgt1, wgt2, wgt3
    integer orders(nsplitorders)

    if (wgt1 .eq. 0d0 .and. wgt2 .eq. 0d0 .and. wgt3 .eq. 0d0) return
! Check for NaN's and INF's. Simply skip the contribution
    if (wgt1 .ne. wgt1) return
    if (wgt2 .ne. wgt2) return
    if (wgt3 .ne. wgt3) return

! Apply user-defined (in FKS_params.dat) contribution type filters if necessary

    if (SelectedContributionTypes(0) .gt. 0) then
      foundIt = .false.
      do i = 1, SelectedContributionTypes(0)
        if (type .eq. SelectedContributionTypes(i)) then
          foundIt = .true.
          exit
        end if
      end do
      if (.not. foundIt) then
! This contribution was not part of the user selection. Skip it.
        return
      end if
    end if

! Apply the user-defined coupling-order filter if present
! First the simple QCD filter
    if (QCD_squared_selected .ne. -1 .and. QCD_squared_selected .ne. orders(qcd_pos)) then
      return
    end if
! Secondly, the more advanced filter
    if (SelectedCouplingOrders(1, 0) .gt. 0) then
      foundIt = .false.
      do j = 1, SelectedCouplingOrders(1, 0)
        foundOrders = .true.
        do i = 1, nsplitorders
          if (SelectedCouplingOrders(i, j) .ne. orders(i)) then
            foundOrders = .false.
            exit
          end if
        end do
        if (foundOrders) then
          foundIt = .true.
          exit
        end if
      end do
      if (.not. foundIt) then
        return
      end if
    end if

    icontr = icontr + 1
    call weight_lines_allocated(nexternal, icontr, max_wgt, max_iproc)
    itype(icontr) = type

    wgt(1, icontr) = wgt1
    wgt(2, icontr) = wgt2
    wgt(3, icontr) = wgt3

    bjx(1, icontr) = xbk(1)
    bjx(2, icontr) = xbk(2)
    scales2(1, icontr) = QES2
    scales2(2, icontr) = scale**2
    scales2(3, icontr) = q2fact(1)
    g_strong(icontr) = g
    nFKS(icontr) = nFKSprocess
    y_bst(icontr) = ybst_til_tolab
    ifold_cnt(icontr) = ifold_counter
    qcdpower(icontr) = QCD_power
    orderstag(icontr) = orders_tag
    amppos(icontr) = amp_pos
    ipr(icontr) = 0
    call set_pdg_impl(icontr, nFKSprocess, idup)

! Compensate for the fact that in the Born matrix elements, we use the
! identical particle symmetry factor of the corresponding real emission
! matrix elements
    wgt_ME_tree(1, icontr) = wgt_me_born
    wgt_ME_tree(2, icontr) = wgt_me_real
    do i = 1, nexternal
      do j = 0, 3
        if (p1_cnt(0, 1, 0) .gt. 0d0 .and. type .ne. 5) then
          momenta_m(j, i, 1, icontr) = p1_cnt(j, i, 0)
        elseif (p1_cnt(0, 1, 1) .gt. 0d0) then
          momenta_m(j, i, 1, icontr) = p1_cnt(j, i, 1)
        elseif (p1_cnt(0, 1, 2) .gt. 0d0) then
          momenta_m(j, i, 1, icontr) = p1_cnt(j, i, 2)
        else
          if (i .lt. fks_i_d(nFKSprocess)) then
            momenta_m(j, i, 1, icontr) = p_born(j, i)
          elseif (i .eq. fks_i_d(nFKSprocess)) then
            momenta_m(j, i, 1, icontr) = 0d0
          else
            momenta_m(j, i, 1, icontr) = p_born(j, i - 1)
          end if
        end if
        momenta_m(j, i, 2, icontr) = p_ev(j, i)
      end do
    end do

    if (type .eq. 1) then
! Real-emission contribution with n+1-body kinematics.
      do i = 1, nexternal
        do j = 0, 3
          momenta(j, i, icontr) = momenta_m(j, i, 2, icontr)
        end do
      end do
      H_event(icontr) = .true.
    elseif (type .ge. 2 .and. type .le. 6 .or. type .eq. 14 .or. type .eq. 15 .or. (type .ge. 20 .and. type .le. 22)) then
! Born, counter term, soft-virtual, or n-body real contributions.
      do i = 1, nexternal
        do j = 0, 3
          momenta(j, i, icontr) = momenta_m(j, i, 1, icontr)
        end do
      end do
      H_event(icontr) = .false.
    else
      write (*, *) 'ERROR: unknown type in add_wgt', type
      stop 1
    end if
    return
  end subroutine add_wgt
  subroutine include_PDF_and_alphas
! Multiply the saved wgt() info by the PDFs, alpha_S and the scale
! dependence and saves the weights in the wgts() array. The weights in
! this array are now correctly normalised to compute the cross section
! or to fill histograms.
    use weight_lines
    use extra_weights
    use mint_module
    use FKSParams
    implicit none
    real :: tBefore, tAfter
    integer orders(nsplitorders)
    integer i, j, iamp, icontr_orig
    logical virt_found
    double precision xlum, dlum, mu2_r, mu2_f, mu2_q, wgt_wo_pdf, conv
    external dlum
    parameter(conv=389379660d0) ! conversion to picobarns
    call cpu_time(tBefore)
    if (icontr .eq. 0) return
    virt_found = .false.
! number of contributions before they are (possibly) increased through a
! call to separate_flavour_config().
    icontr_orig = icontr
    i = 0
    do while (i .lt. icontr)
      i = i + 1
      nFKSprocess = nFKS(i)
      xbk(1) = bjx(1, i)
      xbk(2) = bjx(2, i)
      mu2_q = scales2(1, i)
      mu2_r = scales2(2, i)
      mu2_f = scales2(3, i)
      q2fact(1) = mu2_f
      q2fact(2) = mu2_f
! call the PDFs
      xlum = dlum()
! iwgt=1 is the central value (i.e. no scale/PDF reweighting).
      iwgt = 1
      call weight_lines_allocated(nexternal, max_contr, iwgt, subproc_iproc)
! set_pdg_codes fills the niproc, parton_iproc, parton_pdg and
! parton_pdg_uborn [Do only for the contributions that were already
! available as part of the input -- NOT the ones that are created
! through the call to separate_flavour_config(), since that will
! overwrite the relevant information.]
      if (i .le. icontr_orig) call set_pdg_codes(subproc_iproc, subproc_pd, nFKSprocess, i)
      if (separate_flavour_configs .and. ipr(i) .eq. 0) then
        call separate_flavour_config(i) ! this increases icontr
      end if
      if (separate_flavour_configs .and. ipr(i) .ne. 0) then
        if (nincoming .eq. 2) then
          xlum = subproc_pd(ipr(i))*conv
        else
          xlum = subproc_pd(ipr(i))
        end if
      end if
      wgt_wo_pdf = (wgt(1, i) + wgt(2, i)*log(mu2_r/mu2_q) &
                    + wgt(3, i)*log(mu2_f/mu2_q))*g_strong(i)**QCDpower(i)
      wgts(iwgt, i) = xlum*wgt_wo_pdf
      do j = 1, subproc_iproc
        parton_iproc(j, i) = parton_iproc(j, i)*wgt_wo_pdf
      end do
      if (itype(i) .eq. 14 .and. .not. virt_found) then
        virt_found = .true.
! Special for the soft-virtual needed for the virt-tricks. The
! *_wgt_mint variable should be directly passed to the mint-integrator
! and not be part of the plots nor computation of the cross section.
        virt_wgt_mint(0) = virt_wgt_mint(0)*xlum
        born_wgt_mint(0) = born_wgt_mint(0)*xlum
        do iamp = 1, amp_split_size
          call amp_split_pos_to_orders(iamp, orders)
          QCD_power = orders(qcd_pos)
          virt_wgt_mint(iamp) = virt_wgt_mint(iamp)*xlum
          born_wgt_mint(iamp) = born_wgt_mint(iamp)*xlum
        end do
      end if
    end do
    call cpu_time(tAfter)
    t_as = t_as + (tAfter - tBefore)
    return
  end subroutine include_PDF_and_alphas

  subroutine separate_flavour_config(ict)
    use weight_lines
    implicit none
    integer ict, i_add, i, j, k, ict_new, n
    if ((.not. fixed_order) .or. niproc(ict) .eq. 1) then
      return
    end if
    i_add = niproc(ict) - 1
    call weight_lines_allocated(nexternal, icontr + i_add, max_wgt, max_iproc)
    do i = 1, niproc(ict)
      if (i .eq. 1) then
        niproc(ict) = 1
        ipr(ict) = 1
        cycle
      end if
      ict_new = icontr + (i - 1)
      ipr(ict_new) = i
      itype(ict_new) = itype(ict)
      do j = 1, 3
        wgt(j, ict_new) = wgt(j, ict)
        scales2(j, ict_new) = scales2(j, ict)
      end do
      do j = 1, 2
        bjx(j, ict_new) = bjx(j, ict)
        wgt_ME_tree(j, ict_new) = wgt_ME_tree(j, ict)
      end do
      g_strong(ict_new) = g_strong(ict)
      nFKS(ict_new) = nFKS(ict)
      y_bst(ict_new) = y_bst(ict)
      QCDpower(ict_new) = QCDpower(ict)
      orderstag(ict_new) = orderstag(ict)
      H_event(ict_new) = H_event(ict)
      do k = 1, nexternal
        do j = 0, 3
          momenta(j, k, ict_new) = momenta(j, k, ict)
          do n = 1, 2
            momenta_m(j, k, n, ict_new) = momenta_m(j, k, n, ict)
          end do
        end do
        pdg(k, ict_new) = parton_pdg(k, i, ict)
        pdg_uborn(k, ict_new) = parton_pdg_uborn(k, i, ict)
        parton_pdg(k, 1, ict_new) = parton_pdg(k, i, ict)
        parton_pdg_uborn(k, 1, ict_new) = parton_pdg_uborn(k, i, ict)
      end do
      niproc(ict_new) = 1
      parton_iproc(1, ict_new) = parton_iproc(i, ict)
    end do
    icontr = icontr + i_add
    return
  end subroutine separate_flavour_config

  subroutine set_pdg_codes(iproc, pd, iFKS, ict)
    use weight_lines
    implicit none
    integer j, k, iproc, ict, iFKS
    double precision pd(0:maxproc), conv
    parameter(conv=389379660d0) ! conversion to picobarns
! save also the separate contributions to the PDFs and the corresponding
! PDG codes
    niproc(ict) = iproc
    do j = 1, iproc
      if (nincoming .eq. 2) then
        parton_iproc(j, ict) = pd(j)*conv
      else
!           Keep GeV's for decay processes (no conv. factor needed)
        parton_iproc(j, ict) = pd(j)
      end if
      do k = 1, nexternal
        parton_pdg(k, j, ict) = idup_d(iFKS, k, j)
        if (k .lt. fks_j_d(iFKS)) then
          parton_pdg_uborn(k, j, ict) = idup_d(iFKS, k, j)
        elseif (k .eq. fks_j_d(iFKS)) then
          if (abs(idup_d(iFKS, fks_i_d(iFKS), j)) .eq. &
              abs(idup_d(iFKS, fks_j_d(iFKS), j)) .and. &
              abs(pdg(fks_i_d(iFKS), ict)) .ne. 21) then
! check if any extra cnt is needed
            if (extra_cnt_d(iFKS) .eq. 0) then
! if not, assign a gluon for the QCD splitting
              if (split_type_d(iFKS, qcd_pos)) then
                parton_pdg_uborn(k, j, ict) = 21
              else
                write (*, *) 'set_pdg_codes ', 'ERROR#1 in PDG assigment for underlying Born'
                stop 1
              end if
            else
! if there are extra cnt's, assign the pdg of the
! mother in the born (according to isplitorder_born_d)
              if (isplitorder_born_d(iFKS) .eq. qcd_pos) then
                parton_pdg_uborn(k, j, ict) = 21
              else
                write (*, *) 'set_pdg_codes ', 'ERROR#2 in PDG assigment for underlying Born'
                stop 1
              end if
            end if
          elseif (abs(idup_d(iFKS, fks_i_d(iFKS), j)) .eq. 21) then
            parton_pdg_uborn(k, j, ict) = idup_d(iFKS, fks_j_d(iFKS), j)
          elseif (idup_d(iFKS, fks_j_d(iFKS), j) .eq. 21) then
            parton_pdg_uborn(k, j, ict) = -idup_d(iFKS, fks_i_d(iFKS), j)
          else
            write (*, *) 'set_pdg_codes ', 'ERROR#3 in PDG assigment for underlying Born'
            stop 1
          end if
        elseif (k .lt. fks_i_d(iFKS)) then
          parton_pdg_uborn(k, j, ict) = idup_d(iFKS, k, j)
        elseif (k .eq. nexternal) then
          if (split_type_d(iFKS, qcd_pos)) then
            parton_pdg_uborn(k, j, ict) = 21 ! give the extra particle a gluon PDG code
          end if
        elseif (k .ge. fks_i_d(iFKS)) then
          parton_pdg_uborn(k, j, ict) = idup_d(iFKS, k + 1, j)
        end if
      end do
    end do
    return
  end subroutine set_pdg_codes

  subroutine reweight_scale
! Use the saved weight_lines info to perform scale reweighting. Extends the
! wgts() array to include the weights.
    use weight_lines
    use extra_weights
    use FKSParams
    use alfas_functions_module, only: alphas
    implicit none
    real :: tBefore, tAfter
    integer i, kr, kf, iwgt_save, dd
    double precision xlum(maxscales), dlum, pi, mu2_r(maxscales), c_mu2_r, c_mu2_f, mu2_f(maxscales), mu2_q, g(maxscales), conv
    parameter(pi=3.1415926535897932385d0)
    external dlum
    parameter(conv=389379660d0) ! conversion to picobarns
    call cpu_time(tBefore)
    if (icontr .eq. 0) return
! currently we have 'iwgt' weights in the wgts() array.
    iwgt_save = iwgt
! loop over all the contributions in the weight lines module
    do i = 1, icontr
      iwgt = iwgt_save
      nFKSprocess = nFKS(i)
      xbk(1) = bjx(1, i)
      xbk(2) = bjx(2, i)
      mu2_q = scales2(1, i)
! Loop over the dynamical_scale_choices
      do dd = 1, dyn_scale(0)
! renormalisation scale variation (requires recomputation of the strong
! coupling)
        call set_mu_central(i, dd, c_mu2_r, c_mu2_f)
        do kr = 1, nint(scalevarR(0))
          if ((.not. lscalevar(dd)) .and. kr .ne. 1) exit
          mu2_r(kr) = c_mu2_r*scalevarR(kr)**2
          g(kr) = sqrt(4d0*pi*alphas(sqrt(mu2_r(kr))))
        end do
! factorisation scale variation (require recomputation of the PDFs)
        do kf = 1, nint(scalevarF(0))
          if ((.not. lscalevar(dd)) .and. kf .ne. 1) exit
          mu2_f(kf) = c_mu2_f*scalevarF(kf)**2
          q2fact(1) = mu2_f(kf)
          q2fact(2) = mu2_f(kf)
          xlum(kf) = dlum()
          if (separate_flavour_configs .and. ipr(i) .ne. 0) then
            if (nincoming .eq. 2) then
              xlum(kf) = subproc_pd(ipr(i))*conv
            else
              xlum(kf) = subproc_pd(ipr(i))
            end if
          end if
        end do
        do kf = 1, nint(scalevarF(0))
          if ((.not. lscalevar(dd)) .and. kf .ne. 1) exit
          do kr = 1, nint(scalevarR(0))
            if ((.not. lscalevar(dd)) .and. kr .ne. 1) exit
            iwgt = iwgt + 1   ! increment the iwgt for the wgts() array
            call weight_lines_allocated(nexternal, max_contr, iwgt, max_iproc)
! add the weights to the array
            wgts(iwgt, i) = xlum(kf) &
                            *(wgt(1, i) + wgt(2, i)*log(mu2_r(kr)/mu2_q) &
                              + wgt(3, i)*log(mu2_f(kf)/mu2_q)) &
                            *g(kr)**QCDpower(i)
          end do
        end do
      end do
    end do
    call cpu_time(tAfter)
    tr_s = tr_s + (tAfter - tBefore)
    return
  end subroutine reweight_scale

  subroutine reweight_pdf
! Use the saved weight_lines info to perform PDF reweighting. Extends the
! wgts() array to include the weights.
    use weight_lines
    use extra_weights
    use FKSParams
    use alfas_functions_module, only: alphas
    implicit none
    real :: tBefore, tAfter
    integer n, i, nn
    double precision xlum, dlum, pi, mu2_r, mu2_f, mu2_q, g, conv
    parameter(pi=3.1415926535897932385d0)
    external dlum
    parameter(conv=389379660d0) ! conversion to picobarns
    call cpu_time(tBefore)
    if (icontr .eq. 0) return
    do nn = 1, lhaPDFid(0)
! Use as external loop the one over the PDF sets and as internal the one
! over the icontr. This reduces the number of calls to InitPDF and
! allows for better caching of the PDFs
      do n = 0, nmemPDF(nn)
        iwgt = iwgt + 1
        call weight_lines_allocated(nexternal, max_contr, iwgt, max_iproc)
        call InitPDFm(nn, n)
        do i = 1, icontr
          nFKSprocess = nFKS(i)
          xbk(1) = bjx(1, i)
          xbk(2) = bjx(2, i)
          mu2_q = scales2(1, i)
          mu2_r = scales2(2, i)
          mu2_f = scales2(3, i)
          q2fact(1) = mu2_f
          q2fact(2) = mu2_f
! Compute the luminosity
          xlum = dlum()
          if (separate_flavour_configs .and. ipr(i) .ne. 0) then
            if (nincoming .eq. 2) then
              xlum = subproc_pd(ipr(i))*conv
            else
              xlum = subproc_pd(ipr(i))
            end if
          end if
! Recompute the strong coupling: alpha_s in the PDF might change
          g = sqrt(4d0*pi*alphas(sqrt(mu2_r)))
! add the weights to the array
          wgts(iwgt, i) = xlum*(wgt(1, i) + wgt(2, i)*log(mu2_r/mu2_q) + wgt(3, i)*log(mu2_f/mu2_q))*g**QCDpower(i)
        end do
      end do
    end do
    call InitPDFm(1, 0)
    call cpu_time(tAfter)
    tr_pdf = tr_pdf + (tAfter - tBefore)
    return
  end subroutine reweight_pdf

  subroutine fill_pineappl_weights(vegas_wgt)
! Fill the PineAPPL state owned by fnlo_process_common. This subroutine assumes
! that there is an unique PS configuration: at most one Born, one real
! and one set of counter terms. Among other things, this means that one
! must do MC over FKS directories.
    use weight_lines
    implicit none
    integer i, j
    double precision final_state_rescaling, vegas_wgt
    integer pos
    do i = 1, 4
      do j = 1, amp_split_size
        appl_w0(i, j) = 0d0
        appl_wR(i, j) = 0d0
        appl_wF(i, j) = 0d0
        appl_wB(i, j) = 0d0
      end do
      appl_x1(i) = 0d0
      appl_x2(i) = 0d0
      appl_QES2(i) = 0d0
      appl_muR2(i) = 0d0
      appl_muF2(i) = 0d0
      appl_flavmap(i) = 0
    end do
    appl_event_weight = 0d0
    appl_vegaswgt = vegas_wgt
    if (icontr .eq. 0) return
    do i = 1, icontr
      appl_event_weight = appl_event_weight + wgts(1, i)/vegas_wgt
      final_state_rescaling = dble(iproc_save(nFKS(i)))/dble(appl_nproc(flavour_map(nFKS(i))))
      if (itype(i) .eq. 2) then
        pos = lo_qcd_to_amp_pos(qcdpower(i))
      else
        pos = nlo_qcd_to_amp_pos(qcdpower(i))
      end if
! consistency check
      if (appl_qcdpower(pos) .ne. qcdpower(i)) then
        write (*, *) 'ERROR in fill_pineappl_weights, QCDpower', appl_qcdpower(pos), qcdpower(i)
        stop 1
      end if

      if (itype(i) .eq. 1) then
!     real
        appl_w0(1, pos) = appl_w0(1, pos) + wgt(1, i)*final_state_rescaling
        appl_x1(1) = bjx(1, i)
        appl_x2(1) = bjx(2, i)
        appl_flavmap(1) = flavour_map(nFKS(i))
        appl_QES2(1) = scales2(1, i)
        appl_muR2(1) = scales2(2, i)
        appl_muF2(1) = scales2(3, i)
      elseif (itype(i) .eq. 2) then
!     born
        appl_wB(2, pos) = appl_wB(2, pos) + wgt(1, i)*final_state_rescaling
        appl_x1(2) = bjx(1, i)
        appl_x2(2) = bjx(2, i)
        appl_flavmap(2) = flavour_map(nFKS(i))
        appl_QES2(2) = scales2(1, i)
        appl_muR2(2) = scales2(2, i)
        appl_muF2(2) = scales2(3, i)
      elseif (itype(i) .eq. 3 .or. itype(i) .eq. 4 .or. itype(i) .eq. 14 .or. itype(i) .eq. 15) then
!     virtual, soft-virtual or soft-counter
        appl_w0(2, pos) = appl_w0(2, pos) + wgt(1, i)*final_state_rescaling
        appl_wR(2, pos) = appl_wR(2, pos) + wgt(2, i)*final_state_rescaling
        appl_wF(2, pos) = appl_wF(2, pos) + wgt(3, i)*final_state_rescaling
        appl_x1(2) = bjx(1, i)
        appl_x2(2) = bjx(2, i)
        appl_flavmap(2) = flavour_map(nFKS(i))
        appl_QES2(2) = scales2(1, i)
        appl_muR2(2) = scales2(2, i)
        appl_muF2(2) = scales2(3, i)
      elseif (itype(i) .eq. 5) then
!     collinear counter
        appl_w0(3, pos) = appl_w0(3, pos) + wgt(1, i)*final_state_rescaling
        appl_wF(3, pos) = appl_wF(3, pos) + wgt(3, i)*final_state_rescaling
        appl_x1(3) = bjx(1, i)
        appl_x2(3) = bjx(2, i)
        appl_flavmap(3) = flavour_map(nFKS(i))
        appl_QES2(3) = scales2(1, i)
        appl_muR2(3) = scales2(2, i)
        appl_muF2(3) = scales2(3, i)
      elseif (itype(i) .eq. 6) then
!     soft-collinear counter
        appl_w0(4, pos) = appl_w0(4, pos) + wgt(1, i)*final_state_rescaling
        appl_wF(4, pos) = appl_wF(4, pos) + wgt(3, i)*final_state_rescaling
        appl_x1(4) = bjx(1, i)
        appl_x2(4) = bjx(2, i)
        appl_flavmap(4) = flavour_map(nFKS(i))
        appl_QES2(4) = scales2(1, i)
        appl_muR2(4) = scales2(2, i)
        appl_muF2(4) = scales2(3, i)
      else
        write (*, *) 'ERROR in fill_applgrid_weights', itype(i)
        stop 1
      end if
    end do
    return
  end subroutine fill_pineappl_weights

  subroutine get_wgt_nbody(sig)
! Sums all the central weights that contribution to the nbody cross
! section
    use weight_lines
    implicit none
    double precision sig
    integer i
    sig = 0d0
    if (icontr .eq. 0) return
    do i = 1, icontr
      if (itype(i) .eq. 2 .or. itype(i) .eq. 3 .or. itype(i) .eq. 14 .or. itype(i) .eq. 15) then
        sig = sig + wgts(1, i)
      end if
    end do
    return
  end subroutine get_wgt_nbody

  subroutine get_wgt_no_nbody(sig)
! Sums all the central weights that contribution to the cross section
! excluding the nbody contributions.
    use weight_lines
    implicit none
    double precision sig
    integer i
    sig = 0d0
    if (icontr .eq. 0) return
    do i = 1, icontr
      if (itype(i) .ne. 2 .and. itype(i) .ne. 3 .and. itype(i) .ne. 14 .and. itype(i) .ne. 15) then
        sig = sig + wgts(1, i)
      end if
    end do
    return
  end subroutine get_wgt_no_nbody

  subroutine fill_plots
! Calls the analysis routine (which fill plots) for all the
! contributions in the weight_lines module. Instead of really calling
! it for all, it first checks if weights can be summed (i.e. they have
! the same PDG codes and the same momenta) before calling the analysis
! to greatly reduce the calls to the analysis routines.
    use weight_lines
    use extra_weights
    implicit none
    real :: tBefore, tAfter
    integer i, ii, j, max_weight
    double precision, allocatable :: www(:)
! stuff for plotting the different splitorders
    save max_weight
    call cpu_time(tBefore)
    if (icontr .eq. 0) return
! fill the plots_wgts. Check if we can sum weights together before
! calling the analysis routines. This is the case if the PDG codes and
! the momenta are identical.
    do i = 1, icontr
      do j = 1, iwgt
        plot_wgts(j, i) = 0d0
      end do
! The following if lines have been changed with respect to the
! usual (with just 3 plot ids: 20, 11 and 12):
!  The kinematics of soft and collinear counterterms may
!  be different, for those processes without soft singularities
!  from initial(final)-state configurations when the
!  final(initial) confs are integrated (e.g. a a > e+ e-)
!  This gives no problem for normal histogramming (and in
!  fact plot_id 11 13 and 14 are merged into ibody=2 in
!  outfun, but it gives troubles e.g. with applgrid/pineappl.
!  Note that the separation between soft and soft-virtual
!  may not be needed in reality
      if (itype(i) .eq. 2) then
        plot_id(i) = 20 ! Born
      elseif (itype(i) .eq. 1) then
        plot_id(i) = 11 ! real-emission
      elseif (itype(i) .eq. 5) then
        plot_id(i) = 13 ! collinear counter term
      elseif (itype(i) .eq. 6) then
        plot_id(i) = 14 ! soft collinear counter term
      else
        plot_id(i) = 12 ! soft-virtual and soft counter term
      end if
! Loop over all previous icontr. If the plot_id, PDGs and momenta are
! equal to a previous icountr, add the current weight to the plot_wgts
! of that contribution and exit the do-loop. This loop extends to 'i',
! so if the current weight cannot be summed to a previous one, the ii=i
! contribution makes sure that it is added as a new element.
      do ii = 1, i
        if (orderstag(ii) .ne. orderstag(i)) cycle
        if (plot_id(ii) .ne. plot_id(i)) cycle
        if (plot_id(i) .ne. 11) then
          if (.not. pdg_equal(pdg_uborn(1, ii), pdg_uborn(1, i))) cycle
        else
          if (.not. pdg_equal(pdg(1, ii), pdg(1, i))) cycle
        end if
        if (plot_id(i) .ne. 11) then
          if (.not. momenta_equal_uborn( &
              momenta(0, 1, ii), momenta(0, 1, i), &
              fks_j_d(nFKS(ii)), fks_i_d(nFKS(ii)), &
              fks_j_d(nFKS(i)), fks_i_d(nFKS(i)))) cycle
        else
          if (.not. momenta_equal(momenta(0, 1, ii), momenta(0, 1, i))) cycle
        end if
        do j = 1, iwgt
          plot_wgts(j, ii) = plot_wgts(j, ii) + wgts(j, i)
        end do
        exit
      end do
    end do
    do i = 1, icontr
      if (plot_wgts(1, i) .ne. 0d0) then
        if (.not. allocated(www)) then
          allocate (www(iwgt))
          max_weight = iwgt
        elseif (iwgt .ne. max_weight) then
          write (*, *) 'Error in fill_plots (fks_singular.f): '// &
            'number of weights should not vary between PS points', &
            iwgt, max_weight
          stop
        end if
        do j = 1, iwgt
          www(j) = plot_wgts(j, i)
        end do
! call the analysis/histogramming routines
        orders_tag_plot = orderstag(i)
        amp_pos_plot = amppos(i)
        call outfun_impl(momenta(0, 1, i), y_bst(i), www, pdg(1, i), plot_id(i), &
                         external_masses)
      end if
    end do
    call cpu_time(tAfter)
    t_plot = t_plot + (tAfter - tBefore)
    return
  end subroutine fill_plots

  subroutine fill_mint_function(f)
! Fills the function that is returned to the MINT integrator
    use weight_lines
    use mint_module
    implicit none
    integer i, iamp, ithree, isix
    double precision f(nintegrals), sigint
    sigint = 0d0
    do i = 1, icontr
      sigint = sigint + wgts(1, i)
    end do
    f(1) = abs(sigint)
    f(2) = sigint
    f(4) = virtual_over_born    ! not used for anything
    do iamp = 0, amp_split_size
      if (iamp .eq. 0) then
        f(3) = 0d0
        f(6) = 0d0
        f(5) = 0d0
        do i = 1, amp_split_size
          f(3) = f(3) + virt_wgt_mint(i)
          f(6) = f(6) + born_wgt_mint(i)
        end do
        f(5) = abs(f(3))!v3.5.4, this fixes a wrong behaviour
      else
        ithree = 2*iamp + 5
        isix = 2*iamp + 6
        f(ithree) = virt_wgt_mint(iamp)
        f(isix) = born_wgt_mint(iamp)
      end if
    end do
    return
  end subroutine fill_mint_function

  subroutine rotate_invar(pin, pout, cth, sth, cphi, sphi)
! Given the four momentum pin, returns the four momentum pout (in the same
! Lorentz frame) by performing a three-rotation of an angle theta
! (cos(theta)=cth) around the y axis, followed by a three-rotation of an
! angle phi (cos(phi)=cphi) along the z axis. The components of pin
! and pout are given along these axes
    implicit none
    double precision cth, sth, cphi, sphi, pin(0:3), pout(0:3)
    double precision q1, q2, q3
!
    q1 = pin(1)
    q2 = pin(2)
    q3 = pin(3)
    pout(1) = q1*cphi*cth - q2*sphi + q3*cphi*sth
    pout(2) = q1*sphi*cth + q2*cphi + q3*sphi*sth
    pout(3) = -q1*sth + q3*cth
    pout(0) = pin(0)
    return
  end subroutine rotate_invar

  subroutine getaziangles(p, cphi, sphi)
    implicit none
    double precision p(0:3), cphi, sphi
    double precision xlength, cth, sth
!
    xlength = rho(p)
    if (xlength .ne. 0.d0) then
      cth = p(3)/xlength
      sth = sqrt(1 - cth**2)
      if (sth .ne. 0.d0) then
        cphi = p(1)/(xlength*sth)
        sphi = p(2)/(xlength*sth)
      else
        cphi = 1.d0
        sphi = 0.d0
      end if
    else
      cphi = 1.d0
      sphi = 0.d0
    end if
    return
  end subroutine getaziangles

  subroutine phspncheck_born(ecm, xmass, xmom, pass)
! Checks four-momentum conservation.
! WARNING: works only in the partonic c.m. frame
    implicit none
    double precision ecm, xmass(nexternal - 1), xmom(0:3, nexternal - 1)
    double precision tiny, xm, xsum(0:3), xsuma(0:3), xrat(0:3), ptmp(0:3)
    parameter(tiny=5.d-3)
    integer jflag, npart, i, j, jj
    logical pass
!
    pass = .true.
    jflag = 0
    npart = nexternal - 1
    do i = 0, 3
      xsum(i) = 0.d0
      xsuma(i) = 0.d0
      do j = nincoming + 1, npart
        xsum(i) = xsum(i) + xmom(i, j)
        xsuma(i) = xsuma(i) + abs(xmom(i, j))
      end do
      if (i .eq. 0) xsum(i) = xsum(i) - ecm
      if (xsuma(i) .lt. 1.d0) then
        xrat(i) = abs(xsum(i))
      else
        xrat(i) = abs(xsum(i))/xsuma(i)
      end if
      if (xrat(i) .gt. tiny .and. jflag .eq. 0) then
        write (*, *) 'Momentum is not conserved'
        write (*, *) 'i=', i
        do j = 1, npart
          write (*, '(4(d14.8,1x))') (xmom(jj, j), jj=0, 3)
        end do
        jflag = 1
      end if
    end do
    if (jflag .eq. 1) then
      write (*, '(4(d14.8,1x))') (xsum(jj), jj=0, 3)
      write (*, '(4(d14.8,1x))') (xrat(jj), jj=0, 3)
      pass = .false.
      return
    end if
!
    do j = 1, npart
      do i = 0, 3
        ptmp(i) = xmom(i, j)
      end do
      xm = xlen4(ptmp)
      if (abs(xm - xmass(j))/ptmp(0) .gt. tiny .and. abs(xm - xmass(j)) .gt. tiny) then
        write (*, *) 'Mass shell violation'
        write (*, *) 'j=', j
        write (*, *) 'mass=', xmass(j)
        write (*, *) 'mass computed=', xm
        write (*, '(4(d14.8,1x))') (xmom(jj, j), jj=0, 3)
        pass = .false.
        return
      end if
    end do
    return
  end subroutine phspncheck_born

  subroutine phspncheck_nocms(npart, ecm, xmass, xmom, pass)
! Checks four-momentum conservation. Derived from phspncheck;
! works in any frame
    implicit none
    integer npart
    double precision ecm, xmass(-max_branch:max_particles), xmom(0:3, nexternal)
    double precision tiny, vtiny, xm, den, ecmtmp, xsum(0:3), xsuma(0:3), xrat(0:3), ptmp(0:3)
    parameter(tiny=5.d-3)
    parameter(vtiny=1.d-6)
    integer jflag, i, j, jj
    logical pass
!
    pass = .true.
    jflag = 0
    do i = 0, 3
      if (nincoming .eq. 2) then
        xsum(i) = -xmom(i, 1) - xmom(i, 2)
        xsuma(i) = abs(xmom(i, 1)) + abs(xmom(i, 2))
      elseif (nincoming .eq. 1) then
        xsum(i) = -xmom(i, 1)
        xsuma(i) = abs(xmom(i, 1))
      end if
      do j = nincoming + 1, npart
        xsum(i) = xsum(i) + xmom(i, j)
        xsuma(i) = xsuma(i) + abs(xmom(i, j))
      end do
      if (xsuma(i) .lt. 1.d0) then
        xrat(i) = abs(xsum(i))
      else
        xrat(i) = abs(xsum(i))/xsuma(i)
      end if
      if (xrat(i) .gt. tiny .and. jflag .eq. 0) then
        write (*, *) 'Momentum is not conserved (/nocms/)'
        write (*, *) 'i=', i
        do j = 1, npart
          write (*, '(i2,1x,4(d14.8,1x))') j, (xmom(jj, j), jj=0, 3)
        end do
        jflag = 1
      end if
    end do
    if (jflag .eq. 1) then
      write (*, '(a3,1x,4(d14.8,1x))') 'sum', (xsum(jj), jj=0, 3)
      write (*, '(a3,1x,4(d14.8,1x))') 'rat', (xrat(jj), jj=0, 3)
      pass = .false.
      return
    end if
!
    do j = 1, npart
      do i = 0, 3
        ptmp(i) = xmom(i, j)
      end do
      xm = xlen4(ptmp)
      if (ptmp(0) .ge. 1.d0) then
        den = ptmp(0)
      else
        den = 1.d0
      end if
      if (abs(xm - xmass(j))/den .gt. tiny .and. abs(xm - xmass(j)) .gt. tiny) then
        write (*, *) 'Mass shell violation (/nocms/)'
        write (*, *) 'j=', j
        write (*, *) 'mass=', xmass(j)
        write (*, *) 'mass computed=', xm
        write (*, '(4(d14.8,1x))') (xmom(jj, j), jj=0, 3)
        pass = .false.
        return
      end if
    end do
!
    if (nincoming .eq. 2) then
      ecmtmp = sqrt(2d0*dot(xmom(0, 1), xmom(0, 2)))
    elseif (nincoming .eq. 1) then
      ecmtmp = xmom(0, 1)
    end if
    if (abs(ecm - ecmtmp) .gt. vtiny) then
      write (*, *) 'Inconsistent shat (/nocms/)'
      write (*, *) 'ecm given=   ', ecm
      write (*, *) 'ecm computed=', ecmtmp
      write (*, '(4(d14.8,1x))') (xmom(jj, 1), jj=0, 3)
      write (*, '(4(d14.8,1x))') (xmom(jj, 2), jj=0, 3)
      pass = .false.
      return
    end if

    return
  end subroutine phspncheck_nocms

  function xlen4(v)
    implicit none
    double precision xlen4, tmp, v(0:3)
!
    tmp = v(0)**2 - v(1)**2 - v(2)**2 - v(3)**2
    xlen4 = sign(1.d0, tmp)*sqrt(abs(tmp))
    return
  end function xlen4

  subroutine sreal(pp, xi_i_fks, y_ij_fks, wgt)
! Wrapper for the n+1 contribution. Returns the n+1 matrix element
! squared reduced by the FKS damping factor xi**2*(1-y).
! Close to the soft or collinear limits it calls the corresponding
! Born and multiplies with the AP splitting function or eikonal factors.
    implicit none
    double precision pp(0:3, nexternal), wgt
    double precision xi_i_fks, y_ij_fks

    double precision shattmp




    double precision zero, tiny
    parameter(zero=0d0)

    double precision pmass(nexternal)
    pmass = external_masses
    if (softtest .or. colltest) then
      tiny = 1d-12
    else
      tiny = 1d-6
    end if

    if (pp(0, 1) .le. 0.d0) then
! Unphysical kinematics: set matrix elements equal to zero
      wgt = 0.d0
      return
    end if

! Consistency check -- call to set_cms_stuff() must be done prior to
! entering this function
    if (nincoming .eq. 2) then
      shattmp = 2d0*dot(pp(0, 1), pp(0, 2))
    else
      shattmp = pp(0, 1)**2
    end if
    if (abs(shattmp/shat - 1.d0) .gt. 1.d-5) then
      write (*, *) 'Error in sreal: inconsistent shat'
      write (*, *) shattmp, shat
      stop
    end if

    if (1d0 - y_ij_fks .lt. tiny) then
      if (pmass(j_fks) .eq. zero .and. j_fks .le. nincoming) then
        call sborncol_isr(pp, xi_i_fks, y_ij_fks, wgt)
      elseif (pmass(j_fks) .eq. zero .and. j_fks .ge. nincoming + 1) then
        call sborncol_fsr(pp, y_ij_fks, wgt)
      else
        wgt = 0d0
        amp_split(1:amp_split_size) = 0d0
      end if
    elseif (xi_i_fks .lt. tiny) then
      if (need_color_links) then
! has soft singularities
        call sbornsoft(pp, xi_i_fks, y_ij_fks, wgt)
      else
        wgt = 0d0
        amp_split(1:amp_split_size) = 0d0
      end if
    else
      call smatrix_real(pp, wgt)
      wgt = wgt*xi_i_fks**2*(1d0 - y_ij_fks)
      amp_split(1:amp_split_size) = amp_split(1:amp_split_size)*xi_i_fks**2*(1d0 - y_ij_fks)
    end if

    return
  end subroutine sreal

  subroutine sborncol_fsr(p, y_ij_fks, wgt)
    implicit none
    double precision p(0:3, nexternal), wgt
    double precision y_ij_fks
!




    complex(kind=kind(0d0)) xij_aor


    integer i, imother_fks, iord
    double precision t, z, ap(2), E_j_fks, E_i_fks, Q(2), cphi_mother, sphi_mother, pi(0:3), pj(0:3), wgt_born
    complex(kind=kind(0d0)) W1(6), W2(6), W3(6), W4(6), Wij_angle, Wij_recta
    complex(kind=kind(0d0)) azifact

! Colour representations of i_fks, j_fks and the FKS mother

    double precision zero, vtiny
    parameter(zero=0d0)
    parameter(vtiny=1d-8)
    complex(kind=kind(0d0)) ximag
    parameter(ximag=(0.d0, 1.d0))
    double precision amp_split_local(amp_split_size)
    complex(kind=kind(0d0)) wgt1(2)
    complex(kind=kind(0d0)) ans_extra_cnt(2, nsplitorders)

!
    amp_split_local(1:amp_split_size) = 0d0

!
    if (p_born(0, 1) .le. 0.d0) then
! Unphysical kinematics: set matrix elements equal to zero
      write (*, *) "No born momenta in sborncol_fsr"
      wgt = 0.d0
      return
    end if

    E_j_fks = p(0, j_fks)
    E_i_fks = p(0, i_fks)
    z = 1d0 - E_i_fks/(E_i_fks + E_j_fks)
    t = z*shat/4d0
    call sborn(p_born, wgt_born)
    if (iextra_cnt .gt. 0) call extra_cnt(p_born, iextra_cnt, ans_extra_cnt)
    call AP_reduced(j_type, i_type, t, z, ap)
    call Qterms_reduced_timelike(j_type, i_type, t, z, Q)
    wgt = 0d0
    iord = qcd_pos
    if (.not. split_type(iord)) then
      amp_split(1:amp_split_size) = 0d0
      return
    end if
! check if any extra_cnt is needed
    if (iextra_cnt .gt. 0) then
      if (iord .eq. isplitorder_born) then
        call sborn(p_born, wgt_born)
        wgt1(1) = ans_cnt(1, iord)
        wgt1(2) = ans_cnt(2, iord)
      elseif (iord .eq. isplitorder_cnt) then
! this is the contribution from the extra cnt
        call extra_cnt(p_born, iextra_cnt, ans_extra_cnt)
        wgt1(1) = ans_extra_cnt(1, iord)
        wgt1(2) = ans_extra_cnt(2, iord)
      else
        write (*, *) 'ERROR in sborncol_fsr', iord
        stop
      end if
    else
      call sborn(p_born, wgt_born)
      wgt1(1) = ans_cnt(1, iord)
      wgt1(2) = ans_cnt(2, iord)
    end if
    if (abs(j_type) .eq. 3 .and. i_type .eq. 8) then
      Q(1) = 0d0
      wgt1(2) = 0d0
    elseif (m_type .eq. 8) then
! Insert <ij>/(/ij/) which is not included by sborn()
      if (1d0 - y_ij_fks .lt. vtiny) then
        azifact = xij_aor
      else
        do i = 0, 3
          pi(i) = p_i_fks_ev(i)
          pj(i) = p(i, j_fks)
        end do
        call IXXXSO(pi, ZERO, +1, +1, W1)
        call OXXXSO(pj, ZERO, -1, +1, W2)
        call IXXXSO(pi, ZERO, -1, +1, W3)
        call OXXXSO(pj, ZERO, +1, +1, W4)
        Wij_angle = (0d0, 0d0)
        Wij_recta = (0d0, 0d0)
        do i = 1, 4
          Wij_angle = Wij_angle + W1(i)*W2(i)
          Wij_recta = Wij_recta + W3(i)*W4(i)
        end do
        azifact = Wij_angle/Wij_recta
      end if
! Insert the extra factor due to Madgraph convention for polarization vectors
      imother_fks = min(i_fks, j_fks)
      call getaziangles(p_born(:, imother_fks), cphi_mother, sphi_mother)
      wgt1(2) = -(cphi_mother - ximag*sphi_mother)**2*wgt1(2)*azifact
      amp_split_cnt(1:amp_split_size, 2, iord) = &
        -(cphi_mother - ximag*sphi_mother)**2 &
        *amp_split_cnt(1:amp_split_size, 2, iord)*azifact
    else
      write (*, *) 'FATAL ERROR in sborncol_fsr', i_type, j_type, i_fks, j_fks
      stop 1
    end if
    wgt = wgt + dble(wgt1(1)*ap(1) + wgt1(2)*Q(1))
    amp_split_local(1:amp_split_size) = &
      amp_split_local(1:amp_split_size) + &
      dble(amp_split_cnt(1:amp_split_size, 1, iord)*AP(1) + &
           amp_split_cnt(1:amp_split_size, 2, iord)*Q(1))
    wgt = wgt*iden_comp
    amp_split(1:amp_split_size) = amp_split_local(1:amp_split_size)*iden_comp
    return
  end subroutine sborncol_fsr

  subroutine sborncol_isr(p, xi_i_fks, y_ij_fks, wgt)
    implicit none
    double precision p(0:3, nexternal), wgt
    double precision xi_i_fks, y_ij_fks
!

    double precision p_born_used(0:3, nexternal - 1)




    complex(kind=kind(0d0)) xij_aor

! Colour representations of i_fks, j_fks and the FKS mother

    integer i, iord
    double precision t, z, ap(2), Q(2), cphi_mother, sphi_mother, pi(0:3), pj(0:3), wgt_born
    complex(kind=kind(0d0)) W1(6), W2(6), W3(6), W4(6), Wij_angle, Wij_recta
    complex(kind=kind(0d0)) azifact

    double precision zero, vtiny
    parameter(zero=0d0)
    parameter(vtiny=1d-8)
    complex(kind=kind(0d0)) ximag
    parameter(ximag=(0.d0, 1.d0))
    double precision amp_split_local(amp_split_size)
    complex(kind=kind(0d0)) amp_split_cnt_local(amp_split_size, 2, nsplitorders)
    complex(kind=kind(0d0)) wgt1(2)
    complex(kind=kind(0d0)) ans_extra_cnt(2, nsplitorders)


!
    amp_split_local(1:amp_split_size) = 0d0

! in the case of the collinear CT, use p_born_coll
!  (when not doing event projection).
! For the soft-collinear one, use p_born
    if (xi_i_fks .gt. 0d0 .and. .not. use_evpr) then
      p_born_used(:, :) = p_born_coll(:, :)
    else ! if (xi_i_fks.eq.0d0) then
      p_born_used(:, :) = p_born(:, :)
    end if

    if (p_born_used(0, 1) .le. 0.d0) then
! Unphysical kinematics: set matrix elements equal to zero
      write (*, *) "No born momenta in sborncol_isr"
      wgt = 0.d0
      return
    end if

    z = 1d0 - xi_i_fks
! sreal return {\cal M} of FKS except for the partonic flux 1/(2*s).
! Thus, an extra factor z (implicit in the flux of the reduced Born
! in FKS) has to be inserted here
    t = z*shat/4d0
    call AP_reduced(m_type, i_type, t, z, ap)
    call Qterms_reduced_spacelike(m_type, i_type, t, z, Q)
    wgt = 0d0
    iord = qcd_pos
    if (.not. split_type(iord)) then
      amp_split(1:amp_split_size) = 0d0
      return
    end if
! check if any extra_cnt is needed
    if (iextra_cnt .gt. 0) then
      if (iord .eq. isplitorder_born) then
! this is the contribution from the born ME
        call sborn(p_born_used, wgt_born)
        wgt1(1:2) = ans_cnt(1:2, iord)
      else if (iord .eq. isplitorder_cnt) then
! this is the contribution from the extra cnt
        call extra_cnt(p_born_used, iextra_cnt, ans_extra_cnt)
        wgt1(1:2) = ans_extra_cnt(1:2, iord)
      else
        write (*, *) 'ERROR in sborncol_isr', iord
        stop
      end if
    else
      call sborn(p_born_used, wgt_born)
      wgt1(1:2) = ans_cnt(1:2, iord)
    end if
    amp_split_cnt_local(1:amp_split_size, 1, iord) = amp_split_cnt(1:amp_split_size, 1, iord)
    amp_split_cnt_local(1:amp_split_size, 2, iord) = amp_split_cnt(1:amp_split_size, 2, iord)
    if (abs(m_type) .eq. 3) then
      Q(1) = 0d0
      wgt1(2) = cmplx(0d0, 0d0, kind=kind(0d0))
      amp_split_cnt_local(1:amp_split_size, 2, iord) = cmplx(0d0, 0d0, kind=kind(0d0))
    else
! Insert <ij>/(/ij/) which is not included by sborn()
      if (1d0 - y_ij_fks .lt. vtiny) then
        azifact = xij_aor
      else
        do i = 0, 3
          pi(i) = p_i_fks_ev(i)
          pj(i) = p(i, j_fks)
        end do
        if (j_fks .eq. 2 .and. nincoming .eq. 2) then
! Rotation according to innerpin.m. Use rotate_invar() if a more
! general rotation is needed
          pi(1) = -pi(1)
          pi(3) = -pi(3)
          pj(1) = -pj(1)
          pj(3) = -pj(3)
        end if
        call IXXXSO(pi, ZERO, +1, +1, W1)
        call OXXXSO(pj, ZERO, -1, +1, W2)
        call IXXXSO(pi, ZERO, -1, +1, W3)
        call OXXXSO(pj, ZERO, +1, +1, W4)
        Wij_angle = (0d0, 0d0)
        Wij_recta = (0d0, 0d0)
        do i = 1, 4
          Wij_angle = Wij_angle + W1(i)*W2(i)
          Wij_recta = Wij_recta + W3(i)*W4(i)
        end do
        azifact = Wij_angle/Wij_recta
      end if
! Insert the extra factor due to Madgraph convention for polarization vectors
      cphi_mother = 1.d0
      sphi_mother = 0.d0
      wgt1(2) = -(cphi_mother + ximag*sphi_mother)**2*wgt1(2)*conjg(azifact)
      amp_split_cnt_local(1:amp_split_size, 2, iord) = &
        -(cphi_mother + ximag*sphi_mother)**2 &
        *amp_split_cnt_local(1:amp_split_size, 2, iord)*conjg(azifact)
    end if
    wgt = wgt + dble(wgt1(1)*ap(1) + wgt1(2)*Q(1))
    amp_split_local(1:amp_split_size) = &
      amp_split_local(1:amp_split_size) + &
      dble(amp_split_cnt_local(1:amp_split_size, 1, iord)*AP(1) + &
           amp_split_cnt_local(1:amp_split_size, 2, iord)*Q(1))
    wgt = wgt*iden_comp
    amp_split(1:amp_split_size) = amp_split_local(1:amp_split_size)*iden_comp
    return
  end subroutine sborncol_isr

  subroutine xkplus(PDFscheme, col1, col2, x, xkk)
! This function returns the quantity K^{(+)}_{ab}(x), relevant for
! the MS --> DIS (or any other) scheme change in the factorization scheme.
! It also includes regular terms, multiplied by (1-x).
! There's NO multiplicative (1-x) factor like in the previous functions.
    implicit none
    integer PDFscheme, col1, col2
    double precision x, xkk(2)

    double precision vcf, vtf
    parameter(vcf=4.d0/3.d0)
    parameter(vtf=1.d0/2.d0)

!
    xkk(2) = 0d0
    if (PDFscheme .eq. 0) then
! MSbar, all terms are zero
      xkk(:) = 0d0
    else if (PDFscheme .eq. 1) then
! DIS scheme
      if (col1 .eq. 8 .and. col2 .eq. 8) then ! gg
        xkk(1) = -2*nf*vtf*(1 - x)*(-(x**2 + (1 - x)**2)*log(x) + 8*x*(1 - x) - 1)
        xkk(2) = 0d0
      elseif (abs(col1) .eq. 3 .and. abs(col2) .eq. 3) then ! qq
        xkk(1) = vtf*(1 - x)*(-(x**2 + (1 - x)**2)*log(x) + 8*x*(1 - x) - 1)
      elseif (col1 .eq. 8 .and. abs(col2) .eq. 3) then ! gq
        xkk(1) = -vcf*(-3.d0/2.d0 - (1 + x**2)*log(x) + (1 - x)*(3 + 2*x))
      elseif (abs(col1) .eq. 3 .and. col2 .eq. 8) then ! qg
        xkk(1) = vcf*(-3.d0/2.d0 - (1 + x**2)*log(x) + (1 - x)*(3 + 2*x))
      else
        write (6, *) 'Error in xkplus: wrong colour values', col1, col2
        stop
      end if
    else
      write (6, *) 'Error in xkplus: wrong PDF scheme', PDFscheme
      stop
    end if
    xkk(1) = xkk(1)*g**2
    return
  end subroutine xkplus

  subroutine xklog(PDFscheme, col1, col2, x, xkk)
! This function returns the quantity K^{(l)}_{ab}(x), relevant for
! the MS --> DIS (or any other) scheme change in the factorization scheme.
! There's NO multiplicative (1-x) factor like in the previous functions.
    implicit none
    integer PDFscheme, col1, col2
    double precision x, xkk(2)

    double precision vcf, vtf
    parameter(vcf=4.d0/3.d0)
    parameter(vtf=1.d0/2.d0)
!
    xkk(2) = 0d0
    if (PDFscheme .eq. 0) then
! MSbar, all terms are zero
      xkk(:) = 0d0
    else if (PDFscheme .eq. 1) then
! DIS scheme
      if (col1 .eq. 8 .and. col2 .eq. 8) then ! gg
        xkk(1) = -2*nf*vtf*(1 - x)*(x**2 + (1 - x)**2)
        xkk(2) = 0d0
      elseif (abs(col1) .eq. 3 .and. abs(col2) .eq. 3) then ! qq
        xkk(1) = vtf*(1 - x)*(x**2 + (1 - x)**2)
      elseif (col1 .eq. 8 .and. abs(col2) .eq. 3) then ! gq
        xkk(1) = -vcf*(1 + x**2)
      elseif (abs(col1) .eq. 3 .and. col2 .eq. 8) then ! qg
        xkk(1) = vcf*(1 + x**2)
      else
        write (6, *) 'Error in xklog: wrong colour values', col1, col2
        stop
      end if
    else
      write (6, *) 'Error in xklog: wrong PDF scheme', PDFscheme
      stop
    end if
    xkk(1) = xkk(1)*g**2
    return
  end subroutine xklog

  subroutine xkdelta(PDFscheme, col1, col2, xkk)
! This function returns the quantity K^{(d)}_{ab}, relevant for
! the MS --> DIS (or any other) scheme change in the factorization scheme.
    implicit none
    integer PDFscheme, col1, col2
    double precision xkk(2)

    double precision pi, vcf
    parameter(pi=3.14159265358979312d0)
    parameter(vcf=4.d0/3.d0)
!
    xkk(2) = 0d0
    if (PDFscheme .eq. 0) then
! MSbar, all terms are zero
      xkk(:) = 0d0
    else if (PDFscheme .eq. 1) then
! DIS scheme
      if (col1 .eq. 8 .and. col2 .eq. 8) then ! gg
        xkk(1) = 0.d0
        xkk(2) = 0.d0
      elseif (abs(col1) .eq. 3 .and. abs(col2) .eq. 3) then ! qq
        xkk(1) = 0.d0
      elseif (col1 .eq. 8 .and. abs(col2) .eq. 3) then ! gq
        xkk(1) = vcf*(9.d0/2.d0 + pi**2/3.d0)
      elseif (abs(col1) .eq. 3 .and. col2 .eq. 8) then ! qg
        xkk(1) = -vcf*(9.d0/2.d0 + pi**2/3.d0)
      else
        write (6, *) 'Error in xkdelta: wrong colour values', col1, col2
        stop
      end if
    else
      write (6, *) 'Error in xkdelta: wrong PDF scheme', PDFscheme
      stop
    end if
    xkk(1) = xkk(1)*g**2
    return
  end subroutine xkdelta

  subroutine AP_reduced(col1, col2, t, z, ap)
! Returns Altarelli-Parisi splitting function summed/averaged over helicities
! times prefactors such that |M_n+1|^2 = ap * |M_n|^2. This means
!    AP_reduced = (1-z) P_{S(part1,part2)->part1+part2}(z) * g^2/t
! Therefore, the labeling conventions for particle IDs are not as in FKS:
! part1 and part2 are the two particles emerging from the branching.
! part1 and part2 can be either gluon (8) or (anti-)quark (+-3). z is the
! fraction of the energy of part1 and t is the invariant mass of the mother.
    implicit none

    integer col1, col2
    double precision z, ap(2), t

    double precision CA, TR, CF
    parameter(CA=3d0, TR=1d0/2d0, CF=4d0/3d0)
    ap(2) = 0d0

    if (col1 .eq. 8 .and. col2 .eq. 8) then
! g->gg splitting
      ap(1) = 2d0*CA*((1d0 - z)**2/z + z + z*(1d0 - z)**2)
      ap(2) = 0d0

    elseif (abs(col1) .eq. 3 .and. abs(col2) .eq. 3) then
! g->qqbar splitting
      ap(1) = TR*(z**2 + (1d0 - z)**2)*(1d0 - z)

    elseif (abs(col1) .eq. 3 .and. col2 .eq. 8) then
! q->qg splitting
      ap(1) = CF*(1d0 + z**2)

    elseif (col1 .eq. 8 .and. abs(col2) .eq. 3) then
! q->gq splitting
      ap(1) = CF*(1d0 + (1d0 - z)**2)*(1d0 - z)/z

    else
      write (*, *) 'Fatal Error in AP_reduced', col1, col2
      stop
    end if

    ap(1) = ap(1)*g**2/t
    return
  end subroutine AP_reduced

  subroutine AP_reduced_prime(col1, col2, t, z, apprime)
! Returns (1-z)*P^\prime * gS^2/t, with the same conventions as AP_reduced
    implicit none

    integer col1, col2
    double precision z, apprime(2), t

    double precision TR, CF
    parameter(TR=1d0/2d0, CF=4d0/3d0)
    apprime(2) = 0d0
    if (col1 .eq. 8 .and. col2 .eq. 8) then
! g->gg splitting
      apprime(1) = 0d0
      apprime(2) = 0d0

    elseif (abs(col1) .eq. 3 .and. abs(col2) .eq. 3) then
! g->qqbar splitting
      apprime(1) = -2*TR*z*(1d0 - z)**2

    elseif (abs(col1) .eq. 3 .and. col2 .eq. 8) then
! q->qg splitting
      apprime(1) = -CF*(1d0 - z)**2

    elseif (col1 .eq. 8 .and. abs(col2) .eq. 3) then
! q->gq splitting
      apprime(1) = -CF*z*(1d0 - z)
    else
      write (*, *) 'Fatal error in AP_reduced_prime', col1, col2
      stop
    end if

    apprime(1) = apprime(1)*g**2/t
    return
  end subroutine AP_reduced_prime

  subroutine Qterms_reduced_timelike(col1, col2, t, z, Qterms)
! Eq's B.31 to B.34 of FKS paper, times (1-z)*g^2/t. The labeling
! conventions for particle IDs are the same as those in AP_reduced
    implicit none

    integer col1, col2
    double precision z, Qterms(2), t

    double precision CA, TR
    parameter(CA=3d0, TR=1d0/2d0)
    Qterms(2) = 0d0
    if (col1 .eq. 8 .and. col2 .eq. 8) then
! g->gg splitting
      Qterms(1) = -4d0*CA*z*(1d0 - z)**2
      Qterms(2) = 0d0

    elseif (abs(col1) .eq. 3 .and. abs(col2) .eq. 3) then
! g->qqbar splitting
      Qterms(1) = 4d0*TR*z*(1d0 - z)**2

    elseif (abs(col1) .eq. 3 .and. col2 .eq. 8) then
! q->qg splitting
      Qterms(1) = 0d0

    elseif (col1 .eq. 8 .and. abs(col2) .eq. 3) then
! q->gq splitting
      Qterms(1) = 0d0
    else
      write (*, *) 'Fatal error in Qterms_reduced_timelike', col1, col2
      stop
    end if

    Qterms(1) = Qterms(1)*g**2/t
    return
  end subroutine Qterms_reduced_timelike

  subroutine Qterms_reduced_spacelike(col1, col2, t, z, Qterms)
! Eq's B.42 to B.45 of FKS paper, times (1-z)*gS^2/t. The labeling
! conventions for particle IDs are the same as those in AP_reduced.
! Thus, part1 has momentum fraction z, and it is the one off-shell
! (see (FKS.B.41))
    implicit none

    integer col1, col2
    double precision z, Qterms(2), t

    double precision CA, CF
    parameter(CA=3d0, CF=4d0/3d0)
    Qterms(2) = 0d0
    if (col1 .eq. 8 .and. col2 .eq. 8) then
! g->gg splitting
      Qterms(1) = -4d0*CA*(1d0 - z)**2/z
      Qterms(2) = 0d0

    elseif (abs(col1) .eq. 3 .and. abs(col2) .eq. 3) then
! g->qqbar splitting
      Qterms(1) = 0d0

    elseif (abs(col1) .eq. 3 .and. col2 .eq. 8) then
! q->qg splitting
      Qterms(1) = 0d0

    elseif (col1 .eq. 8 .and. abs(col2) .eq. 3) then
! q->gq splitting
      Qterms(1) = -4d0*CF*(1d0 - z)**2/z
    else
      write (*, *) 'Fatal error in Qterms_reduced_spacelike', col1, col2
      stop
    end if

    Qterms(1) = Qterms(1)*g**2/t
    return
  end subroutine Qterms_reduced_spacelike

  subroutine sbornsoft(pp, xi_i_fks, y_ij_fks, wgt)
    implicit none
!      include "fks.inc"
    integer m, n

    double precision softcontr, pp(0:3, nexternal), wgt, eik, xi_i_fks, y_ij_fks
    double precision wgt1
    integer i, j


    double precision zero, pmass(nexternal)
    parameter(zero=0d0)


    pmass = external_masses
!
! Call the Born to be sure that 'CalculatedBorn' is done correctly. This
! should always be done before calling the color-correlated Borns,
! because of the caching of the diagrams.
!
    call sborn(p_born, wgt1)
!
! Reset the amp_split array
    amp_split(1:amp_split_size) = 0d0

    softcontr = 0d0
    do i = 1, fks_j_from_i(i_fks, 0)
      do j = 1, i
        m = fks_j_from_i(i_fks, i)
        n = fks_j_from_i(i_fks, j)
        if ((m .ne. n .or. (m .eq. n .and. pmass(m) .ne. ZERO)) .and. n .ne. i_fks .and. m .ne. i_fks) then
! wgt includes the gs/w^2
          call sborn_sf(p_born, m, n, wgt)
          if (wgt .ne. 0d0) then
            call eikonal_reduced(pp, m, n, i_fks, j_fks, xi_i_fks, y_ij_fks, eik)
            softcontr = softcontr + wgt*eik*iden_comp
! update the amp_split array
            amp_split(1:amp_split_size) = amp_split(1:amp_split_size) - 2d0*eik*amp_split_soft(1:amp_split_size)*iden_comp
          end if
        end if
      end do
    end do
    wgt = softcontr
! Add minus sign to compensate the minus in the color factor
! of the color-linked Borns (b_sf_0??.f)
! Factor two to fix the limits.
    wgt = -2d0*wgt
    return
  end subroutine sbornsoft

  subroutine eikonal_reduced(pp, m, n, i_fks, j_fks, xi_i_fks, y_ij_fks, eik)
!     Returns the eikonal factor
    implicit none
    double precision eik, pp(0:3, nexternal), xi_i_fks, y_ij_fks
    double precision dotnm, dotni, dotmi, fact
    integer n, m, i_fks, j_fks, i
    integer softcol


    double precision phat_i_fks(0:3)

    double precision zero, pmass(nexternal), tiny
    parameter(zero=0d0)
    parameter(tiny=1d-6)
    pmass = external_masses
! Define the reduced momentum for i_fks
    softcol = 0
    if (1d0 - y_ij_fks .lt. tiny) softcol = 2
    if (p_i_fks_cnt(0, softcol) .lt. 0d0) then
      if (xi_i_fks .eq. 0.d0) then
        write (*, *) 'Error #1 in eikonal_reduced', softcol, xi_i_fks, y_ij_fks
        stop
      end if
      if (pp(0, i_fks) .ne. 0.d0) then
        write (*, *) 'WARNING in eikonal_reduced: no cnt momenta', softcol, xi_i_fks, y_ij_fks
        do i = 0, 3
          phat_i_fks(i) = pp(i, i_fks)/xi_i_fks
        end do
      else
        write (*, *) 'Error #2 in eikonal_reduced', softcol, xi_i_fks, y_ij_fks
        stop
      end if
    else
      do i = 0, 3
        phat_i_fks(i) = p_i_fks_cnt(i, softcol)
      end do
    end if
! Calculate the eikonal factor
    dotnm = dot(pp(0, n), pp(0, m))
    if ((m .ne. j_fks .and. n .ne. j_fks) .or. pmass(j_fks) .ne. ZERO) then
      dotmi = dot(pp(0, m), phat_i_fks)
      dotni = dot(pp(0, n), phat_i_fks)
      fact = 1d0 - y_ij_fks
    elseif (m .eq. j_fks .and. n .ne. j_fks .and. pmass(j_fks) .eq. ZERO) then
      dotni = dot(pp(0, n), phat_i_fks)
      dotmi = sqrtshat/2d0*pp(0, j_fks)
      fact = 1d0
    elseif (m .ne. j_fks .and. n .eq. j_fks .and. pmass(j_fks) .eq. ZERO) then
      dotni = sqrtshat/2d0*pp(0, j_fks)
      dotmi = dot(pp(0, m), phat_i_fks)
      fact = 1d0
    else
      write (*, *) 'Error #3 in eikonal_reduced'
      stop
    end if

    eik = dotnm/(dotni*dotmi)*fact
    return
  end subroutine eikonal_reduced

  subroutine sreal_deg(p, xi_i_fks, collrem_xi, collrem_lxi)
    use extra_weights
    implicit none
    integer iord, iap
    double precision p(0:3, nexternal), collrem_xi, collrem_lxi
    double precision xi_i_fks
    double precision collrem_xi_tmp, collrem_lxi_tmp

    double precision wgt_born

    double precision p_born_used(0:3, nexternal - 1)




    double precision shattmp, oo2pi, z, t, ap(2), apprime(2), xkkernp(2), xkkernd(2), xkkernl(2), xnorm

! Colour representations of i_fks, j_fks and the FKS mother
    complex(kind=kind(0d0)) wgt1(2)

    double precision one, pi
    parameter(one=1.d0)
    parameter(pi=3.1415926535897932385d0)

    complex(kind=kind(0d0)) ans_extra_cnt(2, nsplitorders)

    double precision amp_split_collrem_xi(amp_split_size), amp_split_collrem_lxi(amp_split_size)
! amp_split for the PDF scheme
    double precision prefact_xi


    logical firsttime_pdf
    data firsttime_pdf/.true./

! The fixed-order template supports the MSbar (0) and DIS (1) schemes.
    if (firsttime_pdf) then
      write (*, *) 'PDFscheme', pdfscheme
      firsttime_pdf = .false.
    end if

    amp_split_collrem_xi(1:amp_split_size) = 0d0
    amp_split_collrem_lxi(1:amp_split_size) = 0d0
    amp_split_wgtdegrem_xi(1:amp_split_size) = 0d0
    amp_split_wgtdegrem_lxi(1:amp_split_size) = 0d0
    amp_split_wgtdegrem_muF(1:amp_split_size) = 0d0
    amp_split_wgtpsch_p(1:amp_split_size) = 0d0
    amp_split_wgtpsch_l(1:amp_split_size) = 0d0
    amp_split_wgtpsch_d(1:amp_split_size) = 0d0

! in the case of the collinear CT, use p_born_coll
!  (when not doing event projection).
! For the soft-collinear one, use p_born
    if (xi_i_fks .gt. 0d0 .and. .not. use_evpr) then
      p_born_used(:, :) = p_born_coll(:, :)
    else ! if (xi_i_fks.eq.0d0) then
      p_born_used(:, :) = p_born(:, :)
    end if

    if (j_fks .gt. nincoming) then
! Do not include this contribution for final-state branchings
      collrem_xi = 0.d0
      collrem_lxi = 0.d0
      if (doreweight) then
        wgtdegrem_xi = 0.d0
        wgtdegrem_lxi = 0.d0
        wgtdegrem_muF = 0.d0
      end if
      return
    end if

    if (p_born_used(0, 1) .le. 0.d0) then
! Unphysical kinematics: set matrix elements equal to zero
      write (*, *) "No born momenta in sreal_deg"
      collrem_xi = 0.d0
      collrem_lxi = 0.d0
      if (doreweight) then
        wgtdegrem_xi = 0.d0
        wgtdegrem_lxi = 0.d0
        wgtdegrem_muF = 0.d0
      end if
      return
    end if

! Consistency check -- call to set_cms_stuff() must be done prior to
! entering this function
    if (nincoming .eq. 2) then
      shattmp = 2d0*dot(p(0, 1), p(0, 2))
    else
      shattmp = p(0, 1)**2
    end if
    if (abs(shattmp/shat - 1.d0) .gt. 1.d-5) then
      write (*, *) 'Error in sreal: inconsistent shat'
      write (*, *) shattmp, shat
      stop
    end if

! A factor gS^2 is included in the Altarelli-Parisi kernels
    oo2pi = one/(8d0*PI**2)

    z = 1d0 - xi_i_fks
    t = one
    call AP_reduced(m_type, i_type, t, z, ap)
    call AP_reduced_prime(m_type, i_type, t, z, apprime)

! call the PDF-scheme kernels here
!   p-> (/1/(1-z)/)_+
!   l-> (/log(1-z)/(1-z)/)_+
!   d-> delta(1-z)
    call xkplus(PDFscheme, m_type, i_type, z, xkkernp)
    call xkdelta(PDFscheme, m_type, i_type, xkkernd)
    call xklog(PDFscheme, m_type, i_type, z, xkkernl)

    collrem_xi = 0.d0
    collrem_lxi = 0.d0
    calculatedborn = .false.
    iord = qcd_pos
    iap = 1
    if (.not. split_type(iord)) return

! check if any extra_cnt is needed
    if (iextra_cnt .gt. 0) then
      if (iord .eq. isplitorder_born) then
! this is the contribution from the born ME
        call sborn(p_born_used, wgt_born)
        wgt1(1) = ans_cnt(1, iord)
        wgt1(2) = ans_cnt(2, iord)
      else if (iord .eq. isplitorder_cnt) then
! this is the contribution from the extra cnt
        call extra_cnt(p_born_used, iextra_cnt, ans_extra_cnt)
        wgt1(1) = ans_extra_cnt(1, iord)
        wgt1(2) = ans_extra_cnt(2, iord)
      else
        write (*, *) 'ERROR in sreal_deg', iord
        stop
      end if
    else
      call sborn(p_born_used, wgt_born)
      wgt1(1) = ans_cnt(1, iord)
      wgt1(2) = ans_cnt(2, iord)
    end if

    collrem_xi_tmp = ap(iap)*log(shat*delta_used/(2*q2fact(j_fks))) - apprime(iap)
    collrem_lxi_tmp = 2*ap(iap)

! The partonic flux 1/(2*s) is inserted in genps. Thus, an extra
! factor z (implicit in the flux of the reduced Born in FKS)
! has to be inserted here
    xnorm = 1.d0/z*iden_comp

    collrem_xi = collrem_xi + oo2pi*dble(wgt1(1))*collrem_xi_tmp*xnorm
    collrem_lxi = collrem_lxi + oo2pi*dble(wgt1(1))*collrem_lxi_tmp*xnorm

    amp_split_collrem_xi(1:amp_split_size) = &
      amp_split_collrem_xi(1:amp_split_size) + &
      dble(amp_split_cnt(1:amp_split_size, 1, iord))*oo2pi &
      *collrem_xi_tmp*xnorm
    amp_split_collrem_lxi(1:amp_split_size) = &
      amp_split_collrem_lxi(1:amp_split_size) + &
      dble(amp_split_cnt(1:amp_split_size, 1, iord))*oo2pi &
      *collrem_lxi_tmp*xnorm

    prefact_xi = ap(iap)*log(shat*delta_used/(2*QES2)) - apprime(iap)
    amp_split_wgtdegrem_xi(1:amp_split_size) = &
      amp_split_wgtdegrem_xi(1:amp_split_size) + &
      oo2pi*dble(amp_split_cnt(1:amp_split_size, 1, iord)) &
      *prefact_xi*xnorm
    amp_split_wgtdegrem_lxi(1:amp_split_size) = amp_split_collrem_lxi(1:amp_split_size)
    amp_split_wgtdegrem_muF(1:amp_split_size) = &
      amp_split_wgtdegrem_muF(1:amp_split_size) - &
      oo2pi*dble(amp_split_cnt(1:amp_split_size, 1, iord))*ap(iap)*xnorm
! amp split for the PDF scheme
    if (PDFscheme .ne. 0) then
      amp_split_wgtpsch_p(1:amp_split_size) = &
        amp_split_wgtpsch_p(1:amp_split_size) - &
        dble(amp_split_cnt(1:amp_split_size, 1, iord)) &
        *xkkernp(iap)*oo2pi*xnorm
      amp_split_wgtpsch_l(1:amp_split_size) = &
        amp_split_wgtpsch_l(1:amp_split_size) - &
        dble(amp_split_cnt(1:amp_split_size, 1, iord)) &
        *xkkernl(iap)*oo2pi*xnorm
      amp_split_wgtpsch_d(1:amp_split_size) = &
        amp_split_wgtpsch_d(1:amp_split_size) - &
        dble(amp_split_cnt(1:amp_split_size, 1, iord)) &
        *xkkernd(iap)*oo2pi*xnorm
    end if
    calculatedborn = .false.

    return
  end subroutine sreal_deg

  subroutine set_cms_stuff(icountevts)
    implicit none
    integer icountevts







! rapidity of boost from \tilde{k}_1+\tilde{k}_2 c.m. frame to lab frame --
! same for event and counterevents
! This is the rapidity that enters in the arguments of the sinh() and
! cosh() of the boost, in such a way that
!       y(k)_lab = y(k)_tilde - ybst_til_tolab
! where y(k)_lab and y(k)_tilde are the rapidities computed with a generic
! four-momentum k, in the lab frame and in the \tilde{k}_1+\tilde{k}_2
! c.m. frame respectively
    ybst_til_tolab = -ycm_cnt(0) - 0.5d0*log(ebeam(1)/ebeam(2))
    if (icountevts .eq. -100) then
! set Bjorken x's in run.inc for the computation of PDFs in auto_dsig
      xbk(1) = xbjrk_ev(1)
      xbk(2) = xbjrk_ev(2)
! shat=2*k1.k2 -- consistency of this assignment with momenta checked
! in phspncheck_nocms
      shat = shat_ev
      sqrtshat = sqrtshat_ev
! rapidity of boost from \tilde{k}_1+\tilde{k}_2 c.m. frame to
! k_1+k_2 c.m. frame
      ybst_til_tocm = ycm_ev - ycm_cnt(0)
    else
! do the same as above for the counterevents
      xbk(1) = xbjrk_cnt(1, icountevts)
      xbk(2) = xbjrk_cnt(2, icountevts)
      shat = shat_cnt(icountevts)
      sqrtshat = sqrtshat_cnt(icountevts)
      ybst_til_tocm = ycm_cnt(icountevts) - ycm_cnt(0)
    end if
    return
  end subroutine set_cms_stuff

  subroutine xmom_compare(i_fks, j_fks, jac_cnt, p1_cnt, pass)
    implicit none
    double precision p1_cnt(0:3, nexternal, 0:2), jac_cnt(0:2)
    integer i_fks, j_fks
    integer izero, ione, itwo, iunit, isum
    logical verbose, pass, pass0
    parameter(izero=0)
    parameter(ione=1)
    parameter(itwo=2)
    parameter(iunit=6)
    parameter(verbose=.false.)
!
    isum = 0
    if (jac_cnt(0) .gt. 0.d0) isum = isum + 1
    if (jac_cnt(1) .gt. 0.d0) isum = isum + 2
    if (jac_cnt(2) .gt. 0.d0) isum = isum + 4
    pass = .true.
!
    if (isum .eq. 0 .or. isum .eq. 1 .or. isum .eq. 2 .or. isum .eq. 4) then
! Nothing to be done: 0 or 1 configurations computed
      if (verbose) write (iunit, *) 'none'
    elseif (isum .eq. 3 .or. isum .eq. 5 .or. isum .eq. 7) then
! Soft is taken as reference
      if (isum .eq. 7) then
        if (verbose) then
          write (iunit, *) 'all'
          write (iunit, *) '    '
          write (iunit, *) 'C/S'
        end if
        call xmcompare(verbose, pass0, ione, izero, i_fks, j_fks, p1_cnt)
        pass = pass .and. pass0
        if (verbose) then
          write (iunit, *) '    '
          write (iunit, *) 'SC/S'
        end if
        call xmcompare(verbose, pass0, itwo, izero, i_fks, j_fks, p1_cnt)
        pass = pass .and. pass0
      elseif (isum .eq. 3) then
        if (verbose) then
          write (iunit, *) 'C+S'
          write (iunit, *) '    '
          write (iunit, *) 'C/S'
        end if
        call xmcompare(verbose, pass0, ione, izero, i_fks, j_fks, p1_cnt)
        pass = pass .and. pass0
      elseif (isum .eq. 5) then
        if (verbose) then
          write (iunit, *) 'SC+S'
          write (iunit, *) '    '
          write (iunit, *) 'SC/S'
        end if
        call xmcompare(verbose, pass0, itwo, izero, i_fks, j_fks, p1_cnt)
        pass = pass .and. pass0
      end if
    elseif (isum .eq. 6) then
! Collinear is taken as reference
      if (verbose) then
        write (iunit, *) 'SC+C'
        write (iunit, *) '    '
        write (iunit, *) 'SC/C'
      end if
      call xmcompare(verbose, pass0, itwo, ione, i_fks, j_fks, p1_cnt)
      pass = pass .and. pass0
    else
      write (6, *) 'Fatal error in xmom_compare', isum
      stop
    end if
    if (.not. pass) i_momcmp_count = i_momcmp_count + 1
!
    return
  end subroutine xmom_compare

  subroutine xmcompare(verbose, pass0, inum, iden, i_fks, j_fks, p1_cnt)
    implicit none
    double precision p1_cnt(0:3, nexternal, 0:2)
    logical verbose, pass0
    integer inum, iden, i_fks, j_fks, iunit, ipart, i, j, k
    double precision tiny, vtiny, xnum, xden, xrat
    parameter(iunit=6)
    parameter(tiny=1.d-4)
    parameter(vtiny=1.d-10)
    double precision pmass(nexternal)
    pmass = external_masses
!
    pass0 = .true.
    do ipart = 1, nexternal
      do i = 0, 3
        xnum = p1_cnt(i, ipart, inum)
        xden = p1_cnt(i, ipart, iden)
        if (verbose) then
          if (i .eq. 0) then
            write (iunit, *) ' '
            write (iunit, *) 'part=', ipart
          end if
          call xprintout(iunit, xnum, xden)
        else
          if (ipart .ne. i_fks .and. ipart .ne. j_fks) then
            if (xden .ne. 0.d0) then
              xrat = abs(1 - xnum/xden)
            else
              xrat = abs(xnum)
            end if
            if (abs(xnum) .eq. 0d0 .and. abs(xden) .le. vtiny) xrat = 0d0
! The following line solves some problem as well, but before putting
! it as the standard, one should think a bit about it
            if (abs(xnum) .le. vtiny .and. abs(xden) .le. vtiny) xrat = 0d0
            if (xrat .gt. tiny .and. (pmass(ipart) .eq. 0d0 .or. xnum/pmass(ipart) .gt. vtiny)) then
              write (*, *) 'Kinematics of counterevents'
              write (*, *) inum, iden
              write (*, *) 'is different. Particle:', ipart
              write (*, *) xrat, xnum, xden
              do j = 1, nexternal
                write (*, *) j, (p1_cnt(k, j, inum), k=0, 3)
              end do
              do j = 1, nexternal
                write (*, *) j, (p1_cnt(k, j, iden), k=0, 3)
              end do
              xratmax = max(xratmax, xrat)
              pass0 = .false.
            end if
          end if
        end if
      end do
    end do
    do i = 0, 3
      if (j_fks .gt. nincoming) then
        xnum = p1_cnt(i, i_fks, inum) + p1_cnt(i, j_fks, inum)
        xden = p1_cnt(i, i_fks, iden) + p1_cnt(i, j_fks, iden)
      else
        xnum = -p1_cnt(i, i_fks, inum) + p1_cnt(i, j_fks, inum)
        xden = -p1_cnt(i, i_fks, iden) + p1_cnt(i, j_fks, iden)
      end if
      if (verbose) then
        if (i .eq. 0) then
          write (iunit, *) ' '
          write (iunit, *) 'part=i+j'
        end if
        call xprintout(iunit, xnum, xden)
      else
        if (xden .ne. 0.d0) then
          xrat = abs(1 - xnum/xden)
        else
          xrat = abs(xnum)
        end if
        if (xrat .gt. tiny) then
          write (*, *) 'Kinematics of counterevents'
          write (*, *) inum, iden
          write (*, *) 'is different. Particle i+j'
          xratmax = max(xratmax, xrat)
          pass0 = .false.
        end if
      end if
    end do
    return
  end subroutine xmcompare

  subroutine xprintout(iunit, xv, xlim)
    implicit none
    integer iunit
    double precision xv, xlim
!
    if (abs(xlim) .gt. 1.d-30) then
      write (iunit, *) xv/xlim, xv, xlim
    else
      write (iunit, *) xv, xlim
    end if
    return
  end subroutine xprintout

! The following has been derived with minor modifications from the
! analogous routine written for VBF

! The following has been derived with minor modifications from the
! analogous routine written for VBF
  subroutine checkres2(xsecvc, xseclvc, wgt, wgtl, xp, lxp, iflag, imax, iev, i_fks, j_fks, iret)
!     same as checkres, but also limits are arrays.
    implicit none
    double precision xsecvc(15), xseclvc(15), wgt(15), wgtl(15), lxp(0:3, nexternal + 1), xp(15, 0:3, nexternal + 1)
    double precision ckc(15), rckc(15), rat
    integer iflag, imax, iev, i_fks, j_fks, iret, ithrs, istop, iwrite, i, k, l, imin, icount
    parameter(ithrs=3)
    parameter(istop=0)
    parameter(iwrite=1)
!
    if (imax .gt. 15) then
      write (6, *) 'Error in checkres: imax is too large', imax
      stop
    end if
    do i = 1, imax
      if (xseclvc(i) .eq. 0.d0) then
        ckc(i) = abs(xsecvc(i))
      else
        ckc(i) = abs(xsecvc(i)/xseclvc(i) - 1.d0)
      end if
    end do
    if (iflag .eq. 0) then
      rat = 4.d0
    elseif (iflag .eq. 1) then
      rat = 2.d0
    else
      write (6, *) 'Error in checkres: iflag=', iflag
      write (6, *) ' Must be 0 for soft, 1 for collinear'
      stop
    end if
!
    i = 1
    do while (ckc(i) .gt. 0.1d0 .and. xseclvc(i) .ne. 0d0)
      i = i + 1
    end do
    imin = i
    do i = imin, imax - 1
      if (ckc(i + 1) .ne. 0.d0) then
        rckc(i) = ckc(i)/ckc(i + 1)
      else
        rckc(i) = 1.d8
      end if
    end do
    icount = 0
    i = imin
    do while (icount .lt. ithrs .and. i .lt. imax)
      if (rckc(i) .gt. rat) then
        icount = icount + 1
      else
        icount = 0
      end if
      i = i + 1
    end do
!
    iret = 0
    if (icount .ne. ithrs) then
      iret = 1
      if (istop .eq. 1) then
        write (6, *) 'Test failed', iflag
        write (6, *) 'Event #', iev
        stop
      end if
      if (iwrite .eq. 1) then
        write (77, *) '    '
        if (iflag .eq. 0) then
          write (77, *) 'Soft #', iev
        elseif (iflag .eq. 1) then
          write (77, *) 'Collinear #', iev
        end if
        write (77, *) 'ME*wgt:'
        do i = 1, imax
          call xprintout(77, xsecvc(i), xseclvc(i))
        end do
        write (77, *) 'wgt:'
        do i = 1, imax
          call xprintout(77, wgt(i), wgtl(i))
        end do
!
        write (78, *) '    '
        if (iflag .eq. 0) then
          write (78, *) 'Soft #', iev
        elseif (iflag .eq. 1) then
          write (78, *) 'Collinear #', iev
        end if
        do k = 1, nexternal
          write (78, *) ''
          write (78, *) 'part:', k
          do l = 0, 3
            write (78, *) 'comp:', l
            do i = 1, imax
              call xprintout(78, xp(i, l, k), lxp(l, k))
            end do
          end do
        end do
        if (iflag .eq. 0) then
          write (78, *) ''
          write (78, *) 'part: i_fks reduced'
          do l = 0, 3
            write (78, *) 'comp:', l
            do i = 1, imax
              call xprintout(78, xp(i, l, nexternal + 1), lxp(l, nexternal + 1))
            end do
          end do
          write (78, *) ''
          write (78, *) 'part: i_fks full/reduced'
          do l = 0, 3
            write (78, *) 'comp:', l
            do i = 1, imax
              call xprintout(78, xp(i, l, i_fks), xp(i, l, nexternal + 1))
            end do
          end do
        elseif (iflag .eq. 1) then
          write (78, *) ''
          write (78, *) 'part: i_fks+j_fks'
          do l = 0, 3
            write (78, *) 'comp:', l
            do i = 1, imax
              call xprintout(78, xp(i, l, i_fks) + xp(i, l, j_fks), lxp(l, i_fks) + lxp(l, j_fks))
            end do
          end do
        end if
      end if
    end if
    return
  end subroutine checkres2

  subroutine bornsoftvirtual(p, bsv_wgt, virt_wgt, born_wgt)
    use extra_weights
    use mint_module
    implicit none
    real :: tBefore, tAfter
!      include "fks.inc"
    double precision p(0:3, nexternal), bsv_wgt, born_wgt, avv_wgt
    double precision wgt1
    double precision Q, Ej, wgt, contr, eikIreg
    double precision aso2pi
    double precision shattmp
    integer i, j, aj, m, n, k



    double precision pi
    parameter(pi=3.1415926535897932385d0)

    double precision c_used, gamma_used, gammap_used
    double precision double, single, xmu2
    logical ComputePoles, fksprefact
    parameter(ComputePoles=.false.)
    parameter(fksprefact=.true.)


    double precision virt_wgt



! timing statistics
! For the MINT folding
    double precision virt_wgt_save
    save virt_wgt_save

    double precision pmass(nexternal), zero
    parameter(zero=0d0)
    logical firsttime
    data firsttime/.true./
    logical need_color_links_used
    data need_color_links_used/.false./
    double precision oneo8pi2
    parameter(oneo8pi2=1d0/(8d0*pi**2))
    integer nFKSprocess_save, nFKSprocess_col
    data nFKSprocess_col/0/
    double precision bsv_wgt_mufoqes, bsv_wgt_mufomur
    double precision contr_mufoqes, contr_mufomur
! to keep track of the various split orders
    integer iamp
    integer orders(nsplitorders)
    double precision amp_split_born(amp_split_size)
    double precision amp_split_bsv(amp_split_size)
    pmass = external_masses
    if (firsttime) then
! Check whether any real-emission configuration needs colour links.
      nFKSprocess_save = nFKSprocess
      do nFKSprocess = 1, FKS_configs
        call fks_inc_chooser()
        need_color_links_used = need_color_links_used .or. need_color_links
! Keep track of a configuration that needs colour links.
        if (need_color_links .and. nFKSprocess_col .eq. 0) nFKSprocess_col = nFKSprocess
      end do
      if (need_color_links_used) then
        write (*, *) 'Color-linked born are used'
      else
        write (*, *) 'Color-linked born are not used'
      end if
      firsttime = .false.
      nFKSprocess = nFKSprocess_save
      call fks_inc_chooser()
    end if

    aso2pi = g**2/(8*pi**2)

    amp_split_bsv(1:amp_split_size) = 0d0
    amp_split_virt(1:amp_split_size) = 0d0
    amp_split_avv(1:amp_split_size) = 0d0

    if (.not. need_color_links_used) then
! just return 0
      bsv_wgt = 0d0
      virt_wgt = 0d0
      born_wgt = 0d0
      goto 999
    end if

! Consistency check -- call to set_cms_stuff() must be done prior to
! entering this function
    if (nincoming .eq. 2) then
      shattmp = 2d0*dot(p(0, 1), p(0, 2))
    else
      shattmp = p(0, 1)**2
    end if
    if (abs(shattmp/shat - 1.d0) .gt. 1.d-5) then
      write (*, *) 'Error in bornsoftvirtual: inconsistent shat'
      write (*, *) shattmp, shat
      stop
    end if

    call sborn(p_born, wgt1)

! Born contribution:
    bsv_wgt = wgt1
    born_wgt = wgt1
    virt_wgt = 0d0
    avv_wgt = 0d0
    amp_split_born(1:amp_split_size) = amp_split(1:amp_split_size)
    amp_split_bsv(1:amp_split_size) = amp_split(1:amp_split_size)

    if (abrv .eq. 'born') goto 549
    if (abrv .eq. 'virt') goto 547

! Q contribution eq 5.5 and 5.6 of FKS
    Q = 0d0
    if (split_type_used(qcd_pos)) then
      do i = 1, nexternal
        if (i .ne. i_fks .and. pmass(i) .eq. ZERO) then
! set the colour factors according to the
! type of the leg
          if (particle_type(i) .eq. 8) then
            aj = 0
          elseif (abs(particle_type(i)) .eq. 3) then
            aj = 1
          else
            aj = -1
          end if
          Ej = p(0, i)

          if (aj .eq. -1) cycle
          c_used = c(aj)
          gamma_used = gamma(aj)
          gammap_used = gammap(aj)

          if (i .gt. nincoming) then
! Q terms for final state partons
            if (abrv .ne. 'virt') then
! 1+2+3+4
              Q = Q + gammap_used &
                  - dlog(shat*deltaO/2d0/QES2) &
                  *(gamma_used &
                    - 2d0*c_used*dlog(2d0*Ej/xicut_used/sqrtshat)) &
                  + 2d0*c_used*(dlog(2d0*Ej/sqrtshat)**2 &
                                - dlog(xicut_used)**2) &
                  - 2d0*gamma_used*dlog(2d0*Ej/sqrtshat)
            else
              write (*, *) 'Error in bornsoftvirtual'
              write (*, *) 'abrv in Q:', abrv
              stop
            end if

          else
! Q terms for initial state partons
            if (abrv .ne. 'virt') then
! 1+2+3+4
              Q = Q - dlog(q2fact(i)/QES2)*(gamma_used + 2d0*c_used*dlog(xicut_used))
            else
              write (*, *) 'Error in bornsoftvirtual'
              write (*, *) 'abrv in Q:', abrv
              stop
            end if
          end if
        end if
      end do
! end of the external particle loop
      bsv_wgt = bsv_wgt + aso2pi*Q*dble(ans_cnt(1, qcd_pos))
      amp_split_bsv(1:amp_split_size) = amp_split_bsv(1:amp_split_size) + aso2pi*Q*dble(amp_split_cnt(1:amp_split_size, 1, qcd_pos))
    end if

!     If doing MC over helicities, must sum over the two
!     helicity contributions for the Q-terms of collinear limit.
547 continue
    if (abrv .eq. 'virt') goto 548
!
! I(reg) terms, eq 5.5 of FKS
    nFKSprocess_save = nFKSprocess
    if (need_color_links_used) then
      need_color_links = need_color_links_used
      nFKSprocess = nFKSprocess_col
! setup the fks i/j info
      call fks_inc_chooser()
! the following call to born is to setup the goodhel(nfksprocess)
      call sborn(p_born, wgt1)
      contr = 0d0
      do i = 1, fks_j_from_i(i_fks, 0)
        do j = 1, i
          m = fks_j_from_i(i_fks, i)
          n = fks_j_from_i(i_fks, j)
          if ((m .ne. n .or. (m .eq. n .and. pmass(m) .ne. ZERO)) .and. n .ne. i_fks .and. m .ne. i_fks) then
! To be sure that color-correlated Borns work well, we need to have
! *always* a call to sborn(p_born,wgt) just before. This is okay,
! because there is a call above in this subroutine
! wgt includes the gs/w^2
            call sborn_sf(p_born, m, n, wgt)
            if (wgt .ne. 0d0) then
              call eikonal_Ireg(p, m, n, xicut_used, eikIreg)
              contr = contr + wgt*eikIreg
              do k = 1, amp_split_size
                amp_split_bsv(k) = amp_split_bsv(k) - 2d0*eikIreg*oneo8pi2*amp_split_soft(k)
              end do
            end if
          end if
        end do
      end do

! WARNING: THE FACTOR -2 BELOW COMPENSATES FOR THE MISSING -2 IN THE
! COLOUR LINKED BORN -- SEE ALSO SBORNSOFT().
! If the colour-linked Borns were normalized as reported in the paper
! we should set
!   bsv_wgt=bsv_wgt+ao2pi*contr  <-- DO NOT USE THIS LINE
!
      bsv_wgt = bsv_wgt - 2*oneo8pi2*contr
    end if

! set back the fks i/j info as prior to enter this function
    nFKSprocess = nFKSprocess_save
    call fks_inc_chooser()

548 continue
! Finite part of one-loop corrections
! convert to Binoth Les Houches Accord standards
    virt_wgt = 0d0

    call sborn(p_born, wgt1)
! Use the QCD counterterm Born to approximate the virtual.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     THIS IS DANGEROUS: if these are not always the same for all
!     events, the whole virt_trics doesn't work and gives wrong results!
!     CHECK THIS.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
    do iamp = 1, amp_split_size
      amp_split_virt(iamp) = 0d0
      amp_split_born_for_virt(iamp) = 0d0
      if (dble(amp_split_cnt(iamp, 1, qcd_pos)) .ne. 0d0) then
        amp_split_born_for_virt(iamp) = dble(amp_split_cnt(iamp, 1, qcd_pos))
      end if
    end do

    if (fold .eq. 0) then
      if ((ran2() .le. virtual_fraction(ichan) .and. abrv(1:3) .ne. 'nov') .or. abrv(1:4) .eq. 'virt') then
        call cpu_time(tBefore)
        call BinothLHA(p_born, born_wgt, virt_wgt)
        do iamp = 1, amp_split_size
          amp_split_virt(iamp) = amp_split_finite_ML(iamp)
        end do
        virtual_over_born = virt_wgt/born_wgt
        virt_wgt = 0d0
        do iamp = 1, amp_split_size
          if (amp_split_virt(iamp) .eq. 0d0) cycle
          if (use_poly_virtual) then
            amp_split_virt(iamp) = amp_split_virt(iamp) - polyfit(iamp)*amp_split_born_for_virt(iamp)
          else
            amp_split_virt(iamp) = amp_split_virt(iamp) - average_virtual(iamp, ichan)*amp_split_born_for_virt(iamp)
          end if
          virt_wgt = virt_wgt + amp_split_virt(iamp)
        end do
        if (abrv .ne. 'virt') then
          virt_wgt = virt_wgt/virtual_fraction(ichan)
          do iamp = 1, amp_split_size
            amp_split_virt(iamp) = amp_split_virt(iamp)/virtual_fraction(ichan)
          end do
        end if
        call cpu_time(tAfter)
        tOLP = tOLP + (tAfter - tBefore)
      end if
      virt_wgt_save = virt_wgt
      amp_split_virt_save(1:amp_split_size) = amp_split_virt(1:amp_split_size)
    elseif (fold .eq. 1) then
      virt_wgt = virt_wgt_save
      amp_split_virt(1:amp_split_size) = amp_split_virt_save(1:amp_split_size)
    end if
    if (abrv(1:4) .ne. 'virt') then
      if (use_poly_virtual) then
        avv_wgt = polyfit(0)*born_wgt
        do iamp = 1, amp_split_size
          if (amp_split_born_for_virt(iamp) .eq. 0d0) cycle
          amp_split_avv(iamp) = polyfit(iamp)*amp_split_born_for_virt(iamp)
        end do
      else
        avv_wgt = average_virtual(0, ichan)*born_wgt
        do iamp = 1, amp_split_size
          if (amp_split_born_for_virt(iamp) .eq. 0d0) cycle
          amp_split_avv(iamp) = average_virtual(iamp, ichan)*amp_split_born_for_virt(iamp)
        end do
      end if
    end if

! eq.(MadFKS.C.13)
    if (abrv .ne. 'virt') then
! this is to update the amp_split array
      call sborn(p_born, wgt1)
      bsv_wgt_mufoqes = 0d0
      do iamp = 1, amp_split_size
        if (dble(amp_split_cnt(iamp, 1, qcd_pos)) .eq. 0d0) cycle
        call amp_split_pos_to_orders(iamp, orders)
        contr_mufoqes = pi*beta0*dble(orders(qcd_pos) - 2) &
                        *log(q2fact(1)/QES2)*aso2pi*dble(amp_split_cnt(iamp, 1, qcd_pos))
        amp_split_bsv(iamp) = amp_split_bsv(iamp) + contr_mufoqes
        bsv_wgt_mufoqes = bsv_wgt_mufoqes + contr_mufoqes
      end do
      bsv_wgt = bsv_wgt + bsv_wgt_mufoqes
    end if

!  eq.(MadFKS.C.14)
    if (abrv(1:2) .ne. 'vi') then
      bsv_wgt_mufomur = 0d0
      do iamp = 1, amp_split_size
        if (dble(amp_split_cnt(iamp, 1, qcd_pos)) .eq. 0d0) cycle
        call amp_split_pos_to_orders(iamp, orders)
        contr_mufomur = -pi*beta0*dble(orders(qcd_pos) - 2) &
                        *log(q2fact(1)/scale**2)*aso2pi*dble(amp_split_cnt(iamp, 1, qcd_pos))
        amp_split_bsv(iamp) = amp_split_bsv(iamp) + contr_mufomur
        bsv_wgt_mufomur = bsv_wgt_mufomur + contr_mufomur
      end do
      bsv_wgt = bsv_wgt + bsv_wgt_mufomur
    end if

549 continue

    wgtwnstmpmuf = 0.d0
    wgtnstmp = 0d0
    wgtwnstmpmur = 0.d0
    amp_split_wgtnstmp(1:amp_split_size) = 0d0
    amp_split_wgtwnstmpmuf(1:amp_split_size) = 0d0
    amp_split_wgtwnstmpmur(1:amp_split_size) = 0d0

    if (abrv .ne. 'born' .and. abrv .ne. 'grid') then
      call sborn(p_born, wgt1)
      if (abrv(1:2) .eq. 'vi') then
        wgtwnstmpmur = 0.d0
      else
        if (split_type_used(qcd_pos)) then
          do i = 1, nincoming
            if (particle_type(i) .eq. 8) then
              aj = 0
            elseif (abs(particle_type(i)) .eq. 3) then
              aj = 1
            else
              aj = -1
            end if
            if (aj .eq. -1) cycle
            c_used = c(aj)
            gamma_used = gamma(aj)
            gammap_used = gammap(aj)
            do iamp = 1, amp_split_size
              if (dble(amp_split_cnt(iamp, 1, qcd_pos)) .eq. 0d0) cycle
              amp_split_wgtwnstmpmuf(iamp) = &
                amp_split_wgtwnstmpmuf(iamp) &
                - (gamma_used + 2d0*c_used*dlog(xicut_used)) &
                *dble(amp_split_cnt(iamp, 1, qcd_pos))*aso2pi
            end do
          end do            !end loop i=1,nincoming
        end if
        do iamp = 1, amp_split_size
          if (dble(amp_split_cnt(iamp, 1, qcd_pos)) .eq. 0d0) cycle
          call amp_split_pos_to_orders(iamp, orders)
          amp_split_wgtwnstmpmur(iamp) = dble(amp_split_cnt(iamp, 1, qcd_pos)) &
                                         *pi*beta0*dble(orders(qcd_pos) - 2)*aso2pi
        end do
      end if
! bsv_wgt here always contains the Born; must subtract it, since
! we need the pure NLO terms only
      amp_split_wgtnstmp(1:amp_split_size) = &
        amp_split_bsv(1:amp_split_size) - amp_split_born(1:amp_split_size) &
        - log(q2fact(1)/QES2) &
        *amp_split_wgtwnstmpmuf(1:amp_split_size) &
        - log(scale**2/QES2) &
        *amp_split_wgtwnstmpmur(1:amp_split_size)
    end if

    amp_split(1:amp_split_size) = amp_split_bsv(1:amp_split_size)

    if (abrv(1:2) .eq. 'vi') then
      bsv_wgt = bsv_wgt - born_wgt
      born_wgt = 0d0
    end if

    if (ComputePoles) then
      call sborn(p_born, wgt1)

      print *, "           "
      write (*, 123) ((p(i, j), i=0, 3), j=1, nexternal)
      xmu2 = q2fact(1)
      call getpoles(p, xmu2, double, single, fksprefact)
      print *, "BORN", born_wgt!/conv
      print *, "DOUBLE", double
      print *, "SINGLE", single
!         print*,"LOOP",virt_wgt!/born_wgt/ao2pi*2d0
!         print*,"LOOP2",(virtcor+born_wgt*4d0/3d0-double*pi**2/6d0)
!         stop
123   format(4(1x, d22.16))
    end if

999 continue
    return
  end subroutine bornsoftvirtual

  subroutine eikonal_Ireg(p, m, n, xicut_used, eikIreg)
    implicit none
    double precision pi, pi2
    parameter(pi=3.1415926535897932385d0)
    parameter(pi2=pi**2)
    double precision p(0:3, nexternal), xicut_used, eikIreg
    integer m, n



    double precision Ei, Ej, kikj, rij, tmp, xmj, betaj, betai
    double precision xmi2, xmj2, vij, xi0, alij, tHVvl, tHVv
    double precision arg1, arg2, arg3, arg4, xi1a

    double precision pmass(nexternal)
    pmass = external_masses
    tmp = 0.d0
    if (pmass(m) .eq. 0.d0 .and. pmass(n) .eq. 0.d0) then
      if (m .eq. n) then
        write (*, *) 'Error #2 in eikonal_Ireg', m, n
        stop
      end if
      Ei = p(0, n)
      Ej = p(0, m)
      kikj = dot(p(0, n), p(0, m))
      rij = kikj/(2*Ei*Ej)
      if (abs(rij - 1.d0) .gt. 1.d-6) then
        if (abrv .ne. 'virt') then
! 1+2+3+4
          tmp = 1d0/2d0*dlog(xicut_used**2*shat/QES2)**2 &
                + dlog(xicut_used**2*shat/QES2)*dlog(rij) &
                - ddilog(rij) + 1d0/2d0*dlog(rij)**2 &
                - dlog(1 - rij)*dlog(rij)
        else
          write (*, *) 'Error #11 in eikonal_Ireg', abrv
          stop
        end if
      else
        if (abrv .ne. 'virt') then
! 1+2+3+4
          tmp = 1d0/2d0*dlog(xicut_used**2*shat/QES2)**2 - pi2/6.d0
        else
          write (*, *) 'Error #12 in eikonal_Ireg', abrv
          stop
        end if
      end if
    elseif ((pmass(m) .ne. 0.d0 .and. pmass(n) .eq. 0.d0) .or. (pmass(m) .eq. 0.d0 .and. pmass(n) .ne. 0.d0)) then
      if (m .eq. n) then
        write (*, *) 'Error #3 in eikonal_Ireg', m, n
        stop
      end if
      if (pmass(m) .ne. 0.d0 .and. pmass(n) .eq. 0.d0) then
        Ei = p(0, n)
        Ej = p(0, m)
        xmj = pmass(m)
        betaj = sqrt(1 - xmj**2/Ej**2)
      else
        Ei = p(0, m)
        Ej = p(0, n)
        xmj = pmass(n)
        betaj = sqrt(1 - xmj**2/Ej**2)
      end if
      kikj = dot(p(0, n), p(0, m))
      rij = kikj/(2*Ei*Ej)

      if (abrv .ne. 'virt') then
! 1+2+3+4
        tmp = dlog(xicut_used) &
              *(dlog(xicut_used*shat/QES2) + 2*dlog(kikj/(xmj*Ei))) &
              - ddilog(1 - (1 + betaj)/(2*rij)) &
              + ddilog(1 - 2*rij/(1 - betaj)) &
              + 1/2.d0*log(2*rij/(1 - betaj))**2 &
              + dlog(shat/QES2)*dlog(kikj/(xmj*Ei)) - pi2/12.d0 &
              + 1/4.d0*dlog(shat/QES2)**2 &
              - 1/4.d0*dlog((1 + betaj)/(1 - betaj))**2
      else
        write (*, *) 'Error #13 in eikonal_Ireg', abrv
        stop
      end if
    elseif (pmass(m) .ne. 0.d0 .and. pmass(n) .ne. 0.d0) then
      if (n .eq. m) then
        Ei = p(0, n)
        betai = sqrt(1 - pmass(n)**2/Ei**2)
        if (abrv .ne. 'virt') then
! 1+2+3+4
          if (betai .gt. 1d-6) then
            tmp = dlog(xicut_used**2*shat/QES2) - 1/betai*dlog((1 + betai)/(1 - betai))
          else
            tmp = dlog(xicut_used**2*shat/QES2) - 2d0*(1d0 + betai**2/3d0 + betai**4/5d0)
          end if
        else
          write (*, *) 'Error #14 in eikonal_Ireg', abrv
          stop
        end if
      else
        Ei = p(0, n)
        Ej = p(0, m)
        betai = sqrt(1 - pmass(n)**2/Ei**2)
        betaj = sqrt(1 - pmass(m)**2/Ej**2)
        xmi2 = pmass(n)**2
        xmj2 = pmass(m)**2
        kikj = dot(p(0, n), p(0, m))
        vij = sqrt(1 - xmi2*xmj2/kikj**2)
        alij = kikj*(1 + vij)/xmi2
        tHVvl = (alij**2*xmi2 - xmj2)/2.d0
        tHVv = tHVvl/(alij*Ei - Ej)
        arg1 = alij*Ei
        arg2 = arg1*betai
        arg3 = Ej
        arg4 = arg3*betaj
        if (vij .lt. 1d0) then
          xi0 = 1/vij*log((1 + vij)/(1 - vij))
        else
          xi0 = dlog(4d0*kikj**2/(xmi2*xmj2))
        end if
!          xi0=1/vij*log((1+vij)/(1-vij))
        xi1a = kikj**2*(1 + vij)/xmi2*(xj1a(arg1, arg2, tHVv, tHVvl) - xj1a(arg3, arg4, tHVv, tHVvl))

        if (abrv .ne. 'virt') then
! 1+2+3+4
          tmp = 1/2.d0*xi0*dlog(xicut_used**2*shat/QES2) + 1/2.d0*xi1a
        else
          write (*, *) 'Error #15 in eikonal_Ireg', abrv
          stop
        end if
      end if
    else
      write (*, *) 'Error #4 in eikonal_Ireg', m, n, pmass(m), pmass(n)
      stop
    end if
    eikIreg = tmp
    return
  end subroutine eikonal_Ireg

  function xj1a(x, y, tHVv, tHVvl)
    implicit none
    double precision xj1a, x, y, tHVv, tHVvl
!
    xj1a = 1/(2*tHVvl)*(dlog((x - y)/(x + y))**2 + 4*ddilog(1 - (x + y)/tHVv) + 4*ddilog(1 - (x - y)/tHVv))
    return
  end function xj1a

  function DDILOG(X)
!
! $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
!
! $Log: imp64.inc,v $
! Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
! Mathlib gen
!
!
! imp64.inc
!
    implicit none
    integer i
    double precision ddilog, x, y, s, a, t, h, alfa, b0, b1, b2, c(0:19)
    double precision z1, hf, pi, pi3, pi6, pi12
    parameter(Z1=1, HF=Z1/2)
    parameter(PI=3.14159265358979324d0)
    parameter(PI3=PI**2/3, PI6=PI**2/6, PI12=PI**2/12)
    data C(0)/0.42996693560813697d0/
    data C(1)/0.40975987533077105d0/
    data C(2)/-0.01858843665014592d0/
    data C(3)/0.00145751084062268d0/
    data C(4)/-0.00014304184442340d0/
    data C(5)/0.00001588415541880d0/
    data C(6)/-0.00000190784959387d0/
    data C(7)/0.00000024195180854d0/
    data C(8)/-0.00000003193341274d0/
    data C(9)/0.00000000434545063d0/
    data C(10)/-0.00000000060578480d0/
    data C(11)/0.00000000008612098d0/
    data C(12)/-0.00000000001244332d0/
    data C(13)/0.00000000000182256d0/
    data C(14)/-0.00000000000027007d0/
    data C(15)/0.00000000000004042d0/
    data C(16)/-0.00000000000000610d0/
    data C(17)/0.00000000000000093d0/
    data C(18)/-0.00000000000000014d0/
    data C(19)/+0.00000000000000002d0/
    if (X .eq. 1) then
      H = PI6
    elseif (X .eq. -1) then
      H = -PI12
    else
      T = -X
      if (T .le. -2) then
        Y = -1/(1 + T)
        S = 1
        A = -PI3 + HF*(log(-T)**2 - log(1 + 1/T)**2)
      elseif (T .lt. -1) then
        Y = -1 - T
        S = -1
        A = log(-T)
        A = -PI6 + A*(A + log(1 + 1/T))
      else if (T .le. -HF) then
        Y = -(1 + T)/T
        S = 1
        A = log(-T)
        A = -PI6 + A*(-HF*A + log(1 + T))
      else if (T .lt. 0) then
        Y = -T/(1 + T)
        S = -1
        A = HF*log(1 + T)**2
      else if (T .le. 1) then
        Y = T
        S = 1
        A = 0
      else
        Y = 1/T
        S = -1
        A = PI6 + HF*log(T)**2
      end if
      H = Y + Y - 1
      ALFA = H + H
      B1 = 0
      B2 = 0
      do I = 19, 0, -1
        B0 = C(I) + ALFA*B1 - B2
        B2 = B1
        B1 = B0
      end do
      H = -(S*(B0 - H*B2) + A)
    end if
    DDILOG = H
    return
  end function DDILOG

  subroutine getpoles(p, xmu2, double, single, fksprefact, split_poles)
! Returns the residues of double and single poles according to
! eq.(B.1) and eq.(B.2) if fksprefact=.true.. When fksprefact=.false.,
! the prefactor (mu2/Q2)^ep in eq.(B.1) is expanded, and giving an
! extra contribution to the single pole
    implicit none
!      include "fks.inc"
    double precision p(0:3, nexternal), xmu2, double, single
    logical fksprefact
    double precision, optional, intent(out) :: split_poles(amp_split_size, 2)
    double precision wgt1
    double precision born, wgt, kikj, vij, aso2pi
    double precision contr1, contr2
    integer aj, i, j, m, n, k
    double precision pmass(nexternal), zero, pi
    parameter(pi=3.1415926535897932385d0)
    parameter(zero=0d0)
    double precision oneo8pi2
    parameter(oneo8pi2=1d0/(8d0*pi**2))
    integer nFKSprocess_save, nFKSprocess_col
    logical need_color_links_used
    double precision soft_fact
    pmass = external_masses
    nFKSprocess_col = 0

    need_color_links_used = .false.

! Check whether any real-emission configuration needs colour links.
    nFKSprocess_save = nFKSprocess
    do nFKSprocess = 1, FKS_configs
      call fks_inc_chooser()
      need_color_links_used = need_color_links_used .or. need_color_links
      if (need_color_links .and. nFKSprocess_col .eq. 0) nFKSprocess_col = nFKSprocess
    end do
    nFKSprocess = nFKSprocess_save
    call fks_inc_chooser()

    double = 0.d0
    single = 0.d0
! reset the amp_split_poles_FKS
    do i = 1, amp_split_size
      amp_split_poles_FKS(i, 1) = 0d0
      amp_split_poles_FKS(i, 2) = 0d0
    end do
    aso2pi = g**2/(8d0*pi**2)
    call sborn(p_born, wgt1)
! QCD Born terms
    contr1 = 0d0
    contr2 = 0d0
    born = dble(ans_cnt(1, qcd_pos))
    do i = 1, nexternal
      if (i .ne. i_fks .and. particle_type(i) .ne. 1) then
        if (particle_type(i) .eq. 8) then
          aj = 0
        elseif (abs(particle_type(i)) .eq. 3) then
          aj = 1
        end if
        if (pmass(i) .eq. ZERO) then
          contr2 = contr2 - c(aj)
          contr1 = contr1 - gamma(aj)
        else
          contr1 = contr1 - c(aj)
        end if
      end if
    end do

    double = double + contr2*born*aso2pi
    single = single + contr1*born*aso2pi

    do i = 1, amp_split_size
      amp_split_poles_FKS(i, 1) = amp_split_poles_FKS(i, 1) + dble(amp_split_cnt(i, 1, qcd_pos))*contr1*aso2pi
      amp_split_poles_FKS(i, 2) = amp_split_poles_FKS(i, 2) + dble(amp_split_cnt(i, 1, qcd_pos))*contr2*aso2pi
    end do

! Colour-linked Born terms
    nFKSprocess_save = nFKSprocess
    if (need_color_links_used) then
      need_color_links = .true.
      nFKSprocess = nFKSprocess_col

! setup the fks i/j info
      call fks_inc_chooser()
! the following call to born is to setup the goodhel(nfksprocess)
      call sborn(p_born, wgt1)

      contr1 = 0d0
      do i = 1, fks_j_from_i(i_fks, 0)
        do j = 1, i
          m = fks_j_from_i(i_fks, i)
          n = fks_j_from_i(i_fks, j)
          if (m .ne. n .and. n .ne. i_fks .and. m .ne. i_fks) then
! wgt includes the gs/w^2 factor
            call sborn_sf(p_born, m, n, wgt)
! The factor -2 compensate for that missing in sborn_sf
            wgt = -2d0*wgt
            if (wgt .ne. 0.d0) then
              if (pmass(m) .eq. zero .and. pmass(n) .eq. zero) then
                kikj = dot(p(0, n), p(0, m))
                soft_fact = dlog(2d0*kikj/QES2)
              elseif (pmass(m) .ne. zero .and. pmass(n) .eq. zero) then
                kikj = dot(p(0, n), p(0, m))
                soft_fact = -0.5d0*dlog(pmass(m)**2/QES2) + dlog(2d0*kikj/QES2)
              elseif (pmass(m) .eq. zero .and. pmass(n) .ne. zero) then
                kikj = dot(p(0, n), p(0, m))
                soft_fact = -0.5d0*dlog(pmass(n)**2/QES2) + dlog(2d0*kikj/QES2)
              elseif (pmass(m) .ne. zero .and. pmass(n) .ne. zero) then
                kikj = dot(p(0, n), p(0, m))
                vij = dsqrt(1d0 - (pmass(n)*pmass(m)/kikj)**2)
                if (vij .gt. 1d-6) then
                  soft_fact = 0.5d0*1/vij*log((1 + vij)/(1 - vij))
                else
                  soft_fact = (1d0 + vij**2/3d0 + vij**4/5d0)
                end if
              else
                write (*, *) 'Error in getpoles', i, j, n, m, pmass(n), pmass(m)
                stop
              end if
              contr1 = contr1 + soft_fact*wgt
              do k = 1, amp_split_size
                amp_split_poles_FKS(k, 1) = amp_split_poles_FKS(k, 1) + amp_split_soft(k)*(-2d0)*soft_fact*oneo8pi2
              end do
            end if
          end if
        end do
      end do
      single = single + contr1*oneo8pi2
    end if

! Restore the selected FKS configuration.
    nFKSprocess = nFKSprocess_save
    call fks_inc_chooser()

    if (.not. fksprefact) single = single + double*dlog(xmu2/QES2)
    if (present(split_poles)) split_poles = amp_split_poles_FKS
!
    return
  end subroutine getpoles

  subroutine setfksfactor
    use weight_lines
    use extra_weights
    use mint_module
    implicit none

    double precision CA, CF, PI
    parameter(CA=3d0, CF=4d0/3d0)
    parameter(pi=3.1415926535897932385d0)




    integer i, j, fac1, fac2, kchan, open_status





    double precision dfac1
    integer fac_i, fac_j, i_fks_pdg, j_fks_pdg, iden(nexternal)



! Colour representations of i_fks, j_fks and the FKS mother
    softtest = .false.
    colltest = .false.

    if (j_fks .gt. nincoming) then
      delta_used = deltaO
    else
      delta_used = deltaI
    end if

    xicut_used = xicut
    xiScut_used = xiScut
    if (nbody .or. (abrv .eq. 'born' .or. abrv .eq. 'grid' .or. abrv(1:2) .eq. 'vi')) then
      xiBSVcut_used = 1d0
    else
      xiBSVcut_used = xiBSVcut
    end if

    c(0) = CA
    c(1) = CF
    gamma(0) = (11d0*CA - 2d0*Nf)/6d0
    gamma(1) = CF*3d0/2d0
    gammap(0) = (67d0/9d0 - 2d0*PI**2/3d0)*CA - 23d0/18d0*Nf
    gammap(1) = (13/2d0 - 2d0*PI**2/3d0)*CF

! Beta_0 defined according to (MadFKS.C.5)
    beta0 = gamma(0)/(2*pi)
    if (firsttime_nFKSprocess(nFKSprocess)) then
      firsttime_nFKSprocess(nFKSprocess) = .false.
!---------------------------------------------------------------------
!              Symmetry Factors
!---------------------------------------------------------------------
! fkssymmetryfactor:
! Calculate the FKS symmetry factors to be able to reduce the number
! of directories to (maximum) 4 (neglecting quark flavors):
!     1. i_fks=gluon, j_fks=gluon
!     2. i_fks=gluon, j_fks=quark
!     3. i_fks=gluon, j_fks=anti-quark
!     4. i_fks=quark, j_fks=anti-quark (or vice versa).
! This sets the fkssymmetryfactor (in which the quark flavors are taken
! into account) for the subtracted reals.
!
! fkssymmetryfactorBorn:
! Note that in the Born's included here, the final state identical
! particle factor is set equal to the identical particle factor
! for the real contribution to be able to get the correct limits for the
! subtraction terms and the approximated real contributions.
! However when we want to calculate the Born contributions only, we
! have to correct for this difference. Since we only include the Born
! related to a soft limit (this uniquely defines the Born for a given real)
! the difference is always n!/(n-1)!=n, where n is the number of final state
! gluons in the real contribution.
!
! Furthermore, because we are not integrating all the directories, we also
! have to include a fkssymmetryfactor for the Born contributions. However,
! this factor is not the same as the factor defined above, because in this
! case i_fks is fixed to the extra gluon (which goes soft and defines the
! Born contribution) and should therefore not be taken into account when
! calculating the symmetry factor. Together with the factor n above this
! sets the fkssymmetryfactorBorn equal to the fkssymmetryfactor for the
! subtracted reals.
!
! We set fkssymmetryfactorBorn to zero when i_fks not a gluon
!
      i_fks_pdg = pdg_type(i_fks)
      j_fks_pdg = pdg_type(j_fks)

      fac_i_FKS(nFKSprocess) = 0
      fac_j_FKS(nFKSprocess) = 0
      do i = nincoming + 1, nexternal
        if (i_fks_pdg .eq. pdg_type(i)) fac_i_FKS(nFKSprocess) = fac_i_FKS(nFKSprocess) + 1
        if (j_fks_pdg .eq. pdg_type(i)) fac_j_FKS(nFKSprocess) = fac_j_FKS(nFKSprocess) + 1
      end do
! Overwrite if initial state singularity
      if (j_fks .le. nincoming) fac_j_FKS(nFKSprocess) = 1

! i_fks and j_fks of the same type? -> subtract 1 to avoid double counting
      if (j_fks .gt. nincoming .and. i_fks_pdg .eq. j_fks_pdg) fac_j_FKS(nFKSprocess) = fac_j_FKS(nFKSprocess) - 1

! THESE TESTS WORK ONLY FOR FINAL STATE SINGULARITIES
! MZ the test may be removed sooner or later
      if (j_fks .gt. nincoming) then
        if (i_fks_pdg .eq. j_fks_pdg .and. i_fks_pdg .ne. 21) then
          write (*, *) 'ERROR, if PDG type of i_fks and j_fks '// &
            'are equal, they MUST be gluons', &
            i_fks, j_fks, i_fks_pdg, j_fks_pdg
          stop
        elseif (abs(particle_type(i_fks)) .eq. 3) then
          if (particle_type(i_fks) .ne. -particle_type(j_fks) .or. pdg_type(i_fks) .ne. -pdg_type(j_fks)) then
            write (*, *) 'ERROR, if i_fks is a color triplet,'// &
              ' j_fks must be its anti-particle,'// &
              ' or an initial state gluon.', &
              i_fks, j_fks, particle_type(i_fks), &
              particle_type(j_fks), pdg_type(i_fks), pdg_type(j_fks)
            stop
          end if
        elseif (abs(i_fks_pdg) .ne. 21) then ! if not already above, it must be a gluon
          write (*, *) 'ERROR, i_fks is not a gluon and falls not'// &
            ' in other categories', i_fks, j_fks, i_fks_pdg, j_fks_pdg
          stop
        end if
      end if

      ngluons_FKS(nFKSprocess) = 0
      do i = nincoming + 1, nexternal
        if (pdg_type(i) .eq. 21) ngluons_FKS(nFKSprocess) = ngluons_FKS(nFKSprocess) + 1
      end do

! Set color types of i_fks, j_fks and fks_mother.
      i_type = particle_type(i_fks)
      j_type = particle_type(j_fks)
      call get_mother_colour_impl(i_type, j_type, m_type, i_fks, j_fks)
      i_type_FKS(nFKSprocess) = i_type
      j_type_FKS(nFKSprocess) = j_type
      m_type_FKS(nFKSprocess) = m_type

! Compute the identical particle symmetry factor that is in the
! real-emission matrix elements.
      iden_real_FKS(nFKSprocess) = 1
      do i = 1, nexternal
        iden(i) = 1
      end do
      do i = nincoming + 2, nexternal
        do j = nincoming + 1, i - 1
          if (pdg_type(j) .eq. pdg_type(i)) then
            iden(j) = iden(j) + 1
            iden_real_FKS(nFKSprocess) = iden_real_FKS(nFKSprocess)*iden(j)
            exit
          end if
        end do
      end do
! Compute the identical particle symmetry factor that is in the
! Born matrix elements.
      iden_born_FKS(nFKSprocess) = 1
      call weight_lines_allocated(nexternal, max_contr, max_wgt, max_iproc)
      call set_pdg_impl(0, nFKSprocess, idup)
      do i = 1, nexternal
        iden(i) = 1
      end do
      do i = nincoming + 2, nexternal - 1
        do j = nincoming + 1, i - 1
          if (pdg_uborn(j, 0) .eq. pdg_uborn(i, 0)) then
            iden(j) = iden(j) + 1
            iden_born_FKS(nFKSprocess) = iden_born_FKS(nFKSprocess)*iden(j)
            exit
          end if
        end do
      end do
    end if

    i_type = i_type_FKS(nFKSprocess)
    j_type = j_type_FKS(nFKSprocess)
    m_type = m_type_FKS(nFKSprocess)

! Compensating factor needed in the soft & collinear counterterms for
! the fact that the identical particle symmetry factor in the Born
! matrix elements is not the one that should be used for those terms
! (should be the one in the real instead).
    iden_comp = dble(iden_born_FKS(nFKSprocess))/dble(iden_real_FKS(nFKSprocess))

    fac_i = fac_i_FKS(nFKSprocess)
    fac_j = fac_j_FKS(nFKSprocess)
    ngluons = ngluons_FKS(nFKSprocess)
! Setup the FKS symmetry factors.
    if (nbody .and. pdg_type(i_fks) .eq. 21) then
      fkssymmetryfactor = dble(ngluons)
      fkssymmetryfactorDeg = dble(ngluons)
      fkssymmetryfactorBorn = 1d0
    elseif (pdg_type(i_fks) .eq. -21) then
      fkssymmetryfactor = 1d0
      fkssymmetryfactorDeg = 1d0
      fkssymmetryfactorBorn = 1d0
    else
      fkssymmetryfactor = dble(fac_i*fac_j)
      fkssymmetryfactorDeg = dble(fac_i*fac_j)
      if (pdg_type(i_fks) .eq. 21) then
        fkssymmetryfactorBorn = dble(fac_i*fac_j)
      else
        fkssymmetryfactorBorn = 0d0
      end if
      if (abrv .eq. 'grid') then
        fkssymmetryfactorBorn = 1d0
        fkssymmetryfactor = 0d0
        fkssymmetryfactorDeg = 0d0
      end if
    end if

    if (setfks_firsttime) then
! Check to see if this channel needs to be included in the multi-channeling
      do kchan = 1, nchans
        diagramsymmetryfactor_save(kchan) = 0d0
      end do
      if (multi_channel) then
        open (unit=19, file="symfact.dat", status="old", iostat=open_status)
        if (open_status .eq. 0) then
          i = 0
          do
            i = i + 1
            read (19, *, err=23, end=23) dfac1, fac2
            fac1 = nint(dfac1)
            if (nint(dfac1*10) - fac1*10 .eq. 2) then
              i = i - 1
              cycle
            end if
            do kchan = 1, nchans
              if (i .eq. iconfigs(kchan)) then
                if (config_map(iconfigs(kchan), 0) .ne. fac1) then
                  write (*, *) 'inconsistency in symfact.dat', i, kchan, iconfigs(kchan), config_map(iconfigs(kchan), 0), fac1
                  stop
                end if
                diagramsymmetryfactor_save(kchan) = dble(fac2)
              end if
            end do
          end do
23        continue
          close (19)
        else
          diagramsymmetryfactor_save(1:nchans) = 1d0
        end if
      else                   ! no multi_channel
        do kchan = 1, nchans
          diagramsymmetryfactor_save(kchan) = 1d0
        end do
      end if
      setfks_firsttime = .false.
    end if
    diagramsymmetryfactor = diagramsymmetryfactor_save(ichan)

    return

  end subroutine setfksfactor

  subroutine set_mu_central(ic, dd, c_mu2_r, c_mu2_f)
    use weight_lines
    use extra_weights
    implicit none
    integer ic, dd, i, j
    double precision c_mu2_r, c_mu2_f, muR, muF(2), pp(0:3, nexternal)
    if (dd .eq. 1) then
      c_mu2_r = scales2(2, ic)
      c_mu2_f = scales2(3, ic)
    else
! need to recompute the scales using the momenta
      dynamical_scale_choice = dyn_scale(dd)
      do i = 1, nexternal
        do j = 0, 3
          pp(j, i) = momenta(j, i, ic)
        end do
      end do
      call set_ren_scale(pp, muR)
      c_mu2_r = muR**2
      call set_fac_scale(pp, muF)
      c_mu2_f = muF(1)**2
!     reset the default dynamical_scale_choice
      dynamical_scale_choice = dyn_scale(1)
    end if
    return
  end subroutine set_mu_central

  function ran2()
!     Wrapper for the random numbers; needed for the NLO stuff
    use mint_module
    use ranmar_module, only: ntuple
    implicit none
    double precision ran2, x, a, b
    integer jconfig
    a = 0d0                     ! min allowed value for x
    b = 1d0                     ! max allowed value for x
    jconfig = iconfig           ! integration channel (for off-set)
    call ntuple(x, a, b, jconfig)
    ran2 = x
    return
  end function ran2

  subroutine fill_configurations_common
    implicit none
    integer :: saved_process, configuration
    call require_fks_singular_state()
    config_mass = 0d0
    config_width = 0d0
    config_forest = 0
    config_sprop = 0
    config_tprid = 0
    config_map = 0
    call fill_configurations_born(config_forest(:, :, :, 0), &
                                  config_sprop(:, :, 0), config_tprid(:, :, 0), config_map(:, 0), &
                                  config_mass(:, :, 0), config_width(:, :, 0))
    saved_process = nfksprocess
    do configuration = 1, fks_configs
      nfksprocess = configuration
      call configs_and_props_inc_chooser()
      call fill_configurations_real(config_forest(:, :, :, configuration), &
                                    config_sprop(:, :, configuration), config_tprid(:, :, configuration), &
                                    config_map(:, configuration), config_mass(:, :, configuration), &
                                    config_width(:, :, configuration))
    end do
    nfksprocess = saved_process
  end subroutine fill_configurations_common

  subroutine fill_configurations_born(iforest_in, sprop_in, tprid_in, &
                                      mapconfig_in, pmass_in, pwidth_in)
    implicit none
    integer, intent(inout) :: iforest_in(2, -max_branch:-1, lmaxconfigs)
    integer, intent(inout) :: sprop_in(-max_branch:-1, lmaxconfigs)
    integer, intent(inout) :: tprid_in(-max_branch:-1, lmaxconfigs)
    integer, intent(inout) :: mapconfig_in(0:lmaxconfigs)
    double precision, intent(inout) :: pmass_in(-nexternal:0, lmaxconfigs)
    double precision, intent(inout) :: pwidth_in(-nexternal:0, lmaxconfigs)
    integer :: i, j, k
    call require_fks_singular_state()
    do i = 1, born_lmaxconfigs_used
      do j = -born_max_branch_used, -1
        do k = 1, 2
          iforest_in(k, j, i) = born_forest(k, j, i)
        end do
        sprop_in(j, i) = born_sprop(j, i)
        tprid_in(j, i) = born_tprid(j, i)
        pmass_in(j, i) = born_mass(j, i)
        pwidth_in(j, i) = born_width(j, i)
      end do
      mapconfig_in(i) = born_map(i)
    end do
    mapconfig_in(0) = born_map(0)
  end subroutine fill_configurations_born

  subroutine fill_configurations_real(iforest_in, sprop_in, tprid_in, &
                                      mapconfig_in, pmass_in, pwidth_in)
    implicit none
    integer, intent(inout) :: iforest_in(2, -max_branch:-1, lmaxconfigs)
    integer, intent(inout) :: sprop_in(-max_branch:-1, lmaxconfigs)
    integer, intent(inout) :: tprid_in(-max_branch:-1, lmaxconfigs)
    integer, intent(inout) :: mapconfig_in(0:lmaxconfigs)
    double precision, intent(inout) :: pmass_in(-nexternal:0, lmaxconfigs)
    double precision, intent(inout) :: pwidth_in(-nexternal:0, lmaxconfigs)
    integer :: i, j, k
    call require_fks_singular_state()
    do i = 1, lmaxconfigs
      do j = -max_branch, -1
        do k = 1, 2
          iforest_in(k, j, i) = real_forest(k, j, i)
        end do
        sprop_in(j, i) = real_sprop(j, i)
        tprid_in(j, i) = real_tprid(j, i)
        pmass_in(j, i) = real_mass(j, i)
        pwidth_in(j, i) = real_width(j, i)
      end do
      mapconfig_in(i) = real_map(i)
    end do
    mapconfig_in(0) = real_map(0)
  end subroutine fill_configurations_real

end module fks_singular_module
