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
  use kin_functions_module, only: dot => dot_impl
  use phase_space_kinematics, only: getaziangles
  use fks_qcd_splitting, only: xkplus, xklog, xkdelta, AP_reduced, &
                                AP_reduced_prime, &
                                Qterms_reduced_timelike, &
                                Qterms_reduced_spacelike
  use fks_random_module, only: random_unit_interval
  use fks_soft_kernels, only: eikonal_reduced, eikonal_Ireg
  use fks_sij_module, only: initialize_fks_sij_module, &
                            set_fks_sij_partition_state, fks_sij_impl
  use FKSParams, only: use_poly_virtual
  use chooser_functions_module, only: get_mother_colour_impl, set_pdg_impl
  use madfks_plot_module, only: initplot_impl, outfun_impl
  use fnlo_process_common, only: nfksprocess, i_fks, j_fks, &
                                 soft_counterevent, &
                                 collinear_counterevent, &
                                 soft_collinear_counterevent, real_event, &
                                 event_xi, event_y, event_xi_hat, &
                                 event_fks_momentum, event_xi_max, &
                                 event_xi_norm, event_bjorken_x, &
                                 event_sqrt_shat, event_shat, &
                                 ybst_til_tolab, ybst_til_tocm, &
                                 f_b, f_nb, f_r, f_s, f_c, &
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
                                 isplitorder_cnt, iden_comp, &
                                 i_momcmp_count, xratmax, c, &
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
  double precision, pointer :: stored_event_momenta(:, :, :) => null()
  double precision, pointer :: stored_event_jacobian(:) => null()
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

  interface
    double precision function dlum(bjorken_x)
      implicit none
      double precision, intent(in) :: bjorken_x(2)
    end function dlum
  end interface

  public :: compute_born
  public :: compute_nbody_noborn, compute_real_emission
  public :: compute_soft_counter_term, compute_collinear_counter_term
  public :: compute_soft_collinear_ct_impl
  public :: compute_prefactors_nbody, include_multichannel_enhance
  public :: compute_prefactors_n1body, include_pdf_and_alphas
  public :: reweight_scale
  public :: reweight_pdf, fill_pineappl_weights, get_wgt_nbody
  public :: get_wgt_no_nbody, fill_plots, fill_mint_function
  public :: sreal
  public :: getpoles, setfksfactor, fill_configurations_common
  public :: initialize_fks_model_state, initialize_fks_phase_state
  public :: initialize_fks_amplitude_state, initialize_fks_config_state
  public :: initialize_fks_pineappl_state, initialize_fks_generated_state
  public :: validate_fks_singular_state

contains

  double precision function evaluate_fks_sij(event_slot, p, ii_fks, &
                                             jj_fks, xi_i_fks, y_ij_fks)
    implicit none
    double precision, intent(in) :: p(0:3, nexternal)
    double precision, intent(in) :: xi_i_fks, y_ij_fks
    integer, intent(in) :: event_slot, ii_fks, jj_fks

    call require_fks_singular_state()
    call initialize_fks_sij_module(nexternal, nincoming, fks_a, fks_b, &
                                   a_h_damp, one_h_damp, useenergy, usebeta)
    call set_fks_sij_partition_state(fks_j_from_i, particle_type, is_aorg, &
                                     i_fks, j_fks, &
                                     ybst_til_tocm(event_slot), &
                                     event_sqrt_shat(event_slot), &
                                     event_shat(event_slot), &
                                     event_fks_momentum(:, &
                                       soft_counterevent:soft_collinear_counterevent), &
                                     external_masses)
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
                                        p_born_norad_in, event_momenta_in, &
                                        event_jacobian_in, idup_in, mothup_in, icolup_in, niprocs_in, &
                                        is_aorg_in, amp2_in, jamp2_in, subproc_pd_in, subproc_iproc_in, flavour_map_in, &
                                        iproc_save_in, eto_in, etoi_in, maxproc_found_in)
    implicit none
    double precision, target, intent(inout) :: p_born_in(0:, 1:)
    double precision, target, intent(inout) :: p_born_coll_in(0:, 1:)
    double precision, target, intent(inout) :: p_born_norad_in(0:, 1:)
    double precision, target, intent(inout) :: event_momenta_in(0:, 1:, 0:)
    double precision, target, intent(inout) :: event_jacobian_in(0:)
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
    if (size(event_momenta_in, 1) /= 4 .or. &
        size(event_momenta_in, 2) /= nexternal .or. &
        size(event_momenta_in, 3) /= real_event + 1 .or. &
        size(event_jacobian_in) /= real_event + 1) then
      call fail_fks_singular_state('invalid event storage shape')
    end if
    p_born => p_born_in
    p_born_coll => p_born_coll_in
    p_born_norad => p_born_norad_in
    stored_event_momenta => event_momenta_in
    stored_event_jacobian => event_jacobian_in
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
    if (.not. associated(p_born) .or. &
        .not. associated(stored_event_momenta) .or. &
        .not. associated(stored_event_jacobian) .or. &
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
    if (event_xi_hat(real_event)*event_xi_max(soft_counterevent) .gt. &
        xiBSVcut_used) return
    call sborn(p_born, wgt_c)
    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      wgt1 = amp_split(iamp)*f_b/g**(qcd_power)
      call add_wgt(soft_counterevent, 2, orders, wgt1, 0d0, 0d0)
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
    if (event_xi_hat(real_event)*event_xi_max(soft_counterevent) .gt. &
        xiBSVcut_used) return
    call bornsoftvirtual(soft_counterevent, &
                         stored_event_momenta(:, :, soft_counterevent), &
                         bsv_wgt, virt_wgt, born_wgt)
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
      call add_wgt(soft_counterevent, 3, orders, wgt1, wgt2, wgt3)
      call add_wgt(soft_counterevent, 15, orders, wgt4, 0d0, 0d0)
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
      call add_wgt(soft_counterevent, 14, orders, wgt1, 0d0, 0d0)
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
    s_ev = evaluate_fks_sij(real_event, p, i_fks, j_fks, &
                            event_xi(real_event), event_y(real_event))
    if (s_ev .le. 0.d0) return
    call sreal(real_event, p, event_xi(real_event), &
               event_y(real_event), fx_ev)
    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      wgt1 = amp_split(iamp)*s_ev*f_r/g**(qcd_power)
      call add_wgt(real_event, 1, orders, wgt1, 0d0, 0d0)
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
    if (event_xi_hat(real_event)*event_xi_max(soft_counterevent) .gt. &
        xiScut_used) return
    s_s = evaluate_fks_sij(soft_counterevent, &
            stored_event_momenta(:, :, soft_counterevent), &
            i_fks, j_fks, zero, event_y(real_event))
    if (s_s .le. 0d0) return
    call sreal(soft_counterevent, &
               stored_event_momenta(:, :, soft_counterevent), &
               0d0, event_y(real_event), fx_s)

    do iamp = 1, amp_split_size
      if (amp_split(iamp) .eq. 0d0) cycle
      call amp_split_pos_to_orders(iamp, orders)
      QCD_power = orders(qcd_pos)
      orders_tag = get_orders_tag(orders)
      amp_pos = iamp
      g22 = g**(QCD_power)
      wgt1 = 0d0
      if (event_xi(real_event) .le. xiScut_used) then
        wgt1 = -amp_split(iamp)*s_s*f_s/g22
      end if
      if (wgt1 .ne. 0d0) &
        call add_wgt(soft_counterevent, 4, orders, wgt1, 0d0, 0d0)
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
    if (event_y(real_event) .le. 1d0 - deltaS .or. &
        pmass(j_fks) .ne. 0.d0) return
    s_c = evaluate_fks_sij(collinear_counterevent, &
            stored_event_momenta(:, :, collinear_counterevent), &
            i_fks, j_fks, event_xi(collinear_counterevent), one)
    if (s_c .le. 0d0) return
! sreal_deg should be called **BEFORE** sreal
! in order not to overwrtie the amp_split array
    call sreal_deg(collinear_counterevent, &
                   stored_event_momenta(:, :, collinear_counterevent), &
                   event_xi(collinear_counterevent), deg_xi_c, deg_lxi_c)
    call sreal(collinear_counterevent, &
               stored_event_momenta(:, :, collinear_counterevent), &
               event_xi(collinear_counterevent), one, fx_c)

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
                      amp_split_wgtpsch_l(iamp))* &
                     log(event_xi(collinear_counterevent)))*f_dc/g22
      wgt3 = amp_split_wgtdegrem_muF(iamp)*f_dc/g22
      if (wgt1 .ne. 0d0 .or. wgt3 .ne. 0d0) &
        call add_wgt(collinear_counterevent, 5, orders, wgt1, 0d0, wgt3)
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
    if (event_xi_hat(real_event)*event_xi_max(collinear_counterevent) &
        .ge. xiScut_used .or. event_y(real_event) .le. 1d0 - deltaS &
        .or. pmass(j_fks) .ne. 0.d0) return
    s_sc = evaluate_fks_sij(soft_collinear_counterevent, &
             stored_event_momenta(:, :, soft_collinear_counterevent), &
             i_fks, j_fks, zero, one)
    if (s_sc .le. 0d0) return
! sreal_deg should be called **BEFORE** sreal
! in order not to overwrtie the amp_split array
    call sreal_deg(soft_collinear_counterevent, &
      stored_event_momenta(:, :, soft_collinear_counterevent), &
      zero, deg_xi_sc, deg_lxi_sc)
    call sreal(soft_collinear_counterevent, &
               stored_event_momenta(:, :, soft_collinear_counterevent), &
               zero, one, fx_sc)

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
      if (event_xi(collinear_counterevent) .lt. xiScut_used) then
        wgt1 = amp_split(iamp)*s_sc*f_sc/g22
        wgt1 = wgt1 + ( &
               -(amp_split_wgtdegrem_xi(iamp) + amp_split_wgtpsch_p(iamp) + &
                 (amp_split_wgtdegrem_lxi(iamp) + amp_split_wgtpsch_l(iamp))* &
                 log(event_xi(collinear_counterevent)))*f_dsc(1) &
               - (amp_split_wgtdegrem_xi(iamp)*f_dsc(2) + &
                  amp_split_wgtdegrem_lxi(iamp)*f_dsc(3)) &
               + amp_split_wgtpsch_d(iamp)*f_pdfsch_d &
               + amp_split_wgtpsch_p(iamp)*f_pdfsch_p &
               + amp_split_wgtpsch_l(iamp)*f_pdfsch_l)/g22
        wgt3 = -amp_split_wgtdegrem_muF(iamp)*f_dsc(4)/g22
      end if
      if (wgt1 .ne. 0d0 .or. wgt3 .ne. 0d0) &
        call add_wgt(soft_collinear_counterevent, 6, orders, wgt1, 0d0, wgt3)
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
        rndec(i) = random_unit_interval(iconfig)
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
! f_* multiplication factors for Born and nbody
    f_b = stored_event_jacobian(soft_counterevent)* &
          event_xi_norm(real_event)/ &
          (min(event_xi_max(real_event), xiBSVcut_used)* &
           event_shat(soft_counterevent)/(16*pi**2))* &
          fkssymmetryfactorBorn*vegas_wgt
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
        call set_alphas(stored_event_momenta(:, :, real_event))
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
    prefact = event_xi_norm(real_event)/event_xi(real_event)/ &
              (1 - event_y(real_event))
    f_r = prefact*jac_ev*fkssymmetryfactor*vegas_wgt
    if (.not. nocntevents) then
      prefact_cnt_ssc = event_xi_norm(real_event)/ &
                        min(event_xi_max(real_event), xiScut_used)* &
                        log(xicut_used/min(event_xi_max(real_event), &
                                          xiScut_used))/ &
                        (1 - event_y(real_event))
      f_s = (prefact + prefact_cnt_ssc)* &
            stored_event_jacobian(soft_counterevent)* &
            fkssymmetryfactor*vegas_wgt
      if (pmass(j_fks) .eq. 0d0) then
! For the soft-collinear, these should be itwo. But they are always
! equal to ione, so no need to define separate factors.
        prefact_c = event_xi_norm(collinear_counterevent)/ &
                    event_xi(collinear_counterevent)/ &
                    (1 - event_y(real_event))
        prefact_coll = event_xi_norm(collinear_counterevent)/ &
                       event_xi(collinear_counterevent)* &
                       log(delta_used/deltaS)/deltaS
        f_c = (prefact_c + prefact_coll)* &
              stored_event_jacobian(collinear_counterevent)* &
              fkssymmetryfactor*vegas_wgt
        prefact_deg = event_xi_norm(collinear_counterevent)/ &
                      event_xi(collinear_counterevent)/deltaS
        prefact_cnt_ssc_c = event_xi_norm(collinear_counterevent)/ &
                            min(event_xi_max(collinear_counterevent), &
                                xiScut_used) &
                            *log(xicut_used/ &
                                 min(event_xi_max(collinear_counterevent), &
                                     xiScut_used)) &
                            /(1 - event_y(real_event))
        prefact_coll_c = event_xi_norm(collinear_counterevent)/ &
                         min(event_xi_max(collinear_counterevent), &
                             xiScut_used) &
                         *log(xicut_used/ &
                              min(event_xi_max(collinear_counterevent), &
                                  xiScut_used)) &
                         *log(delta_used/deltaS)/deltaS
        f_dc = stored_event_jacobian(collinear_counterevent)*prefact_deg/ &
               (event_shat(collinear_counterevent)/(32*pi**2))* &
               fkssymmetryfactorDeg*vegas_wgt
        f_sc = (prefact_c + prefact_coll + prefact_cnt_ssc_c + &
                prefact_coll_c)* &
               stored_event_jacobian(soft_collinear_counterevent)* &
               fkssymmetryfactorDeg*vegas_wgt
        prefact_deg_sxi = event_xi_norm(collinear_counterevent)/ &
                          min(event_xi_max(collinear_counterevent), &
                              xiScut_used)* &
                          log(xicut_used/ &
                              min(event_xi_max(collinear_counterevent), &
                                  xiScut_used))*1/deltaS
        prefact_deg_slxi = event_xi_norm(collinear_counterevent)/ &
                           min(event_xi_max(collinear_counterevent), &
                               xiScut_used) &
                           *(log(xicut_used)**2 &
                             - log(min(event_xi_max(collinear_counterevent), &
                                       xiScut_used))**2) &
                           /(2.d0*deltaS)
        f_dsc(1) = prefact_deg* &
                   stored_event_jacobian(soft_collinear_counterevent)/ &
                   (event_shat(soft_collinear_counterevent)/(32*pi**2))* &
                   fkssymmetryfactorDeg*vegas_wgt
        f_dsc(2) = prefact_deg_sxi* &
                   stored_event_jacobian(soft_collinear_counterevent)/ &
                   (event_shat(soft_collinear_counterevent)/(32*pi**2))* &
                   fkssymmetryfactorDeg*vegas_wgt
        f_dsc(3) = prefact_deg_slxi* &
                   stored_event_jacobian(soft_collinear_counterevent)/ &
                   (event_shat(soft_collinear_counterevent)/(32*pi**2))* &
                   fkssymmetryfactorDeg*vegas_wgt
        f_dsc(4) = (prefact_deg + prefact_deg_sxi)* &
                   stored_event_jacobian(soft_collinear_counterevent)/ &
                   (event_shat(soft_collinear_counterevent)/(32*pi**2))* &
                   fkssymmetryfactorDeg*vegas_wgt
! prefactor for the PDF scheme
        prefact_pdfsch_d = event_xi_norm(collinear_counterevent)/ &
                           xiScut_used/deltaS
        f_pdfsch_d = prefact_pdfsch_d* &
                     stored_event_jacobian(soft_collinear_counterevent)/ &
                     (event_shat(soft_collinear_counterevent)/(32*pi**2))* &
                     fkssymmetryfactorDeg*vegas_wgt
        prefact_pdfsch_p = event_xi_norm(collinear_counterevent)* &
                           dlog(xiScut_used)/xiScut_used/deltaS
        f_pdfsch_p = prefact_pdfsch_p* &
                     stored_event_jacobian(soft_collinear_counterevent)/ &
                     (event_shat(soft_collinear_counterevent)/(32*pi**2))* &
                     fkssymmetryfactorDeg*vegas_wgt
        prefact_pdfsch_l = event_xi_norm(collinear_counterevent)* &
                           dlog(xiScut_used)**2/2d0/xiScut_used/deltaS
        f_pdfsch_l = prefact_pdfsch_l* &
                     stored_event_jacobian(soft_collinear_counterevent)/ &
                     (event_shat(soft_collinear_counterevent)/(32*pi**2))* &
                     fkssymmetryfactorDeg*vegas_wgt
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

  subroutine add_wgt(event_slot, type, orders, wgt1, wgt2, wgt3)
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
    integer event_slot, type, i, j
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

    bjx(1, icontr) = event_bjorken_x(1, event_slot)
    bjx(2, icontr) = event_bjorken_x(2, event_slot)
    scales2(1, icontr) = QES2
    scales2(2, icontr) = scale**2
    scales2(3, icontr) = q2fact(1)
    g_strong(icontr) = g
    nFKS(icontr) = nFKSprocess
    y_bst(icontr) = ybst_til_tolab(event_slot)
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
        if (stored_event_momenta(0, 1, soft_counterevent) .gt. 0d0 &
            .and. type .ne. 5) then
          momenta_m(j, i, 1, icontr) = &
            stored_event_momenta(j, i, soft_counterevent)
        elseif (stored_event_momenta(0, 1, collinear_counterevent) &
                .gt. 0d0) then
          momenta_m(j, i, 1, icontr) = &
            stored_event_momenta(j, i, collinear_counterevent)
        elseif (stored_event_momenta(0, 1, soft_collinear_counterevent) &
                .gt. 0d0) then
          momenta_m(j, i, 1, icontr) = &
            stored_event_momenta(j, i, soft_collinear_counterevent)
        else
          if (i .lt. fks_i_d(nFKSprocess)) then
            momenta_m(j, i, 1, icontr) = p_born(j, i)
          elseif (i .eq. fks_i_d(nFKSprocess)) then
            momenta_m(j, i, 1, icontr) = 0d0
          else
            momenta_m(j, i, 1, icontr) = p_born(j, i - 1)
          end if
        end if
        momenta_m(j, i, 2, icontr) = &
          stored_event_momenta(j, i, real_event)
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
    double precision xlum, mu2_r, mu2_f, mu2_q, wgt_wo_pdf, conv
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
      mu2_q = scales2(1, i)
      mu2_r = scales2(2, i)
      mu2_f = scales2(3, i)
      q2fact(1) = mu2_f
      q2fact(2) = mu2_f
! call the PDFs
      xlum = dlum(bjx(:, i))
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
    double precision xlum(maxscales), pi, mu2_r(maxscales), c_mu2_r, c_mu2_f, mu2_f(maxscales), mu2_q, g(maxscales), conv
    parameter(pi=3.1415926535897932385d0)
    parameter(conv=389379660d0) ! conversion to picobarns
    call cpu_time(tBefore)
    if (icontr .eq. 0) return
! currently we have 'iwgt' weights in the wgts() array.
    iwgt_save = iwgt
! loop over all the contributions in the weight lines module
    do i = 1, icontr
      iwgt = iwgt_save
      nFKSprocess = nFKS(i)
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
          xlum(kf) = dlum(bjx(:, i))
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
    double precision xlum, pi, mu2_r, mu2_f, mu2_q, g, conv
    parameter(pi=3.1415926535897932385d0)
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
          mu2_q = scales2(1, i)
          mu2_r = scales2(2, i)
          mu2_f = scales2(3, i)
          q2fact(1) = mu2_f
          q2fact(2) = mu2_f
! Compute the luminosity
          xlum = dlum(bjx(:, i))
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


  subroutine sreal(event_slot, pp, xi_i_fks, y_ij_fks, wgt)
! Wrapper for the n+1 contribution. Returns the n+1 matrix element
! squared reduced by the FKS damping factor xi**2*(1-y).
! Close to the soft or collinear limits it calls the corresponding
! Born and multiplies with the AP splitting function or eikonal factors.
    implicit none
    integer, intent(in) :: event_slot
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

! Check that the requested event slot matches the supplied momenta.
    if (nincoming .eq. 2) then
      shattmp = 2d0*dot(pp(0, 1), pp(0, 2))
    else
      shattmp = pp(0, 1)**2
    end if
    if (abs(shattmp/event_shat(event_slot) - 1.d0) .gt. 1.d-5) then
      write (*, *) 'Error in sreal: inconsistent shat'
      write (*, *) shattmp, event_shat(event_slot)
      stop
    end if

    if (1d0 - y_ij_fks .lt. tiny) then
      if (pmass(j_fks) .eq. zero .and. j_fks .le. nincoming) then
        call sborncol_isr(pp, xi_i_fks, y_ij_fks, &
                          event_shat(event_slot), wgt)
      elseif (pmass(j_fks) .eq. zero .and. j_fks .ge. nincoming + 1) then
        call sborncol_fsr(pp, y_ij_fks, event_shat(event_slot), wgt)
      else
        wgt = 0d0
        amp_split(1:amp_split_size) = 0d0
      end if
    elseif (xi_i_fks .lt. tiny) then
      if (need_color_links) then
! has soft singularities
        call sbornsoft(pp, xi_i_fks, y_ij_fks, &
                       event_sqrt_shat(event_slot), wgt)
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

  subroutine sborncol_fsr(p, y_ij_fks, partonic_shat, wgt)
    implicit none
    double precision p(0:3, nexternal), wgt
    double precision y_ij_fks, partonic_shat
!
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
    t = z*partonic_shat/4d0
    call sborn(p_born, wgt_born)
    if (iextra_cnt .gt. 0) call extra_cnt(p_born, iextra_cnt, ans_extra_cnt)
    call AP_reduced(j_type, i_type, t, z, g, ap)
    call Qterms_reduced_timelike(j_type, i_type, t, z, g, Q)
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
          pi(i) = event_fks_momentum(i, real_event)
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

  subroutine sborncol_isr(p, xi_i_fks, y_ij_fks, partonic_shat, wgt)
    implicit none
    double precision p(0:3, nexternal), wgt
    double precision xi_i_fks, y_ij_fks, partonic_shat
!

    double precision p_born_used(0:3, nexternal - 1)
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
    t = z*partonic_shat/4d0
    call AP_reduced(m_type, i_type, t, z, g, ap)
    call Qterms_reduced_spacelike(m_type, i_type, t, z, g, Q)
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
          pi(i) = event_fks_momentum(i, real_event)
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


  subroutine sbornsoft(pp, xi_i_fks, y_ij_fks, partonic_sqrt_shat, wgt)
    implicit none
!      include "fks.inc"
    integer m, n

    double precision softcontr, pp(0:3, nexternal), wgt, eik
    double precision xi_i_fks, y_ij_fks, partonic_sqrt_shat
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
            call eikonal_reduced( &
              pp, m, n, i_fks, j_fks, xi_i_fks, y_ij_fks, &
              event_fks_momentum(:, &
                soft_counterevent:soft_collinear_counterevent), &
              external_masses, partonic_sqrt_shat, eik)
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


  subroutine sreal_deg(event_slot, p, xi_i_fks, collrem_xi, collrem_lxi)
    use extra_weights
    implicit none
    integer event_slot, iord, iap
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

! Check that the requested event slot matches the supplied momenta.
    if (nincoming .eq. 2) then
      shattmp = 2d0*dot(p(0, 1), p(0, 2))
    else
      shattmp = p(0, 1)**2
    end if
    if (abs(shattmp/event_shat(event_slot) - 1.d0) .gt. 1.d-5) then
      write (*, *) 'Error in sreal: inconsistent shat'
      write (*, *) shattmp, event_shat(event_slot)
      stop
    end if

! A factor gS^2 is included in the Altarelli-Parisi kernels
    oo2pi = one/(8d0*PI**2)

    z = 1d0 - xi_i_fks
    t = one
    call AP_reduced(m_type, i_type, t, z, g, ap)
    call AP_reduced_prime(m_type, i_type, t, z, g, apprime)

! call the PDF-scheme kernels here
!   p-> (/1/(1-z)/)_+
!   l-> (/log(1-z)/(1-z)/)_+
!   d-> delta(1-z)
    call xkplus(PDFscheme, m_type, i_type, z, g, nf, xkkernp)
    call xkdelta(PDFscheme, m_type, i_type, g, xkkernd)
    call xklog(PDFscheme, m_type, i_type, z, g, nf, xkkernl)

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

    collrem_xi_tmp = ap(iap)*log( &
                       event_shat(event_slot)*delta_used/ &
                       (2*q2fact(j_fks))) - apprime(iap)
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

    prefact_xi = ap(iap)*log( &
                   event_shat(event_slot)*delta_used/(2*QES2)) - &
                 apprime(iap)
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



  subroutine bornsoftvirtual(event_slot, p, bsv_wgt, virt_wgt, born_wgt)
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
    integer event_slot, i, j, aj, m, n, k



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

! Check that the requested event slot matches the supplied momenta.
    if (nincoming .eq. 2) then
      shattmp = 2d0*dot(p(0, 1), p(0, 2))
    else
      shattmp = p(0, 1)**2
    end if
    if (abs(shattmp/event_shat(event_slot) - 1.d0) .gt. 1.d-5) then
      write (*, *) 'Error in bornsoftvirtual: inconsistent shat'
      write (*, *) shattmp, event_shat(event_slot)
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
                  - dlog(event_shat(event_slot)*deltaO/2d0/QES2) &
                  *(gamma_used &
                    - 2d0*c_used*dlog( &
                        2d0*Ej/xicut_used/event_sqrt_shat(event_slot))) &
                  + 2d0*c_used*(dlog( &
                        2d0*Ej/event_sqrt_shat(event_slot))**2 &
                                - dlog(xicut_used)**2) &
                  - 2d0*gamma_used*dlog( &
                      2d0*Ej/event_sqrt_shat(event_slot))
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
              call eikonal_Ireg(p, m, n, xicut_used, external_masses, &
                                 event_shat(event_slot), qes2, abrv, &
                                 eikIreg)
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
      if ((random_unit_interval(iconfig) .le. virtual_fraction(ichan) &
           .and. abrv(1:3) .ne. 'nov') .or. &
          abrv(1:4) .eq. 'virt') then
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
