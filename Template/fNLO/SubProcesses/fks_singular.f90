module fks_singular_module
  use process_dimensions, only: nexternal, nincoming, max_branch, &
                                lmaxconfigs, fks_configs, nsplitorders, &
                                qcd_pos, amp_split_size, &
                                validate_process_dimensions
  use run_state, only: q2fact, scale
  use timing_state, only: tOLP
  use mint_module, only: maxchannels
  use fks_metadata, only: validate_fks_metadata
  use split_orders, only: amp_split_pos_to_orders
  use kin_functions_module, only: dot => dot_impl
  use phase_space_kinematics, only: getaziangles
  use factorized_phase_space, only: factorized_radiation_state, &
       fetch_factorized_kernel_momenta, &
       fetch_factorized_radiation_state
  use decay_chain_metadata, only: has_decay_chains, context_for_fks, &
       born_context, context_core_count, core_leg_pdg
  use decay_chain_kinematics, only: get_core_born_momenta, &
       get_core_mass_buffer, active_core_count, &
       map_core_color_pair
  use decay_chain_scales, only: corrected_born_qcd_squared_order
  use nlo_decay_metadata, only: has_nlo_decay, &
       nlo_decay_born_context, nlo_decay_context_for_fks, &
       nlo_decay_local_count, nlo_decay_local_pdg, &
       nlo_decay_local_is_final, nlo_decay_fks_i, nlo_decay_fks_j, &
       nlo_decay_partner_count, nlo_decay_partner_local, &
       nlo_decay_map_color_link, nlo_decay_corrected_node
  use nlo_decay_kinematics, only: get_nlo_decay_born_kernel, &
       get_nlo_decay_mass_buffer, nlo_decay_parent_mass
  use nlo_contribution_bundle, only: has_nlo_contribution_bundle, &
       active_nlo_contribution, &
       active_contribution_fks_first, active_contribution_fks_last, &
       active_contribution_has_virtual, active_virtual_grid_index
  use fks_qcd_splitting, only: AP_reduced, AP_reduced_prime, &
                                Qterms_reduced_timelike, &
                                Qterms_reduced_spacelike
  use fks_random_module, only: random_unit_interval
  use fks_soft_kernels, only: eikonal_reduced, eikonal_Ireg
  use fks_sij_module, only: initialize_fks_sij_module, &
                            set_fks_sij_partition_state, fks_sij_impl
  use FKSParams, only: use_poly_virtual
  use chooser_functions_module, only: get_mother_colour_impl, set_pdg_impl
  use fks_model_state_module, only: g => strong_coupling, active_flavours, &
                                    external_masses, validate_fks_model_state
  use spin_density_fks_matrices, only: &
       spin_density_fks_collection_enabled, &
       reset_spin_density_born_matrix, &
       reset_spin_density_real_matrix, &
       reset_spin_density_color_matrix, &
       spin_density_born_matrix_available, &
       spin_density_real_matrix_available, &
       spin_density_color_matrix_available, &
       reset_spin_density_reduced_matrix, &
       set_spin_density_reduced_from_real, &
       set_spin_density_reduced_from_born, &
       add_spin_density_reduced_color, &
       reset_spin_density_degenerate_matrix, &
       set_spin_density_degenerate_from_born, &
       reset_spin_density_integrated_matrix, &
       add_spin_density_integrated_born, &
       add_spin_density_integrated_color, &
       reset_spin_density_virtual_matrix, &
       reduce_spin_density_virtual_matrix, &
       spin_density_virtual_matrix_available, &
       spin_density_nlo_amp_position
  use fnlo_process_common, only: nfksprocess, i_fks, j_fks, &
                                 soft_counterevent, &
                                 soft_collinear_counterevent, real_event, &
                                 event_fks_momentum, &
                                 event_sqrt_shat, event_shat, &
                                 ybst_til_tocm, xiscut_used, &
                                 xibsvcut_used, delta_used, xicut_used, &
                                 fkssymmetryfactor, &
                                 fkssymmetryfactorborn, &
                                 fkssymmetryfactordeg, ngluons, &
                                 diagramsymmetryfactor, &
                                 calculatedborn => calculated_born, &
                                 virtual_over_born, softtest, colltest, &
                                 need_color_links, xij_aor, &
                                 i_type, j_type, m_type, &
                                 iden_comp, &
                                 c, gamma, gammap, beta0, abrv, &
                                 multi_channel, nbody, qes2, amp_split, &
                                 amp_split_cnt, p_born, &
                                 stored_event_momenta => event_momenta, &
                                 idup, &
                                 fks_j_from_i, particle_type, pdg_type, &
                                 is_aorg, ans_cnt, &
                                 amp_split_virt, &
                                 amp_split_born_for_virt, amp_split_avv, &
                                 amp_split_wgtnstmp, &
                                 amp_split_wgtwnstmpmuf, &
                                 amp_split_wgtwnstmpmur, &
                                 amp_split_wgtdegrem_xi, &
                                 amp_split_wgtdegrem_lxi, &
                                 amp_split_wgtdegrem_muf, &
                                 amp_split_soft, &
                                 amp_split_finite_ml, amp_split_poles_fks, &
                                 config_mass, &
                                 config_width, config_forest, config_sprop, &
                                 config_tprid, config_map, real_forest, &
                                 real_sprop, real_tprid, real_map, &
                                 real_mass, real_width
  implicit none
  private

  interface
    integer function get_color(ipdg)
      integer, intent(in) :: ipdg
    end function get_color
  end interface

  double precision, parameter :: fks_a = 1.5d0, fks_b = 1.5d0
  double precision, parameter :: a_h_damp = 1d0, one_h_damp = 1d-2
  double precision, parameter :: deltao = 1d0, deltai = 1d0, xicut = 0.5d0
  double precision, parameter :: deltas = 1d0, xiscut = 0.5d0, xibsvcut = 1d0
  double precision, parameter :: deltaminy = 0.95d0, skewy = 10d0, alphay = 2d0

  integer, allocatable, save :: born_forest(:, :, :), born_sprop(:, :)
  integer, allocatable, save :: born_tprid(:, :), born_map(:)
  double precision, allocatable, save :: born_mass(:, :), born_width(:, :)
  integer, save :: born_max_branch_used = 0, born_lmaxconfigs_used = 0

  logical, allocatable, save :: firsttime_nfksprocess(:)
  double precision, allocatable, save :: diagramsymmetryfactor_save(:)
  integer, allocatable, save :: fac_i_fks(:), fac_j_fks(:)
  integer, allocatable, save :: i_type_fks(:), j_type_fks(:), m_type_fks(:)
  integer, allocatable, save :: ngluons_fks(:)
  integer, allocatable, save :: iden_real_fks(:), iden_born_fks(:)
  logical, save :: setfks_firsttime = .true.
  logical, save :: fks_singular_state_initialized = .false.

  public :: evaluate_fks_sij, sreal, sreal_deg, bornsoftvirtual
  public :: evaluate_born_matrix, evaluate_virtual_matrix
  public :: getpoles, setfksfactor, fill_configurations_common
  public :: fks_subtraction_shat
  public :: initialize_fks_generated_state
  public :: validate_fks_singular_state

contains

  double precision function evaluate_fks_sij(event_slot, ii_fks, jj_fks, &
                                             xi_i_fks, y_ij_fks)
    implicit none
    double precision, intent(in) :: xi_i_fks, y_ij_fks
    integer, intent(in) :: event_slot, ii_fks, jj_fks
    double precision :: kernel_momenta(0:3, nexternal)
    double precision :: kernel_masses(nexternal)
    double precision :: kernel_fks_momenta(0:3, &
         soft_counterevent:soft_collinear_counterevent)
    type(factorized_radiation_state) :: radiation

    call require_fks_singular_state()
    if (has_nlo_decay()) then
      evaluate_fks_sij = evaluate_nlo_decay_fks_sij( &
           event_slot, ii_fks, jj_fks, xi_i_fks, y_ij_fks)
      return
    end if
    call select_kernel_event(event_slot, kernel_momenta, kernel_masses)
    call load_kernel_radiation(event_slot, radiation)
    call load_counterevent_fks_momenta(kernel_fks_momenta)
    call initialize_fks_sij_module(nexternal, nincoming, fks_a, fks_b, &
                                   a_h_damp, one_h_damp)
    call set_fks_sij_partition_state(fks_j_from_i, particle_type, is_aorg, &
                                     i_fks, j_fks, &
                                     radiation%y_to_cm, &
                                     radiation%sqrt_shat, &
                                     radiation%shat, &
                                     kernel_fks_momenta, &
                                     kernel_masses)
    evaluate_fks_sij = fks_sij_impl(kernel_momenta, ii_fks, jj_fks, &
                                    xi_i_fks, y_ij_fks)
  end function evaluate_fks_sij


  double precision function evaluate_nlo_decay_fks_sij(event_slot, &
       ii_fks, jj_fks, xi_i_fks, y_ij_fks)
    integer, intent(in) :: event_slot, ii_fks, jj_fks
    double precision, intent(in) :: xi_i_fks, y_ij_fks
    integer :: context, local_count, local_i, local_j
    integer :: configuration, emitted, partner_count, position, leg
    integer :: local_partner_map(nexternal, 0:nexternal)
    integer :: local_particle_type(nexternal)
    logical :: local_is_aorg(nexternal)
    double precision :: local_event(0:3, nexternal)
    double precision :: local_masses(nexternal)
    double precision :: local_fks_momenta(0:3, 0:2)
    double precision :: decay_mass
    type(factorized_radiation_state) :: radiation
    logical :: available

    if (ii_fks /= i_fks .or. jj_fks /= j_fks) then
      call fail_fks_singular_state( &
           'NLO-decay S function requested for inactive visible indices')
    end if
    context = nlo_decay_context_for_fks(nfksprocess)
    local_count = nlo_decay_local_count(context)
    local_i = nlo_decay_fks_i(nfksprocess)
    local_j = nlo_decay_fks_j(nfksprocess)
    local_partner_map = 0
    local_particle_type = 1
    local_is_aorg = .false.

    ! Rebuild the same decay-local partner table written for the standalone
    ! 1 -> n process.  Configurations sharing a real matrix element share a
    ! row, including any massive incoming decay parent.
    do configuration = active_contribution_fks_first(), &
                       active_contribution_fks_last()
      if (nlo_decay_context_for_fks(configuration) /= context) cycle
      emitted = nlo_decay_fks_i(configuration)
      partner_count = nlo_decay_partner_count(configuration)
      if (local_partner_map(emitted, 0) == 0) then
        local_partner_map(emitted, 0) = partner_count
        do position = 1, partner_count
          local_partner_map(emitted, position) = &
               nlo_decay_partner_local(configuration, position)
        end do
      else
        if (local_partner_map(emitted, 0) /= partner_count) then
          call fail_fks_singular_state( &
               'inconsistent NLO-decay partner rows')
        end if
        do position = 1, partner_count
          if (local_partner_map(emitted, position) /= &
              nlo_decay_partner_local(configuration, position)) then
            call fail_fks_singular_state( &
                 'inconsistent NLO-decay partner ordering')
          end if
        end do
      end if
    end do

    do leg = 1, local_count
      local_particle_type(leg) = get_color(nlo_decay_local_pdg(context, leg))
      local_is_aorg(leg) = abs(nlo_decay_local_pdg(context, leg)) == 21
    end do
    local_event = 0d0
    call fetch_factorized_kernel_momenta( &
         event_slot, active_factorized_block(), local_count, &
         local_event(:, 1:local_count), available)
    if (.not. available) then
      call fail_fks_singular_state( &
           'decay-local S-function momenta are unavailable')
    end if
    call get_nlo_decay_mass_buffer(nfksprocess, local_masses)
    call load_counterevent_fks_momenta(local_fks_momenta)
    call load_kernel_radiation(event_slot, radiation)
    decay_mass = radiation%sqrt_shat
    call initialize_fks_sij_module(local_count, 1, fks_a, fks_b, &
                                   a_h_damp, one_h_damp, .true.)
    call set_fks_sij_partition_state( &
         local_partner_map(1:local_count, 0:local_count), &
         local_particle_type(1:local_count), &
         local_is_aorg(1:local_count), &
         local_i, local_j, 0d0, decay_mass, decay_mass**2, &
         local_fks_momenta, local_masses(1:local_count))
    evaluate_nlo_decay_fks_sij = fks_sij_impl( &
         local_event(:, 1:local_count), local_i, local_j, &
         xi_i_fks, y_ij_fks)
  end function evaluate_nlo_decay_fks_sij

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
    call allocate_fks_singular_caches()
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
    if (.not. allocated(firsttime_nfksprocess)) allocate (firsttime_nfksprocess(fks_configs))
    if (.not. allocated(diagramsymmetryfactor_save)) allocate (diagramsymmetryfactor_save(maxchannels))
    if (.not. allocated(fac_i_fks)) then
      allocate (fac_i_fks(fks_configs), fac_j_fks(fks_configs))
      allocate (i_type_fks(fks_configs), j_type_fks(fks_configs), m_type_fks(fks_configs))
      allocate (ngluons_fks(fks_configs))
      allocate (iden_real_fks(fks_configs), iden_born_fks(fks_configs))
    end if
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
    call validate_fks_model_state()
    if (.not. allocated(born_forest)) then
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


  logical function uses_factorized_kernel_state()
    if (has_nlo_decay()) then
      uses_factorized_kernel_state = .true.
    else
      uses_factorized_kernel_state = has_decay_chains()
    end if
  end function uses_factorized_kernel_state


  integer function active_factorized_block()
    if (has_nlo_decay()) then
      active_factorized_block = nlo_decay_corrected_node()
    else
      active_factorized_block = 0
    end if
  end function active_factorized_block


  subroutine evaluate_born_matrix(event_slot, weight)
    integer, intent(in) :: event_slot
    double precision, intent(out) :: weight
    double precision :: legacy_momenta(0:3, nexternal - 1)

    legacy_momenta = 0d0
    if (uses_factorized_kernel_state()) then
      if (event_slot /= soft_counterevent) then
        call fail_fks_singular_state( &
             'a factorized Born matrix requested a non-Born event slot')
      end if
    else
      legacy_momenta = p_born
    end if
    if (uses_factorized_kernel_state()) then
      if (has_nlo_contribution_bundle() .and. &
          spin_density_fks_collection_enabled()) &
        call reset_spin_density_born_matrix()
      call sborn_factorized( &
           active_nlo_contribution(), event_slot, weight)
      if (multi_channel) then
        call sborn_factorized_channel_weights( &
             p_born)
      end if
      if (has_nlo_contribution_bundle() .and. &
          spin_density_fks_collection_enabled()) then
        if (.not. spin_density_born_matrix_available()) then
          call fail_fks_singular_state( &
               'the generated Born did not publish its block density')
        end if
      end if
    else
      call sborn(legacy_momenta, weight)
    end if
  end subroutine evaluate_born_matrix


  subroutine evaluate_born_color_matrix(event_slot, first, second, weight)
    integer, intent(in) :: event_slot, first, second
    double precision, intent(out) :: weight
    double precision :: legacy_momenta(0:3, nexternal - 1)

    legacy_momenta = 0d0
    if (uses_factorized_kernel_state()) then
      if (event_slot /= soft_counterevent) then
        call fail_fks_singular_state( &
             'a factorized color matrix requested a non-Born event slot')
      end if
    else
      legacy_momenta = p_born
    end if
    if (uses_factorized_kernel_state()) then
      if (has_nlo_contribution_bundle() .and. &
          spin_density_fks_collection_enabled()) &
        call reset_spin_density_color_matrix()
      call sborn_sf_factorized(active_nlo_contribution(), event_slot, &
                               first, second, weight)
      if (has_nlo_contribution_bundle() .and. &
          spin_density_fks_collection_enabled()) then
        if (.not. spin_density_color_matrix_available()) then
          call fail_fks_singular_state( &
               'the generated color Born did not publish its block density')
        end if
      end if
    else
      call sborn_sf(legacy_momenta, first, second, weight)
    end if
  end subroutine evaluate_born_color_matrix


  subroutine evaluate_real_matrix(event_slot, weight)
    integer, intent(in) :: event_slot
    double precision, intent(out) :: weight
    double precision :: legacy_momenta(0:3, nexternal)

    legacy_momenta = 0d0
    if (uses_factorized_kernel_state()) then
      if (event_slot /= real_event) then
        call fail_fks_singular_state( &
             'a factorized real matrix requested a counterevent slot')
      end if
    else
      legacy_momenta = stored_event_momenta(:, :, event_slot)
    end if
    if (uses_factorized_kernel_state()) then
      if (has_nlo_contribution_bundle() .and. &
          spin_density_fks_collection_enabled()) &
        call reset_spin_density_real_matrix()
      call smatrix_real_factorized(nfksprocess, event_slot, weight)
      if (has_nlo_contribution_bundle() .and. &
          spin_density_fks_collection_enabled()) then
        if (.not. spin_density_real_matrix_available()) then
          call fail_fks_singular_state( &
               'the generated real did not publish its block density')
        end if
      end if
    else
      call smatrix_real(legacy_momenta, weight)
    end if
  end subroutine evaluate_real_matrix


  subroutine evaluate_virtual_matrix(event_slot, born_weight, virtual_weight)
    integer, intent(in) :: event_slot
    double precision, intent(inout) :: born_weight
    double precision, intent(out) :: virtual_weight
    double precision :: legacy_momenta(0:3, nexternal - 1)

    legacy_momenta = 0d0
    if (uses_factorized_kernel_state()) then
      if (event_slot /= soft_counterevent) then
        call fail_fks_singular_state( &
             'a factorized virtual matrix requested a non-Born event slot')
      end if
    else
      legacy_momenta = p_born
    end if
    if (uses_factorized_kernel_state()) then
      call BinothLHA_factorized( &
           active_nlo_contribution(), event_slot, born_weight, &
           virtual_weight)
    else
      call BinothLHA(legacy_momenta, born_weight, virtual_weight)
    end if
  end subroutine evaluate_virtual_matrix


  subroutine load_kernel_radiation(event_slot, radiation)
    integer, intent(in) :: event_slot
    type(factorized_radiation_state), intent(out) :: radiation
    logical :: available

    if (uses_factorized_kernel_state()) then
      call fetch_factorized_radiation_state( &
           event_slot, active_factorized_block(), radiation, available)
      if (.not. available) then
        call fail_fks_singular_state( &
             'block-local radiation state is unavailable')
      end if
    else
      radiation = factorized_radiation_state()
      radiation%fks_momentum = event_fks_momentum(:, event_slot)
      radiation%sqrt_shat = event_sqrt_shat(event_slot)
      radiation%shat = event_shat(event_slot)
      radiation%y_to_cm = ybst_til_tocm(event_slot)
    end if
  end subroutine load_kernel_radiation


  subroutine load_counterevent_fks_momenta(momenta)
    double precision, intent(out) :: momenta(0:3, &
         soft_counterevent:soft_collinear_counterevent)
    type(factorized_radiation_state) :: radiation
    integer :: event_slot
    logical :: available

    if (.not. uses_factorized_kernel_state()) then
      momenta = event_fks_momentum(:, &
           soft_counterevent:soft_collinear_counterevent)
      return
    end if
    momenta = -1d0
    do event_slot = soft_counterevent, soft_collinear_counterevent
      call fetch_factorized_radiation_state( &
           event_slot, active_factorized_block(), radiation, available)
      if (available) momenta(:, event_slot) = radiation%fks_momentum
    end do
  end subroutine load_counterevent_fks_momenta


  subroutine select_kernel_event(event_slot, kernel_momenta, kernel_masses, &
                                 kernel_count)
    integer, intent(in) :: event_slot
    double precision, intent(out) :: kernel_momenta(0:3, nexternal)
    double precision, intent(out) :: kernel_masses(nexternal)
    integer, intent(out), optional :: kernel_count

    logical :: available
    integer :: selected_count

    if (has_nlo_decay()) then
      selected_count = nlo_decay_local_count( &
           nlo_decay_context_for_fks(nfksprocess))
      if (present(kernel_count)) then
        kernel_count = selected_count
      end if
      kernel_momenta = 0d0
      call fetch_factorized_kernel_momenta( &
           event_slot, active_factorized_block(), selected_count, &
           kernel_momenta(:, 1:selected_count), available)
      if (.not. available) then
        call fail_fks_singular_state( &
             'decay-local kernel momenta are unavailable')
      end if
      call get_nlo_decay_mass_buffer(nfksprocess, kernel_masses)
    else if (has_decay_chains()) then
      selected_count = active_core_count(nfksprocess)
      if (present(kernel_count)) then
        kernel_count = selected_count
      end if
      kernel_momenta = 0d0
      call fetch_factorized_kernel_momenta( &
           event_slot, 0, selected_count, &
           kernel_momenta(:, 1:selected_count), available)
      if (.not. available) then
        call fail_fks_singular_state( &
             'production-kernel momenta are unavailable')
      end if
      call get_core_mass_buffer(context_for_fks(nfksprocess), &
                                kernel_masses)
    else
      if (present(kernel_count)) kernel_count = nexternal
      kernel_momenta = stored_event_momenta(:, :, event_slot)
      kernel_masses = external_masses
    end if
  end subroutine select_kernel_event


  double precision function fks_subtraction_shat(event_slot)
    integer, intent(in) :: event_slot
    type(factorized_radiation_state) :: radiation
    logical :: available
    if (uses_factorized_kernel_state()) then
      call fetch_factorized_radiation_state( &
           event_slot, active_factorized_block(), radiation, available)
      if (available) then
        fks_subtraction_shat = radiation%shat
      else if (has_nlo_decay()) then
        ! Massive decay emitters do not have collinear counterevents.
        fks_subtraction_shat = nlo_decay_parent_mass()**2
      else
        ! A massive production mapping can omit counterevents as well.  Its
        ! physical real event is still owned by the production block.
        call fetch_factorized_radiation_state( &
             real_event, active_factorized_block(), radiation, available)
        if (.not. available) then
          call fail_fks_singular_state( &
               'production radiation scale is unavailable')
        end if
        fks_subtraction_shat = radiation%shat
      end if
    else
      fks_subtraction_shat = event_shat(event_slot)
    end if
  end function fks_subtraction_shat


  subroutine select_kernel_born(kernel_momenta)
    double precision, intent(out) :: kernel_momenta(0:3, nexternal - 1)

    if (has_decay_chains()) then
      call get_core_born_momenta(kernel_momenta)
    else
      kernel_momenta = p_born
    end if
  end subroutine select_kernel_born


  subroutine select_kernel_masses(kernel_masses, kernel_count)
    double precision, intent(out) :: kernel_masses(nexternal)
    integer, intent(out), optional :: kernel_count

    if (has_decay_chains()) then
      if (present(kernel_count)) then
        kernel_count = active_core_count(nfksprocess)
      end if
      call get_core_mass_buffer(context_for_fks(nfksprocess), &
                                kernel_masses)
    else
      if (present(kernel_count)) kernel_count = nexternal
      kernel_masses = external_masses
    end if
  end subroutine select_kernel_masses


  subroutine select_visible_color_pair(core_first, core_second, &
                                       visible_first, visible_second, &
                                       multiplier)
    integer, intent(in) :: core_first, core_second
    integer, intent(out) :: visible_first, visible_second
    double precision, intent(out), optional :: multiplier
    double precision :: local_multiplier

    local_multiplier = 1d0
    if (has_nlo_decay()) then
      call nlo_decay_map_color_link(nfksprocess, core_first, core_second, &
                                    visible_first, visible_second, &
                                    local_multiplier)
    else if (has_decay_chains()) then
      call map_core_color_pair(core_first, core_second, visible_first, &
                               visible_second)
    else
      visible_first = core_first
      visible_second = core_second
    end if
    if (present(multiplier)) multiplier = local_multiplier
  end subroutine select_visible_color_pair


  subroutine select_kernel_properties(kernel_particle_type, kernel_i, &
                                      kernel_initial_count)
    integer, intent(out) :: kernel_particle_type(nexternal)
    integer, intent(out) :: kernel_i, kernel_initial_count
    integer :: context, leg

    if (has_nlo_decay()) then
      context = nlo_decay_context_for_fks(nfksprocess)
      kernel_particle_type = 1
      do leg = 1, nlo_decay_local_count(context)
        kernel_particle_type(leg) = &
             get_color(nlo_decay_local_pdg(context, leg))
      end do
      kernel_i = nlo_decay_fks_i(nfksprocess)
      kernel_initial_count = 1
    else
      kernel_particle_type = particle_type
      kernel_i = i_fks
      kernel_initial_count = nincoming
    end if
  end subroutine select_kernel_properties


  integer function selected_partner_count(kernel_i)
    integer, intent(in) :: kernel_i
    if (has_nlo_decay()) then
      selected_partner_count = nlo_decay_partner_count(nfksprocess)
    else
      selected_partner_count = fks_j_from_i(kernel_i, 0)
    end if
  end function selected_partner_count


  integer function selected_partner(kernel_i, position)
    integer, intent(in) :: kernel_i, position
    if (has_nlo_decay()) then
      selected_partner = &
           nlo_decay_partner_local(nfksprocess, position)
    else
      selected_partner = fks_j_from_i(kernel_i, position)
    end if
  end function selected_partner

  subroutine sreal(event_slot, xi_i_fks, y_ij_fks, wgt)
! Wrapper for the n+1 contribution. Returns the n+1 matrix element
! squared reduced by the FKS damping factor xi**2*(1-y).
! Close to the soft or collinear limits it calls the corresponding
! Born and multiplies with the AP splitting function or eikonal factors.
    implicit none
    integer, intent(in) :: event_slot
    double precision, intent(in) :: xi_i_fks, y_ij_fks
    double precision, intent(out) :: wgt

    double precision shattmp, partonic_shat, partonic_sqrt_shat
    integer :: kernel_i, kernel_j, kernel_initial_count




    double precision zero, tiny
    parameter(zero=0d0)

    double precision pmass(nexternal)
    double precision kernel_momenta(0:3, nexternal)
    if (spin_density_fks_collection_enabled()) then
      call reset_spin_density_reduced_matrix()
    end if
    call select_kernel_event(event_slot, kernel_momenta, pmass)
    partonic_shat = fks_subtraction_shat(event_slot)
    partonic_sqrt_shat = sqrt(partonic_shat)
    if (has_nlo_decay()) then
      kernel_i = nlo_decay_fks_i(nfksprocess)
      kernel_j = nlo_decay_fks_j(nfksprocess)
      kernel_initial_count = 1
    else
      kernel_i = i_fks
      kernel_j = j_fks
      kernel_initial_count = nincoming
    end if
    if (softtest .or. colltest) then
      tiny = 1d-12
    else
      tiny = 1d-6
    end if

    if (kernel_momenta(0, 1) .le. 0.d0) then
! Unphysical kinematics: set matrix elements equal to zero
      wgt = 0.d0
      return
    end if

! Check that the requested event slot matches the supplied momenta.
    if (kernel_initial_count .eq. 2) then
      shattmp = 2d0*dot(kernel_momenta(0, 1), &
                        kernel_momenta(0, 2))
    else
      shattmp = kernel_momenta(0, 1)**2
    end if
    if (abs(shattmp/partonic_shat - 1.d0) .gt. 1.d-5) then
      write (*, *) 'Error in sreal: inconsistent shat'
      write (*, *) shattmp, partonic_shat
      stop
    end if

    if (1d0 - y_ij_fks .lt. tiny) then
      if (pmass(kernel_j) .eq. zero .and. &
          kernel_j .le. kernel_initial_count) then
        call sborncol_isr(kernel_momenta, xi_i_fks, y_ij_fks, &
                          partonic_shat, wgt)
      elseif (pmass(kernel_j) .eq. zero .and. &
              kernel_j .ge. kernel_initial_count + 1) then
        call sborncol_fsr(kernel_momenta, y_ij_fks, &
                          partonic_shat, wgt)
      else
        wgt = 0d0
        amp_split(1:amp_split_size) = 0d0
      end if
    elseif (xi_i_fks .lt. tiny) then
      if (need_color_links) then
! has soft singularities
        call sbornsoft(event_slot, kernel_momenta, xi_i_fks, y_ij_fks, &
                       partonic_sqrt_shat, wgt)
      else
        wgt = 0d0
        amp_split(1:amp_split_size) = 0d0
      end if
    else
      call evaluate_real_matrix(event_slot, wgt)
      if (spin_density_fks_collection_enabled()) then
        call set_spin_density_reduced_from_real( &
             xi_i_fks**2*(1d0 - y_ij_fks))
      end if
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
    integer i, imother_fks, iord, selected_i_fks, selected_j_fks
    integer kernel_count
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
    double precision kernel_born(0:3, nexternal)
    double precision kernel_masses(nexternal)
    type(factorized_radiation_state) :: radiation
    complex(kind=kind(0d0)) wgt1(2)
    complex(kind=8) :: density_ordinary, density_correlated
!
    amp_split_local(1:amp_split_size) = 0d0
    call load_kernel_radiation(real_event, radiation)
    kernel_born = 0d0
    if (has_nlo_decay()) then
      selected_i_fks = nlo_decay_fks_i(nfksprocess)
      selected_j_fks = nlo_decay_fks_j(nfksprocess)
      call get_nlo_decay_born_kernel( &
           nfksprocess, kernel_born, kernel_masses, kernel_count)
    else
      selected_i_fks = i_fks
      selected_j_fks = j_fks
      call select_kernel_born(kernel_born(:, 1:nexternal - 1))
    end if

!
    if (p_born(0, 1) .le. 0.d0) then
! Unphysical kinematics: set matrix elements equal to zero
      write (*, *) "No born momenta in sborncol_fsr"
      wgt = 0.d0
      return
    end if

    E_j_fks = p(0, selected_j_fks)
    E_i_fks = p(0, selected_i_fks)
    z = 1d0 - E_i_fks/(E_i_fks + E_j_fks)
    t = z*partonic_shat/4d0
    call evaluate_born_matrix(soft_counterevent, wgt_born)
    call AP_reduced(j_type, i_type, t, z, g, ap)
    call Qterms_reduced_timelike(j_type, i_type, t, z, g, Q)
    density_ordinary = cmplx(ap(1)*iden_comp, 0d0, kind=8)
    density_correlated = cmplx(Q(1)*iden_comp, 0d0, kind=8)
    wgt = 0d0
    iord = qcd_pos
    wgt1(1) = ans_cnt(1, iord)
    wgt1(2) = ans_cnt(2, iord)
    if (abs(j_type) .eq. 3 .and. i_type .eq. 8) then
      Q(1) = 0d0
      wgt1(2) = 0d0
      density_correlated = (0d0, 0d0)
    elseif (m_type .eq. 8) then
! Insert <ij>/(/ij/) which is not included by sborn()
      if (1d0 - y_ij_fks .lt. vtiny) then
        azifact = xij_aor
      else
        do i = 0, 3
          pi(i) = radiation%fks_momentum(i)
          pj(i) = p(i, selected_j_fks)
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
      imother_fks = min(selected_i_fks, selected_j_fks)
      call getaziangles(kernel_born(:, imother_fks), cphi_mother, &
                        sphi_mother)
      wgt1(2) = -(cphi_mother - ximag*sphi_mother)**2*wgt1(2)*azifact
      density_correlated = density_correlated* &
           (-(cphi_mother - ximag*sphi_mother)**2*azifact)
      amp_split_cnt(1:amp_split_size, 2, iord) = &
        -(cphi_mother - ximag*sphi_mother)**2 &
        *amp_split_cnt(1:amp_split_size, 2, iord)*azifact
    else
      write (*, *) 'FATAL ERROR in sborncol_fsr', i_type, j_type, &
                   selected_i_fks, selected_j_fks
      stop 1
    end if
    if (spin_density_fks_collection_enabled()) then
      call set_spin_density_reduced_from_born( &
           density_ordinary, density_correlated)
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
    type(factorized_radiation_state) :: radiation
    complex(kind=kind(0d0)) amp_split_cnt_local(amp_split_size, 2, nsplitorders)
    complex(kind=kind(0d0)) wgt1(2)
    complex(kind=8) :: density_ordinary, density_correlated
!
    amp_split_local(1:amp_split_size) = 0d0
    call load_kernel_radiation(real_event, radiation)

    p_born_used(:, :) = p_born(:, :)

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
    density_ordinary = cmplx(ap(1)*iden_comp, 0d0, kind=8)
    density_correlated = cmplx(Q(1)*iden_comp, 0d0, kind=8)
    wgt = 0d0
    iord = qcd_pos
    call evaluate_born_matrix(soft_counterevent, wgt_born)
    wgt1(1:2) = ans_cnt(1:2, iord)
    amp_split_cnt_local(1:amp_split_size, 1, iord) = amp_split_cnt(1:amp_split_size, 1, iord)
    amp_split_cnt_local(1:amp_split_size, 2, iord) = amp_split_cnt(1:amp_split_size, 2, iord)
    if (abs(m_type) .eq. 3) then
      Q(1) = 0d0
      wgt1(2) = cmplx(0d0, 0d0, kind=kind(0d0))
      density_correlated = (0d0, 0d0)
      amp_split_cnt_local(1:amp_split_size, 2, iord) = cmplx(0d0, 0d0, kind=kind(0d0))
    else
! Insert <ij>/(/ij/) which is not included by sborn()
      if (1d0 - y_ij_fks .lt. vtiny) then
        azifact = xij_aor
      else
        do i = 0, 3
          pi(i) = radiation%fks_momentum(i)
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
      density_correlated = density_correlated* &
           (-(cphi_mother + ximag*sphi_mother)**2*conjg(azifact))
      amp_split_cnt_local(1:amp_split_size, 2, iord) = &
        -(cphi_mother + ximag*sphi_mother)**2 &
        *amp_split_cnt_local(1:amp_split_size, 2, iord)*conjg(azifact)
    end if
    if (spin_density_fks_collection_enabled()) then
      call set_spin_density_reduced_from_born( &
           density_ordinary, density_correlated)
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


  subroutine sbornsoft(event_slot, pp, xi_i_fks, y_ij_fks, &
                       partonic_sqrt_shat, wgt)
    implicit none
!      include "fks.inc"
    integer, intent(in) :: event_slot
    integer m, n, visible_m, visible_n

    double precision softcontr, pp(0:3, nexternal), wgt, eik
    double precision xi_i_fks, y_ij_fks, partonic_sqrt_shat
    double precision wgt1
    integer i, j


    double precision zero, pmass(nexternal)
    double precision kernel_fks_momenta(0:3, &
         soft_counterevent:soft_collinear_counterevent)
    parameter(zero=0d0)


    if (has_nlo_decay()) then
      call sbornsoft_nlo_decay(event_slot, xi_i_fks, y_ij_fks, wgt)
      return
    end if

    call select_kernel_masses(pmass)
    call load_counterevent_fks_momenta(kernel_fks_momenta)
!
! Call the Born to be sure that 'CalculatedBorn' is done correctly. This
! should always be done before calling the color-correlated Borns,
! because of the caching of the diagrams.
!
    call evaluate_born_matrix(soft_counterevent, wgt1)
    if (spin_density_fks_collection_enabled()) then
      call reset_spin_density_reduced_matrix()
    end if
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
          call select_visible_color_pair(m, n, visible_m, visible_n)
          call evaluate_born_color_matrix( &
               soft_counterevent, visible_m, visible_n, wgt)
          if (wgt .ne. 0d0 .or. &
              spin_density_fks_collection_enabled()) then
            call eikonal_reduced( &
              pp, m, n, i_fks, j_fks, xi_i_fks, y_ij_fks, &
              kernel_fks_momenta, &
              pmass, partonic_sqrt_shat, eik)
            softcontr = softcontr + wgt*eik*iden_comp
            if (spin_density_fks_collection_enabled()) then
              ! The generated colour-linked Born returns its matrix before
              ! the explicit strong-coupling factor applied to the scalar
              ! WGT (see SBORN_SF_FACTORIZED).  Restore the same convention
              ! before accumulating the reduced density.
              call add_spin_density_reduced_color( &
                   cmplx(-2d0*eik*iden_comp*g**2, 0d0, kind=8))
            end if
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


  subroutine sbornsoft_nlo_decay(event_slot, xi_i_fks, y_ij_fks, wgt)
    integer, intent(in) :: event_slot
    double precision, intent(in) :: xi_i_fks, y_ij_fks
    double precision, intent(out) :: wgt
    integer :: first_position, second_position, m, n
    integer :: visible_m, visible_n, local_i, local_j
    integer :: partner_count, local_count
    double precision :: softcontr, eik, link_weight, born_weight
    double precision :: link_multiplier, decay_mass
    double precision :: local_momenta(0:3, nexternal)
    double precision :: local_masses(nexternal)
    double precision :: local_fks_momenta(0:3, 0:2)
    double precision, parameter :: zero = 0d0
    logical :: available

    local_count = nlo_decay_local_count( &
         nlo_decay_context_for_fks(nfksprocess))
    local_momenta = 0d0
    call fetch_factorized_kernel_momenta( &
         event_slot, active_factorized_block(), local_count, &
         local_momenta(:, 1:local_count), available)
    if (.not. available) then
      call fail_fks_singular_state( &
           'decay-local soft-kernel momenta are unavailable')
    end if
    call get_nlo_decay_mass_buffer(nfksprocess, local_masses)
    call load_counterevent_fks_momenta(local_fks_momenta)
    local_i = nlo_decay_fks_i(nfksprocess)
    local_j = nlo_decay_fks_j(nfksprocess)
    decay_mass = nlo_decay_parent_mass()

    ! Prime the generated Born cache before evaluating its color-linked forms.
    call evaluate_born_matrix(soft_counterevent, born_weight)
    if (spin_density_fks_collection_enabled()) then
      call reset_spin_density_reduced_matrix()
    end if
    amp_split(1:amp_split_size) = 0d0
    softcontr = 0d0
    partner_count = nlo_decay_partner_count(nfksprocess)
    do first_position = 1, partner_count
      do second_position = 1, first_position
        m = nlo_decay_partner_local(nfksprocess, first_position)
        n = nlo_decay_partner_local(nfksprocess, second_position)
        if ((m /= n .or. local_masses(m) /= zero) .and. &
            m /= local_i .and. n /= local_i) then
          call nlo_decay_map_color_link( &
               nfksprocess, m, n, visible_m, visible_n, link_multiplier)
          call evaluate_born_color_matrix( &
               soft_counterevent, visible_m, visible_n, link_weight)
          if (link_weight /= 0d0 .or. &
              spin_density_fks_collection_enabled()) then
            call eikonal_reduced( &
                 local_momenta, m, n, local_i, local_j, xi_i_fks, &
                 y_ij_fks, local_fks_momenta, local_masses, decay_mass, eik)
            softcontr = softcontr + &
                 link_multiplier*link_weight*eik*iden_comp
            if (spin_density_fks_collection_enabled()) then
              call add_spin_density_reduced_color(cmplx( &
                   -2d0*link_multiplier*eik*iden_comp*g**2, &
                   0d0, kind=8))
            end if
            amp_split(1:amp_split_size) = &
                 amp_split(1:amp_split_size) - &
                 2d0*link_multiplier*eik* &
                 amp_split_soft(1:amp_split_size)*iden_comp
          end if
        end if
      end do
    end do

    ! Match the sign and normalization convention of the standalone decay.
    wgt = -2d0*softcontr
  end subroutine sbornsoft_nlo_decay


  subroutine sreal_deg(event_slot, xi_i_fks, collrem_xi, collrem_lxi)
    use extra_weights
    implicit none
    integer event_slot, iord, iap
    double precision collrem_xi, collrem_lxi
    double precision xi_i_fks
    double precision collrem_xi_tmp, collrem_lxi_tmp

    double precision wgt_born

    double precision shattmp, oo2pi, z, t, ap(2), apprime(2), xnorm
    double precision kernel_momenta(0:3, nexternal)
    double precision kernel_masses(nexternal)
    integer :: kernel_initial_count

! Colour representations of i_fks, j_fks and the FKS mother
    complex(kind=kind(0d0)) wgt1(2)

    double precision one, pi
    parameter(one=1.d0)
    parameter(pi=3.1415926535897932385d0)

    double precision amp_split_collrem_xi(amp_split_size), amp_split_collrem_lxi(amp_split_size)
    double precision prefact_xi
    double precision subtraction_shat
    complex(kind=8) :: density_multipliers(3)

    amp_split_collrem_xi(1:amp_split_size) = 0d0
    amp_split_collrem_lxi(1:amp_split_size) = 0d0
    amp_split_wgtdegrem_xi(1:amp_split_size) = 0d0
    amp_split_wgtdegrem_lxi(1:amp_split_size) = 0d0
    amp_split_wgtdegrem_muF(1:amp_split_size) = 0d0
    if (spin_density_fks_collection_enabled()) then
      call reset_spin_density_degenerate_matrix()
    end if

    subtraction_shat = fks_subtraction_shat(event_slot)
    call select_kernel_event(event_slot, kernel_momenta, kernel_masses)

    if (j_fks .gt. nincoming) then
! Do not include this contribution for final-state branchings
      collrem_xi = 0.d0
      collrem_lxi = 0.d0
      return
    end if

    if (kernel_momenta(0, 1) .le. 0.d0) then
! Unphysical kinematics: set matrix elements equal to zero
      write (*, *) "No born momenta in sreal_deg"
      collrem_xi = 0.d0
      collrem_lxi = 0.d0
      return
    end if

! Check that the requested event slot matches the supplied momenta.
    kernel_initial_count = merge(1, nincoming, has_nlo_decay())
    if (kernel_initial_count .eq. 2) then
      shattmp = 2d0*dot(kernel_momenta(0, 1), kernel_momenta(0, 2))
    else
      shattmp = kernel_momenta(0, 1)**2
    end if
    if (abs(shattmp/subtraction_shat - 1.d0) .gt. 1.d-5) then
      write (*, *) 'Error in sreal: inconsistent shat'
      write (*, *) shattmp, subtraction_shat
      stop
    end if

! A factor gS^2 is included in the Altarelli-Parisi kernels
    oo2pi = one/(8d0*PI**2)

    z = 1d0 - xi_i_fks
    t = one
    call AP_reduced(m_type, i_type, t, z, g, ap)
    call AP_reduced_prime(m_type, i_type, t, z, g, apprime)

    collrem_xi = 0.d0
    collrem_lxi = 0.d0
    calculatedborn = .false.
    iord = qcd_pos
    iap = 1
    call evaluate_born_matrix(soft_counterevent, wgt_born)
    wgt1(1) = ans_cnt(1, iord)
    wgt1(2) = ans_cnt(2, iord)

    collrem_xi_tmp = ap(iap)*log( &
                       subtraction_shat*delta_used/ &
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
                   subtraction_shat*delta_used/(2*QES2)) - &
                 apprime(iap)
    amp_split_wgtdegrem_xi(1:amp_split_size) = &
      amp_split_wgtdegrem_xi(1:amp_split_size) + &
      oo2pi*dble(amp_split_cnt(1:amp_split_size, 1, iord)) &
      *prefact_xi*xnorm
    amp_split_wgtdegrem_lxi(1:amp_split_size) = amp_split_collrem_lxi(1:amp_split_size)
    amp_split_wgtdegrem_muF(1:amp_split_size) = &
      amp_split_wgtdegrem_muF(1:amp_split_size) - &
      oo2pi*dble(amp_split_cnt(1:amp_split_size, 1, iord))*ap(iap)*xnorm
    if (spin_density_fks_collection_enabled()) then
      density_multipliers(1) = cmplx( &
           oo2pi*prefact_xi*xnorm, 0d0, kind=8)
      density_multipliers(2) = cmplx( &
           oo2pi*collrem_lxi_tmp*xnorm, 0d0, kind=8)
      density_multipliers(3) = cmplx( &
           -oo2pi*ap(iap)*xnorm, 0d0, kind=8)
      call set_spin_density_degenerate_from_born(density_multipliers)
    end if
    calculatedborn = .false.

    return
  end subroutine sreal_deg



  subroutine bornsoftvirtual(event_slot, bsv_wgt, virt_wgt, born_wgt)
    use extra_weights
    use mint_module
    implicit none
    real :: tBefore, tAfter
!      include "fks.inc"
    double precision bsv_wgt, born_wgt, avv_wgt
    double precision wgt1
    double precision Q, Ej, wgt, contr, eikIreg
    double precision aso2pi
    double precision shattmp
    integer event_slot, i, j, aj, m, n, k, virtual_grid
    integer kernel_count, visible_m, visible_n
    integer kernel_i, kernel_initial_count, partner_count
    integer kernel_particle_type(nexternal)



    double precision pi
    parameter(pi=3.1415926535897932385d0)

    double precision c_used, gamma_used, gammap_used
    double precision virt_wgt



    double precision pmass(nexternal), zero
    parameter(zero=0d0)
    logical firsttime
    data firsttime/.true./
    logical need_color_links_used
    data need_color_links_used/.false./
    double precision oneo8pi2
    parameter(oneo8pi2=1d0/(8d0*pi**2))
    integer nFKSprocess_save, nFKSprocess_col, scan_first, scan_last
    data nFKSprocess_col/0/
    double precision bsv_wgt_mufoqes, bsv_wgt_mufomur
    double precision contr_mufoqes, contr_mufomur
! to keep track of the various split orders
    integer iamp
    integer orders(nsplitorders), production_born_qcd_order
    double precision amp_split_born(amp_split_size)
    double precision amp_split_bsv(amp_split_size)
    double precision kernel_momenta(0:3, nexternal)
    double precision kernel_shat, kernel_sqrt_shat, link_multiplier
    complex(kind=8) :: density_coefficients(3)
    double precision :: density_muf, density_mur, density_factor
    double precision :: density_virtual_average, density_sampling_fraction
    double precision :: virtual_sampling_fraction
    logical, external :: sdm_virtual_uses_analytic_provider
    type(factorized_radiation_state) :: radiation
    if (has_nlo_contribution_bundle()) then
      need_color_links_used = .false.
      nFKSprocess_col = 0
      nFKSprocess_save = nFKSprocess
      scan_first = active_contribution_fks_first()
      scan_last = active_contribution_fks_last()
      do nFKSprocess = scan_first, scan_last
        call fks_inc_chooser()
        need_color_links_used = need_color_links_used .or. need_color_links
        if (need_color_links .and. nFKSprocess_col .eq. 0) then
          nFKSprocess_col = nFKSprocess
        end if
      end do
      nFKSprocess = nFKSprocess_save
      call fks_inc_chooser()
    else if (firsttime) then
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
    call select_kernel_event(event_slot, kernel_momenta, pmass, &
                             kernel_count)
    call load_kernel_radiation(event_slot, radiation)
    kernel_sqrt_shat = radiation%sqrt_shat
    kernel_shat = radiation%shat
    call select_kernel_properties(kernel_particle_type, kernel_i, &
                                  kernel_initial_count)

    aso2pi = g**2/(8*pi**2)

    amp_split_bsv(1:amp_split_size) = 0d0
    amp_split_virt(1:amp_split_size) = 0d0
    amp_split_avv(1:amp_split_size) = 0d0
    virtual_sampling_fraction = virtual_fraction(ichan)
    if (sdm_virtual_uses_analytic_provider( &
        active_nlo_contribution())) virtual_sampling_fraction = 1d0
    if (spin_density_fks_collection_enabled()) then
      call reset_spin_density_integrated_matrix()
      call reset_spin_density_virtual_matrix()
    end if

    if (.not. need_color_links_used) then
! just return 0
      bsv_wgt = 0d0
      virt_wgt = 0d0
      born_wgt = 0d0
      goto 999
    end if

! Check that the requested event slot matches the supplied momenta.
    if (kernel_initial_count .eq. 2) then
      shattmp = 2d0*dot(kernel_momenta(0, 1), &
                        kernel_momenta(0, 2))
    else
      shattmp = kernel_momenta(0, 1)**2
    end if
    if (abs(shattmp/kernel_shat - 1.d0) .gt. 1.d-5) then
      write (*, *) 'Error in bornsoftvirtual: inconsistent shat'
      write (*, *) shattmp, kernel_shat
      stop
    end if

    call evaluate_born_matrix(soft_counterevent, wgt1)
    if (spin_density_fks_collection_enabled()) then
      density_coefficients = (0d0, 0d0)
      density_coefficients(1) = (1d0, 0d0)
      call add_spin_density_integrated_born(density_coefficients)
    end if

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
    do i = 1, kernel_count
      if (i .ne. kernel_i .and. pmass(i) .eq. ZERO) then
! set the colour factors according to the
! type of the leg
        if (kernel_particle_type(i) .eq. 8) then
          aj = 0
        elseif (abs(kernel_particle_type(i)) .eq. 3) then
          aj = 1
        else
          aj = -1
        end if
        Ej = kernel_momenta(0, i)

        if (aj .eq. -1) cycle
        c_used = c(aj)
        gamma_used = gamma(aj)
        gammap_used = gammap(aj)

        if (i .gt. kernel_initial_count) then
! Q terms for final state partons
          if (abrv .ne. 'virt') then
! 1+2+3+4
            Q = Q + gammap_used &
                - dlog(kernel_shat*deltaO/2d0/QES2) &
                *(gamma_used &
                  - 2d0*c_used*dlog( &
                      2d0*Ej/xicut_used/kernel_sqrt_shat)) &
                + 2d0*c_used*(dlog( &
                      2d0*Ej/kernel_sqrt_shat)**2 &
                              - dlog(xicut_used)**2) &
                - 2d0*gamma_used*dlog( &
                    2d0*Ej/kernel_sqrt_shat)
          else
            write (*, *) 'Error in bornsoftvirtual'
            write (*, *) 'abrv in Q:', abrv
            stop
          end if

        else
! Q terms for initial state partons
          if (abrv .ne. 'virt') then
! 1+2+3+4
            Q = Q - dlog(q2fact(i)/QES2) &
                *(gamma_used + 2d0*c_used*dlog(xicut_used))
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
    amp_split_bsv(1:amp_split_size) = amp_split_bsv(1:amp_split_size) + &
      aso2pi*Q*dble(amp_split_cnt(1:amp_split_size, 1, qcd_pos))
    if (spin_density_fks_collection_enabled()) then
      density_coefficients = (0d0, 0d0)
      density_coefficients(1) = cmplx(aso2pi*Q, 0d0, kind=8)
      call add_spin_density_integrated_born(density_coefficients)
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
      call evaluate_born_matrix(soft_counterevent, wgt1)
      contr = 0d0
      partner_count = selected_partner_count(kernel_i)
      do i = 1, partner_count
        do j = 1, i
          m = selected_partner(kernel_i, i)
          n = selected_partner(kernel_i, j)
          if ((m .ne. n .or. (m .eq. n .and. &
              pmass(m) .ne. ZERO)) .and. n .ne. kernel_i .and. &
              m .ne. kernel_i) then
! To be sure that color-correlated Borns work well, we need to have
! *always* a call to sborn(p_born,wgt) just before. This is okay,
! because there is a call above in this subroutine
! wgt includes the gs/w^2
            call select_visible_color_pair(m, n, visible_m, visible_n, &
                                           link_multiplier)
            call evaluate_born_color_matrix( &
                 soft_counterevent, visible_m, visible_n, wgt)
            if (wgt .ne. 0d0 .or. &
                spin_density_fks_collection_enabled()) then
              call eikonal_Ireg(kernel_momenta, m, n, xicut_used, pmass, &
                                 kernel_shat, qes2, abrv, &
                                 eikIreg)
              contr = contr + link_multiplier*wgt*eikIreg
              if (spin_density_fks_collection_enabled()) then
                density_coefficients = (0d0, 0d0)
                density_coefficients(1) = cmplx( &
                     -2d0*link_multiplier*eikIreg*oneo8pi2*g**2, &
                     0d0, kind=8)
                call add_spin_density_integrated_color( &
                     density_coefficients)
              end if
              do k = 1, amp_split_size
                amp_split_bsv(k) = amp_split_bsv(k) - &
                  2d0*link_multiplier*eikIreg*oneo8pi2*amp_split_soft(k)
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

    call evaluate_born_matrix(soft_counterevent, wgt1)
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

    if (active_contribution_has_virtual() .and. &
        ((random_unit_interval(iconfig) .le. virtual_sampling_fraction &
         .and. abrv(1:3) .ne. 'nov') .or. &
        abrv(1:4) .eq. 'virt')) then
      call cpu_time(tBefore)
      call evaluate_virtual_matrix( &
           soft_counterevent, born_wgt, virt_wgt)
      do iamp = 1, amp_split_size
        amp_split_virt(iamp) = amp_split_finite_ML(iamp)
      end do
      virtual_over_born = virt_wgt/born_wgt
      virt_wgt = 0d0
      do iamp = 1, amp_split_size
        if (amp_split_virt(iamp) .eq. 0d0) cycle
        virtual_grid = active_virtual_grid_index(iamp, amp_split_size)
        if (virtual_grid == 0) then
          write (*,*) 'ERROR: a virtual split order has no bundle grid'
          stop 1
        end if
        if (use_poly_virtual) then
          amp_split_virt(iamp) = amp_split_virt(iamp) - &
            polyfit(virtual_grid)*amp_split_born_for_virt(iamp)
        else
          amp_split_virt(iamp) = amp_split_virt(iamp) - &
            average_virtual(virtual_grid, ichan)* &
            amp_split_born_for_virt(iamp)
        end if
        virt_wgt = virt_wgt + amp_split_virt(iamp)
      end do
      if (abrv .ne. 'virt') then
        virt_wgt = virt_wgt/virtual_sampling_fraction
        do iamp = 1, amp_split_size
          amp_split_virt(iamp) = &
          amp_split_virt(iamp)/virtual_sampling_fraction
        end do
      end if
      if (spin_density_fks_collection_enabled() .and. &
          spin_density_virtual_matrix_available()) then
        virtual_grid = active_virtual_grid_index( &
             spin_density_nlo_amp_position(), amp_split_size)
        if (virtual_grid == 0) then
          write (*,*) 'ERROR: the density virtual has no bundle grid'
          stop 1
        end if
        if (use_poly_virtual) then
          density_virtual_average = polyfit(virtual_grid)
        else
          density_virtual_average = average_virtual(virtual_grid, ichan)
        end if
        density_sampling_fraction = 1d0
        if (abrv .ne. 'virt') then
          density_sampling_fraction = virtual_sampling_fraction
        end if
        call reduce_spin_density_virtual_matrix( &
             density_virtual_average, density_sampling_fraction)
      end if
      call cpu_time(tAfter)
      tOLP = tOLP + (tAfter - tBefore)
    end if
    if (abrv(1:4) .ne. 'virt' .and. &
        active_contribution_has_virtual()) then
      if (use_poly_virtual) then
        avv_wgt = 0d0
        do iamp = 1, amp_split_size
          if (amp_split_born_for_virt(iamp) .eq. 0d0) cycle
          virtual_grid = active_virtual_grid_index(iamp, amp_split_size)
          if (virtual_grid == 0) then
            write (*,*) 'ERROR: an averaged virtual has no bundle grid'
            stop 1
          end if
          amp_split_avv(iamp) = polyfit(virtual_grid)* &
               amp_split_born_for_virt(iamp)
          avv_wgt = avv_wgt + amp_split_avv(iamp)
        end do
      else
        avv_wgt = 0d0
        do iamp = 1, amp_split_size
          if (amp_split_born_for_virt(iamp) .eq. 0d0) cycle
          virtual_grid = active_virtual_grid_index(iamp, amp_split_size)
          if (virtual_grid == 0) then
            write (*,*) 'ERROR: an averaged virtual has no bundle grid'
            stop 1
          end if
          amp_split_avv(iamp) = average_virtual(virtual_grid, ichan)* &
               amp_split_born_for_virt(iamp)
          avv_wgt = avv_wgt + amp_split_avv(iamp)
        end do
      end if
    end if

! eq.(MadFKS.C.13)
    if (abrv .ne. 'virt') then
! this is to update the amp_split array
      call evaluate_born_matrix(soft_counterevent, wgt1)
      bsv_wgt_mufoqes = 0d0
      do iamp = 1, amp_split_size
        if (dble(amp_split_cnt(iamp, 1, qcd_pos)) .eq. 0d0) cycle
        call amp_split_pos_to_orders(iamp, orders)
        production_born_qcd_order = &
             corrected_born_qcd_squared_order(orders(qcd_pos))
        contr_mufoqes = pi*beta0*dble(production_born_qcd_order) &
             *log(q2fact(1)/QES2)*aso2pi &
             *dble(amp_split_cnt(iamp, 1, qcd_pos))
        amp_split_bsv(iamp) = amp_split_bsv(iamp) + contr_mufoqes
        bsv_wgt_mufoqes = bsv_wgt_mufoqes + contr_mufoqes
        if (spin_density_fks_collection_enabled() .and. &
            iamp == spin_density_nlo_amp_position()) then
          density_factor = pi*beta0*dble(production_born_qcd_order)* &
               log(q2fact(1)/QES2)*aso2pi
          density_coefficients = (0d0, 0d0)
          density_coefficients(1) = cmplx(density_factor, 0d0, kind=8)
          call add_spin_density_integrated_born(density_coefficients)
        end if
      end do
      bsv_wgt = bsv_wgt + bsv_wgt_mufoqes
    end if

!  eq.(MadFKS.C.14)
    if (abrv(1:2) .ne. 'vi') then
      bsv_wgt_mufomur = 0d0
      do iamp = 1, amp_split_size
        if (dble(amp_split_cnt(iamp, 1, qcd_pos)) .eq. 0d0) cycle
        call amp_split_pos_to_orders(iamp, orders)
        production_born_qcd_order = &
             corrected_born_qcd_squared_order(orders(qcd_pos))
        contr_mufomur = -pi*beta0*dble(production_born_qcd_order) &
             *log(q2fact(1)/scale**2)*aso2pi &
             *dble(amp_split_cnt(iamp, 1, qcd_pos))
        amp_split_bsv(iamp) = amp_split_bsv(iamp) + contr_mufomur
        bsv_wgt_mufomur = bsv_wgt_mufomur + contr_mufomur
        if (spin_density_fks_collection_enabled() .and. &
            iamp == spin_density_nlo_amp_position()) then
          density_factor = -pi*beta0*dble(production_born_qcd_order)* &
               log(q2fact(1)/scale**2)*aso2pi
          density_coefficients = (0d0, 0d0)
          density_coefficients(1) = cmplx(density_factor, 0d0, kind=8)
          call add_spin_density_integrated_born(density_coefficients)
        end if
      end do
      bsv_wgt = bsv_wgt + bsv_wgt_mufomur
    end if

549 continue

    amp_split_wgtnstmp(1:amp_split_size) = 0d0
    amp_split_wgtwnstmpmuf(1:amp_split_size) = 0d0
    amp_split_wgtwnstmpmur(1:amp_split_size) = 0d0

    if (abrv .ne. 'born' .and. abrv .ne. 'grid') then
      call evaluate_born_matrix(soft_counterevent, wgt1)
      density_muf = 0d0
      density_mur = 0d0
      if (abrv(1:2) .ne. 'vi') then
        do i = 1, kernel_initial_count
          if (pmass(i) .ne. zero) cycle
          if (kernel_particle_type(i) .eq. 8) then
            aj = 0
          elseif (abs(kernel_particle_type(i)) .eq. 3) then
            aj = 1
          else
            aj = -1
          end if
          if (aj .eq. -1) cycle
          c_used = c(aj)
          gamma_used = gamma(aj)
          density_muf = density_muf - &
               (gamma_used + 2d0*c_used*dlog(xicut_used))*aso2pi
          do iamp = 1, amp_split_size
            if (dble(amp_split_cnt(iamp, 1, qcd_pos)) .eq. 0d0) cycle
            amp_split_wgtwnstmpmuf(iamp) = &
              amp_split_wgtwnstmpmuf(iamp) &
              - (gamma_used + 2d0*c_used*dlog(xicut_used)) &
              *dble(amp_split_cnt(iamp, 1, qcd_pos))*aso2pi
          end do
        end do
        do iamp = 1, amp_split_size
          if (dble(amp_split_cnt(iamp, 1, qcd_pos)) .eq. 0d0) cycle
          call amp_split_pos_to_orders(iamp, orders)
          production_born_qcd_order = &
               corrected_born_qcd_squared_order(orders(qcd_pos))
          amp_split_wgtwnstmpmur(iamp) = dble(amp_split_cnt(iamp, 1, qcd_pos)) &
               *pi*beta0*dble(production_born_qcd_order)*aso2pi
          if (iamp == spin_density_nlo_amp_position()) then
            density_mur = pi*beta0*dble(production_born_qcd_order)*aso2pi
          end if
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
      if (spin_density_fks_collection_enabled()) then
        density_coefficients(1) = cmplx( &
             -1d0 - log(q2fact(1)/QES2)*density_muf - &
             log(scale**2/QES2)*density_mur, 0d0, kind=8)
        density_coefficients(2) = cmplx(density_mur, 0d0, kind=8)
        density_coefficients(3) = cmplx(density_muf, 0d0, kind=8)
        call add_spin_density_integrated_born(density_coefficients)
      end if
    end if

    amp_split(1:amp_split_size) = amp_split_bsv(1:amp_split_size)

    if (abrv(1:2) .eq. 'vi') then
      bsv_wgt = bsv_wgt - born_wgt
      born_wgt = 0d0
    end if

999 continue
    return
  end subroutine bornsoftvirtual


  subroutine getpoles(p, double, single, split_poles)
! Returns the residues of double and single poles in the FKS convention
! of eqs. (B.1) and (B.2).
    implicit none
    double precision p(0:3, nexternal), double, single
    double precision, optional, intent(out) :: split_poles(amp_split_size, 2)
    double precision wgt1
    double precision born, wgt, kikj, vij, aso2pi
    double precision contr1, contr2
    integer aj, i, j, m, n, k
    integer kernel_count, visible_m, visible_n
    integer kernel_i, kernel_initial_count, partner_count
    integer kernel_particle_type(nexternal)
    double precision pmass(nexternal), zero, pi
    parameter(pi=3.1415926535897932385d0)
    parameter(zero=0d0)
    double precision oneo8pi2
    parameter(oneo8pi2=1d0/(8d0*pi**2))
    integer nFKSprocess_save, nFKSprocess_col, scan_first, scan_last
    logical need_color_links_used
    double precision soft_fact
    double precision link_multiplier
    double precision kernel_p(0:3, nexternal)
    double precision kernel_born(0:3, nexternal - 1)
    nFKSprocess_col = 0

    need_color_links_used = .false.

! Check whether any real-emission configuration needs colour links.
    nFKSprocess_save = nFKSprocess
    if (has_nlo_contribution_bundle()) then
      scan_first = active_contribution_fks_first()
      scan_last = active_contribution_fks_last()
    else
      scan_first = 1
      scan_last = FKS_configs
    end if
    do nFKSprocess = scan_first, scan_last
      call fks_inc_chooser()
      need_color_links_used = need_color_links_used .or. need_color_links
      if (need_color_links .and. nFKSprocess_col .eq. 0) nFKSprocess_col = nFKSprocess
    end do
    nFKSprocess = nFKSprocess_save
    call fks_inc_chooser()
    if (has_nlo_decay()) then
      call get_nlo_decay_born_kernel( &
           nfksprocess, kernel_p, pmass, kernel_count)
    else
      call select_kernel_masses(pmass, kernel_count)
    end if
    call select_kernel_properties(kernel_particle_type, kernel_i, &
                                  kernel_initial_count)
    if (has_decay_chains()) then
      kernel_p = 0d0
      call select_kernel_born(kernel_born)
      kernel_p(:, 1:nexternal - 1) = kernel_born
    else if (.not. has_nlo_decay()) then
      kernel_p = p
    end if

    double = 0.d0
    single = 0.d0
! reset the amp_split_poles_FKS
    do i = 1, amp_split_size
      amp_split_poles_FKS(i, 1) = 0d0
      amp_split_poles_FKS(i, 2) = 0d0
    end do
    aso2pi = g**2/(8d0*pi**2)
    call evaluate_born_matrix(soft_counterevent, wgt1)
! QCD Born terms
    contr1 = 0d0
    contr2 = 0d0
    born = dble(ans_cnt(1, qcd_pos))
    do i = 1, kernel_count
      if (i .ne. kernel_i .and. kernel_particle_type(i) .ne. 1) then
        if (kernel_particle_type(i) .eq. 8) then
          aj = 0
        elseif (abs(kernel_particle_type(i)) .eq. 3) then
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
      call evaluate_born_matrix(soft_counterevent, wgt1)

      contr1 = 0d0
      partner_count = selected_partner_count(kernel_i)
      do i = 1, partner_count
        do j = 1, i
          m = selected_partner(kernel_i, i)
          n = selected_partner(kernel_i, j)
          if (m .ne. n .and. n .ne. kernel_i .and. &
              m .ne. kernel_i) then
! wgt includes the gs/w^2 factor
            call select_visible_color_pair(m, n, visible_m, visible_n, &
                                           link_multiplier)
            call evaluate_born_color_matrix( &
                 soft_counterevent, visible_m, visible_n, wgt)
! The factor -2 compensate for that missing in sborn_sf
            wgt = -2d0*link_multiplier*wgt
            if (wgt .ne. 0.d0) then
              if (pmass(m) .eq. zero .and. pmass(n) .eq. zero) then
                kikj = dot(kernel_p(0, n), kernel_p(0, m))
                soft_fact = dlog(2d0*kikj/QES2)
              elseif (pmass(m) .ne. zero .and. pmass(n) .eq. zero) then
                kikj = dot(kernel_p(0, n), kernel_p(0, m))
                soft_fact = -0.5d0*dlog(pmass(m)**2/QES2) + dlog(2d0*kikj/QES2)
              elseif (pmass(m) .eq. zero .and. pmass(n) .ne. zero) then
                kikj = dot(kernel_p(0, n), kernel_p(0, m))
                soft_fact = -0.5d0*dlog(pmass(n)**2/QES2) + dlog(2d0*kikj/QES2)
              elseif (pmass(m) .ne. zero .and. pmass(n) .ne. zero) then
                kikj = dot(kernel_p(0, n), kernel_p(0, m))
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
                amp_split_poles_FKS(k, 1) = &
                     amp_split_poles_FKS(k, 1) + &
                     link_multiplier*amp_split_soft(k)*(-2d0)* &
                     soft_fact*oneo8pi2
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
    integer real_particle_count, born_particle_count, born_pdg(nexternal)
    integer decay_context, decay_born_context, decay_i, decay_j
    integer decay_leg, decay_other



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
    gamma(0) = (11d0*CA - 2d0*active_flavours)/6d0
    gamma(1) = CF*3d0/2d0
    gammap(0) = (67d0/9d0 - 2d0*PI**2/3d0)*CA - &
                23d0/18d0*active_flavours
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
      if (has_nlo_decay()) then
        decay_context = nlo_decay_context_for_fks(nfksprocess)
        decay_born_context = nlo_decay_born_context()
        decay_i = nlo_decay_fks_i(nfksprocess)
        decay_j = nlo_decay_fks_j(nfksprocess)
        i_fks_pdg = nlo_decay_local_pdg(decay_context, decay_i)
        j_fks_pdg = nlo_decay_local_pdg(decay_context, decay_j)
      else
        i_fks_pdg = pdg_type(i_fks)
        j_fks_pdg = pdg_type(j_fks)
      end if
      if (has_decay_chains() .and. .not. has_nlo_decay()) then
        real_particle_count = active_core_count(nfksprocess)
        born_particle_count = context_core_count(born_context())
        born_pdg = 0
        do i = 1, born_particle_count
          born_pdg(i) = core_leg_pdg(born_context(), i)
        end do
      else
        real_particle_count = nexternal
        born_particle_count = nexternal - 1
      end if

      fac_i_FKS(nFKSprocess) = 0
      fac_j_FKS(nFKSprocess) = 0
      if (has_nlo_decay()) then
        ! The NWA labels particles by their production/decay history.  Only
        ! permutations inside the corrected decay belong to this FKS family;
        ! equal-PDG legs from production or another decay are spectators.
        do decay_leg = 1, nlo_decay_local_count(decay_context)
          if (.not. nlo_decay_local_is_final(decay_context, decay_leg)) cycle
          if (i_fks_pdg .eq. &
              nlo_decay_local_pdg(decay_context, decay_leg)) then
            fac_i_FKS(nFKSprocess) = fac_i_FKS(nFKSprocess) + 1
          end if
          if (j_fks_pdg .eq. &
              nlo_decay_local_pdg(decay_context, decay_leg)) then
            fac_j_FKS(nFKSprocess) = fac_j_FKS(nFKSprocess) + 1
          end if
        end do
      else
        do i = nincoming + 1, real_particle_count
          if (i_fks_pdg .eq. pdg_type(i)) then
            fac_i_FKS(nFKSprocess) = fac_i_FKS(nFKSprocess) + 1
          end if
          if (j_fks_pdg .eq. pdg_type(i)) then
            fac_j_FKS(nFKSprocess) = fac_j_FKS(nFKSprocess) + 1
          end if
        end do
      end if
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
      if (has_nlo_decay()) then
        do decay_leg = 1, nlo_decay_local_count(decay_context)
          if (.not. nlo_decay_local_is_final(decay_context, decay_leg)) cycle
          if (nlo_decay_local_pdg(decay_context, decay_leg) .eq. 21) then
            ngluons_FKS(nFKSprocess) = ngluons_FKS(nFKSprocess) + 1
          end if
        end do
      else
        do i = nincoming + 1, real_particle_count
          if (pdg_type(i) .eq. 21) then
            ngluons_FKS(nFKSprocess) = ngluons_FKS(nFKSprocess) + 1
          end if
        end do
      end if

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
      if (has_nlo_decay()) then
        do decay_leg = 2, nlo_decay_local_count(decay_context)
          if (.not. nlo_decay_local_is_final(decay_context, decay_leg)) cycle
          do decay_other = 1, decay_leg - 1
            if (.not. nlo_decay_local_is_final( &
                decay_context, decay_other)) cycle
            if (nlo_decay_local_pdg(decay_context, decay_other) .eq. &
                nlo_decay_local_pdg(decay_context, decay_leg)) then
              iden(decay_other) = iden(decay_other) + 1
              iden_real_FKS(nFKSprocess) = &
                   iden_real_FKS(nFKSprocess)*iden(decay_other)
              exit
            end if
          end do
        end do
      else
        do i = nincoming + 2, real_particle_count
          do j = nincoming + 1, i - 1
            if (pdg_type(j) .eq. pdg_type(i)) then
              iden(j) = iden(j) + 1
              iden_real_FKS(nFKSprocess) = &
                   iden_real_FKS(nFKSprocess)*iden(j)
              exit
            end if
          end do
        end do
      end if
! Compute the identical particle symmetry factor that is in the
! Born matrix elements.
      iden_born_FKS(nFKSprocess) = 1
      call weight_lines_allocated(nexternal, max_contr, max_wgt, max_iproc)
      call set_pdg_impl(0, nFKSprocess, idup)
      if (.not. has_decay_chains() .and. .not. has_nlo_decay()) then
        born_pdg = pdg_uborn(:, 0)
      end if
      do i = 1, nexternal
        iden(i) = 1
      end do
      if (has_nlo_decay()) then
        do decay_leg = 2, nlo_decay_local_count(decay_born_context)
          if (.not. nlo_decay_local_is_final( &
              decay_born_context, decay_leg)) cycle
          do decay_other = 1, decay_leg - 1
            if (.not. nlo_decay_local_is_final( &
                decay_born_context, decay_other)) cycle
            if (nlo_decay_local_pdg(decay_born_context, decay_other) .eq. &
                nlo_decay_local_pdg(decay_born_context, decay_leg)) then
              iden(decay_other) = iden(decay_other) + 1
              iden_born_FKS(nFKSprocess) = &
                   iden_born_FKS(nFKSprocess)*iden(decay_other)
              exit
            end if
          end do
        end do
      else
        do i = nincoming + 2, born_particle_count
          do j = nincoming + 1, i - 1
            if (born_pdg(j) .eq. born_pdg(i)) then
              iden(j) = iden(j) + 1
              iden_born_FKS(nFKSprocess) = &
                   iden_born_FKS(nFKSprocess)*iden(j)
              exit
            end if
          end do
        end do
      end if
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
