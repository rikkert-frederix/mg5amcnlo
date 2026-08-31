module driver_mintfo_module
  use run_printout_module, only: write_run_summary
  use extra_weights, only: doreweight
  use mint_module, only: maxchannels, n_ave_virt, average_virtual, &
                         virtual_fraction, min_virt_fraction_mint, n_ord_virt, ncalls0, &
                         itmax, imode, ndim, ndimmax, nintegrals, nchans, &
                         iconfig, ichan, &
                         iconfigs, accuracy, wgt_mult, new_point, pass_cuts_check, &
                         virt_wgt_mint, born_wgt_mint, mint, &
                         first_bundle_component_integral, &
                         configure_decay_folding, &
                         decay_fold_group_active, &
                         decay_fold_reuses_production, &
                         record_decay_fold_production_cache, &
                         record_decay_fold_conditional_sector
  use mint_module, only: ans_result => ans, unc_result => unc
  use FKSParams, only: min_virt_fraction, virt_fraction, FKSParamReader
  use weight_lines, only: icontr, iwgt, deallocate_weight_lines
  use process_dimensions, only: nexternal, nincoming, fks_configs, &
                                amp_split_size, lmaxconfigs, &
                                validate_process_and_born_dimensions
  use fks_metadata, only: fks_i_d, pdg_type_d, validate_fks_metadata
  use fks_channel_map, only: fks_channel_count, &
                             fks_channel_configuration, &
                             contribution_channel_count, &
                             contribution_channel_configuration, &
                             print_fks_channel_map, get_born_fks_process
  use fks_random_module, only: random_unit_interval
  use run_state, only: lpp, fixed_fac_scale, muf1_over_ref, &
                       muf2_over_ref, muf1_ref_fixed, muf2_ref_fixed, &
                       do_rwgt_scale, do_rwgt_decay_scale, do_rwgt_pdf
  use genps_fks, only: generate_momenta, generate_momenta_reusing_born
  use decay_chain_metadata, only: real_phase_space_dimension, &
                                  decay_random_dimension
  use decay_chain_kinematics, only: decay_variable_start
  use fnlo_scale_variations, only: configure_fnlo_scale_variations
  use nlo_contribution_bundle, only: has_nlo_contribution_bundle, &
       nlo_contribution_count, contribution_representative_fks, &
       contribution_for_fks, contribution_is_nlo_decay, &
       factorized_integration_dimension, &
       factorized_radiation_start, &
       nlo_virtual_grid_count, bundle_component_count, &
       bundle_component_label, bundle_nlo_component
  use setscales_module, only: set_alphas
  use split_orders, only: check_amp_split
  use cuts_module, only: passcuts, passcuts_multiplicative
  use mc_integer_module, only: get_mc_integer, fill_mc_integer, &
       mc_integer_interval_volume
  use timing_state, only: reset_timing_state, tBorn, tIS, tReal, &
                          tCount, tf_nb, tf_all, t_as, tr_s, tr_pdf, t_plot, &
                          t_cuts, t_isum, tOLP, tGenPS, t_coupl
  use fks_contributions_module, only: compute_prefactors_nbody, &
       compute_prefactors_n1body, include_multichannel_enhance, &
       compute_born, compute_decay_width_counterterm, &
       record_multiplicative_spectator_born_densities, &
       compute_nbody_noborn, compute_soft_counter_term, &
       compute_soft_collinear_ct_impl, compute_collinear_counter_term, &
       compute_real_emission
  use fks_weights_module, only: include_pdf_and_alphas, reweight_scale, &
       reweight_pdf, get_wgt_no_nbody, fill_plots, fill_mint_function, &
       begin_bundle_virtual_tricks, finish_bundle_virtual_tricks, &
       reset_luminosity_cache
  use fks_singular_module, only: fill_configurations_common, setfksfactor
  use madfks_plot_module, only: topout_impl, outfun_multiplicative_impl
  use decay_chain_parameters, only: multiplicative_nlo_enabled, &
       use_decayed_production_ren_scale_momenta
  use spin_density_fks_matrices, only: &
       set_spin_density_fks_collection, reset_spin_density_fks_matrices
  use spin_density_weight_lines, only: clear_spin_density_weight_lines, &
       aggregate_spin_density_weight_lines
  use spin_density_matrix_results, only: spin_density_bornlike_branch, &
       spin_density_real_branch
  use multiplicative_nlo_decay, only: multiplicative_nlo_workspace, &
       multiplicative_component_cache, &
       prepare_generated_multiplicative_workspace, &
       set_multiplicative_weight_count, reset_multiplicative_leaf_iterator, &
       next_multiplicative_leaf, capture_multiplicative_snapshot, &
       set_multiplicative_real_configuration, &
       multiplicative_leaf_has_snapshots, &
       complete_multiplicative_zero_branches, &
       invalidate_multiplicative_component_cache, &
       store_multiplicative_component_cache, &
       restore_multiplicative_component_cache, &
       contract_multiplicative_leaf
  use multiplicative_event_materialization, only: &
       multiplicative_leaf_particle_capacity, &
       materialize_multiplicative_leaf
  use fnlo_process_common, only: nfksprocess, soft_counterevent, &
                                 collinear_counterevent, &
                                 soft_collinear_counterevent, &
                                 real_event, &
                                 ybst_til_tolab, force_polecheck, &
                                 orders_tag_plot
  implicit none
  private

  integer, allocatable, save :: generated_mapconfig(:, :)
  logical, save :: generated_data_initialized = .false.
  public :: run_mintfo_driver
  public :: init_driver_generated_data
  public :: sigint_impl
  public :: get_user_params_impl

  interface
    subroutine init_process_dimensions_bridge()
    end subroutine init_process_dimensions_bridge

    subroutine init_born_dimensions_bridge()
    end subroutine init_born_dimensions_bridge

    subroutine init_fks_metadata_bridge()
    end subroutine init_fks_metadata_bridge

    subroutine init_fks_singular_bridge()
    end subroutine init_fks_singular_bridge

    integer function drv_getpid()
    end function drv_getpid

    subroutine setrun()
    end subroutine setrun

    subroutine setpara(parameter_card)
      character(len=*), intent(in) :: parameter_card
    end subroutine setpara

    subroutine setcuts()
    end subroutine setcuts

    subroutine printout()
    end subroutine printout

    subroutine get_user_params(ncall, nitmax, irestart)
      integer, intent(out) :: ncall, nitmax, irestart
    end subroutine get_user_params

    double precision function sigint(xx, vegas_wgt, ifl, f)
      use mint_module, only: ndimmax, nintegrals
      double precision, intent(in) :: xx(ndimmax), vegas_wgt
      integer, intent(in) :: ifl
      double precision, intent(out) :: f(nintegrals)
    end function sigint

    subroutine fks_inc_chooser()
    end subroutine fks_inc_chooser

    subroutine leshouche_inc_chooser()
    end subroutine leshouche_inc_chooser

  end interface

contains

  subroutine init_driver_generated_data(mapconfig_in)
    implicit none
    integer, intent(in) :: mapconfig_in(0:, 0:)

    call validate_process_and_born_dimensions()
    if (ubound(mapconfig_in, 1) /= lmaxconfigs .or. &
        ubound(mapconfig_in, 2) /= fks_configs) then
      call fail_driver('generated configuration map has the wrong shape')
    end if

    if (generated_data_initialized) then
      if (any(generated_mapconfig /= mapconfig_in)) then
        call fail_driver('generated driver data changed after initialization')
      end if
      return
    end if

    allocate (generated_mapconfig(0:lmaxconfigs, 0:fks_configs))
    generated_mapconfig = mapconfig_in
    generated_data_initialized = .true.
  end subroutine init_driver_generated_data

  subroutine run_mintfo_driver(nndim, flat_grid, momcmp_count, &
                               xratmax, ntot, nsun, nsps, nups, neps, n100, nddp, nqdp, &
                               nini, n10, n1)
    implicit none
    integer, intent(inout) :: nndim, momcmp_count
    logical, intent(inout) :: flat_grid
    double precision, intent(inout) :: xratmax
    integer, intent(inout) :: ntot, nsun, nsps, nups, neps, n100
    integer, intent(inout) :: nddp, nqdp, nini, n10, n1(0:9)
    integer :: i, j, kchan
    integer :: restart_mode
    real :: time_before, time_after, time_other, time_total
    logical :: unequal_factorization_scales

    call init_process_dimensions_bridge()
    call init_born_dimensions_bridge()
    call init_fks_metadata_bridge()
    call validate_process_and_born_dimensions()
    call validate_fks_metadata()
    force_polecheck = .false.

    write (*, *) drv_getpid()
    call reset_timing_state()
    call cpu_time(time_before)

    call FKSParamReader()
    min_virt_fraction_mint = min_virt_fraction
    do kchan = 1, maxchannels
      do i = 0, n_ave_virt
        average_virtual(i, kchan) = 0d0
      end do
      virtual_fraction(kchan) = max(virt_fraction, min_virt_fraction)
    end do
    n_ord_virt = nlo_virtual_grid_count(amp_split_size)

    ntot = 0
    nsun = 0
    nsps = 0
    nups = 0
    neps = 0
    n100 = 0
    nddp = 0
    nqdp = 0
    nini = 0
    n10 = 0
    n1 = 0

    call setrun()
    call setpara('param_card.dat')
    call configure_fnlo_scale_variations()
    call init_fks_singular_bridge()
    call setcuts()
    call printout()
    call write_run_summary()
    call fill_configurations_common()
    call check_amp_split()

    write (*, *) 'getting user params'
    call get_user_params(ncalls0, itmax, restart_mode)
    imode = restart_mode
    flat_grid = imode == 0

    nndim = real_phase_space_dimension()
    if (abs(lpp(1)) >= 1) nndim = nndim + 1
    if (abs(lpp(2)) >= 1) nndim = nndim + 1
    ! ``nndim`` remains the canonical phase-space layout consumed by genps:
    ! common Born variables followed by one (xi,y,phi) triple.  MINT owns an
    ! independent copy of that triple for every factorized NLO block.
    ndim = factorized_integration_dimension(nndim)
    if (has_nlo_contribution_bundle()) then
      call configure_decay_folding( &
           nndim, decay_variable_start(), decay_random_dimension())
    end if

    if (fixed_fac_scale) then
      unequal_factorization_scales = abs( &
                                     muf1_over_ref*muf1_ref_fixed - &
                                     muf2_over_ref*muf2_ref_fixed) > 0d0
    else
      unequal_factorization_scales = &
        abs(muf1_over_ref - muf2_over_ref) > 0d0
    end if
    if (unequal_factorization_scales) then
      write (*, *) 'NLO computations require muF1=muF2'
      stop
    end if

    write (*, *) 'about to integrate ', ndim, ncalls0, itmax
    momcmp_count = 0
    xratmax = 0d0
    if (imode == -1 .or. imode == 0) then
      if (imode == 0) then
        doreweight = .false.
        do_rwgt_scale = .false.
        do_rwgt_decay_scale = .false.
        do_rwgt_pdf = .false.
      else
        doreweight = do_rwgt_scale .or. do_rwgt_decay_scale .or. &
             do_rwgt_pdf
      end if
      write (*, *) 'imode is ', imode
      call mint(sigint)
      call topout_impl()
      call write_bundle_contribution_results()
      call deallocate_weight_lines()
    else
      write (*, *) 'Unknown imode', imode
      stop
    end if

    if (ntot /= 0) then
      write (*, *) 'Satistics from MadLoop:'
      write (*, *)
      write (*, *) &
        '  Total points tried:                              ', ntot
      write (*, *) &
        '  Stability unknown:                               ', nsun
      write (*, *) &
        '  Stable PS point:                                 ', nsps
      write (*, *) &
        '  Unstable PS point (and rescued):                 ', nups
      write (*, *) &
        '  Exceptional PS point (unstable and not rescued): ', neps
      write (*, *) &
        '  Double precision used:                           ', nddp
      write (*, *) &
        '  Quadruple precision used:                        ', nqdp
      write (*, *) &
        '  Initialization phase-space points:               ', nini
      write (*, *) &
        '  Unknown return code (100):                       ', n100
      write (*, *) &
        '  Unknown return code (10):                        ', n10
      write (*, *) &
        '  Unit return code distribution (1):               '
      do j = 0, 9
        if (n1(j) /= 0) write (*, *) '#Unit ', j, ' = ', n1(j)
      end do
    end if

    call cpu_time(time_after)
    time_total = time_after - time_before
    time_other = time_total - (tBorn + tGenPS + tReal + tCount + tIS + &
                               tf_nb + tf_all + t_as + tr_s + tr_pdf + t_plot + &
                               t_cuts + t_isum + t_coupl)
    write (*, *) 'Time spent in Born : ', tBorn
    write (*, *) 'Time spent in PS_Generation : ', tGenPS
    write (*, *) 'Time spent in Reals_evaluation: ', tReal
    write (*, *) 'Time spent in Counter_terms : ', tCount
    write (*, *) 'Time spent in Integrated_CT : ', tIS - tOLP
    write (*, *) 'Time spent in Virtuals : ', tOLP
    write (*, *) 'Time spent in Nbody_prefactor : ', tf_nb
    write (*, *) 'Time spent in N1body_prefactor : ', tf_all
    write (*, *) 'Time spent in Adding_alphas_pdf : ', t_as
    write (*, *) 'Time spent in Reweight_scale : ', tr_s
    write (*, *) 'Time spent in Reweight_pdf : ', tr_pdf
    write (*, *) 'Time spent in Filling_plots : ', t_plot
    write (*, *) 'Time spent in Applying_cuts : ', t_cuts
    write (*, *) 'Time spent in Sum_ident_contr : ', t_isum
    write (*, *) 'Time spent in AlphaS_dependencies : ', t_coupl
    write (*, *) 'Time spent in Other_tasks : ', time_other
    write (*, *) 'Time spent in Total : ', time_total

    open (unit=12, file='res.dat', status='unknown')
    do kchan = 0, nchans
      write (12, *) ans_result(1, kchan), unc_result(1, kchan), &
        ans_result(2, kchan), unc_result(2, kchan), itmax, ncalls0, &
        time_total
    end do
    close (12)

    if (momcmp_count /= 0) then
      write (*, *) '     '
      write (*, *) 'WARNING: genps_fks code 555555'
      write (*, *) momcmp_count, xratmax
    end if
  end subroutine run_mintfo_driver


  subroutine write_bundle_contribution_results()
    implicit none
    integer :: component, component_count, integral, unit_number, ios
    double precision :: component_sum
    character(len=96) :: label

    if (.not. has_nlo_contribution_bundle()) return
    ! The expanded component bins are not a decomposition of the
    ! multiplicative product.  Omitting this file is less misleading than
    ! reporting zero components and a spurious non-zero closure failure.
    if (multiplicative_nlo_enabled()) return
    component_count = bundle_component_count()
    open(newunit=unit_number, file='contribution_results.dat', &
         status='replace', action='write', iostat=ios)
    if (ios /= 0) call fail_driver('cannot open contribution_results.dat')
    write(unit_number, '(a)') 'FORMAT 1'
    write(unit_number, '(a,1x,i0)') 'COUNT', component_count
    component_sum = 0d0
    do component = 1, component_count
      integral = first_bundle_component_integral + component - 1
      call bundle_component_label(component, label)
      write(unit_number, '(a,1x,i0,1x,a,1x,es24.16,1x,es24.16)') &
           'COMPONENT', component, trim(label), ans_result(integral, 0), &
           unc_result(integral, 0)
      component_sum = component_sum + ans_result(integral, 0)
    end do
    write(unit_number, '(a,1x,es24.16)') 'SUM_COMPONENTS', component_sum
    write(unit_number, '(a,1x,es24.16,1x,es24.16)') &
         'TOTAL', ans_result(2, 0), unc_result(2, 0)
    write(unit_number, '(a,1x,es24.16)') &
         'CLOSURE', component_sum - ans_result(2, 0)
    write(unit_number, '(a)') 'END'
    close(unit_number)
  end subroutine write_bundle_contribution_results

  double precision function sigint_impl(xx, vegas_wgt, ifl, f, &
                                        ini_fin_fks, nndim, nbody, event_momenta, p_born, virtual_over_born, &
                                        calculated_born, abrv, wgt_me_born, wgt_me_real)
    implicit none
    double precision, intent(in) :: xx(ndimmax), vegas_wgt
    integer, intent(in) :: ifl, ini_fin_fks(maxchannels), nndim
    double precision, intent(out) :: f(nintegrals)
    logical, intent(inout) :: nbody, calculated_born
    double precision, intent(inout) :: &
      event_momenta(0:3, nexternal, soft_counterevent:real_event)
    double precision, intent(inout) :: p_born(0:3, nexternal - 1)
    double precision, intent(inout) :: virtual_over_born
    character(len=4), intent(in) :: abrv
    double precision, intent(inout) :: wgt_me_born, wgt_me_real
    double precision :: jacobian, momentum(0:3, nexternal)
    double precision :: reweight, volume, sampled_weight
    double precision :: vegas_variables(99)
    integer :: nfks_born, nfks_picked
    integer :: amplitude_order, picked_integer, radiation_block
    integer :: conditional_integer, conditional_dimension
    integer :: nbody_contribution, nbody_contribution_max
    integer :: selected_contribution
    logical :: passcuts_nbody, conditional_decay_integer
    logical :: nplusone_evaluated, skip_nplusone
    integer, save :: folded_picked_integer = 0
    double precision, save :: folded_integer_volume = 0d0
    logical, save :: folded_integer_is_valid = .false.
    integer, save :: folded_selected_contribution = 0
    double precision, save :: folded_contribution_volume = 0d0

    ! Bjorken x and CORE factorization scales are production data.  Preserve
    ! their PDF working set across decay replicas; genuinely DECAYED scales
    ! miss the numerical cache key and are evaluated normally.
    if (.not. decay_fold_reuses_production()) &
         call reset_luminosity_cache()
    if (multiplicative_nlo_enabled()) then
      sigint_impl = sigint_multiplicative_impl( &
           xx, vegas_wgt, ifl, f, ini_fin_fks, nndim, nbody, &
           event_momenta, p_born, virtual_over_born, calculated_born, &
           abrv, wgt_me_born, wgt_me_real)
      return
    end if

    if (new_point .and. ifl /= 2) pass_cuts_check = .false.
    call print_fks_channel_map()
    if (ifl /= 0) then
      write (*, *) 'ERROR ifl not equal to zero in sigint', ifl
      stop 1
    end if
    sigint_impl = 0d0
    icontr = 0
    do amplitude_order = 0, n_ave_virt
      virt_wgt_mint(amplitude_order) = 0d0
      born_wgt_mint(amplitude_order) = 0d0
    end do
    virtual_over_born = 0d0
    wgt_me_born = 0d0
    wgt_me_real = 0d0
    nplusone_evaluated = .false.
    skip_nplusone = .false.
    conditional_decay_integer = .false.
    conditional_integer = 0
    conditional_dimension = 0

    if (decay_fold_reuses_production() .and. &
        folded_integer_is_valid) then
      selected_contribution = folded_selected_contribution
      if (contribution_is_nlo_decay(selected_contribution)) then
        conditional_dimension = additive_mc_dimension( &
             selected_contribution, ini_fin_fks(ichan))
        call get_mc_integer(conditional_dimension, &
             contribution_channel_count( &
                  selected_contribution, ini_fin_fks(ichan)), &
             conditional_integer, volume)
        volume = volume*folded_contribution_volume
        nfks_picked = contribution_channel_configuration( &
             selected_contribution, ini_fin_fks(ichan), &
             conditional_integer)
        conditional_decay_integer = .true.
        call record_decay_fold_conditional_sector()
      else
        picked_integer = folded_picked_integer
        volume = folded_integer_volume
        nfks_picked = fks_channel_configuration( &
             ini_fin_fks(ichan), picked_integer)
      end if
    else
      call get_mc_integer(max(ini_fin_fks(ichan), 1), &
                          fks_channel_count(ini_fin_fks(ichan)), &
                          picked_integer, volume)
      nfks_picked = fks_channel_configuration(ini_fin_fks(ichan), &
                                              picked_integer)
      selected_contribution = 1
      if (has_nlo_contribution_bundle()) &
           selected_contribution = contribution_for_fks(nfks_picked)
      if (decay_fold_group_active()) then
        folded_picked_integer = picked_integer
        folded_integer_volume = volume
        folded_selected_contribution = selected_contribution
        if (contribution_is_nlo_decay(selected_contribution)) then
          folded_contribution_volume = &
               additive_contribution_probability( &
                    selected_contribution, ini_fin_fks(ichan))
        else
          folded_contribution_volume = volume
        end if
        folded_integer_is_valid = .true.
      else
        folded_integer_is_valid = .false.
        folded_selected_contribution = 0
        folded_contribution_volume = 0d0
      end if
    end if

    if (abrv /= 'real') then
      nbody = .true.
      if (has_nlo_contribution_bundle()) then
        nbody_contribution_max = nlo_contribution_count()
        if (abrv == 'born') nbody_contribution_max = 1
      else
        nbody_contribution_max = 1
      end if

      do nbody_contribution = 1, nbody_contribution_max
        calculated_born = .false.
        if (has_nlo_contribution_bundle()) then
          nfks_born = contribution_representative_fks(nbody_contribution)
        else
          call get_born_fks_process(nfks_picked, nfks_born)
        end if
        radiation_block = 1
        if (has_nlo_contribution_bundle()) &
             radiation_block = nbody_contribution
        call update_vegas_x_impl(xx, vegas_variables, nndim, abrv, &
                                 radiation_block)
        call update_fks_dir_impl(nfks_born)
        if (ini_fin_fks(ichan) == 0) then
          jacobian = 1d0
        else
          jacobian = 0.5d0
        end if
        call generate_momenta(nndim, iconfig, jacobian, vegas_variables, &
                              momentum)
        if (p_born(0, 1) >= 0d0) then
          call compute_prefactors_nbody(vegas_wgt)
          passcuts_nbody = passcuts( &
                             event_momenta(0, 1, soft_counterevent), &
                             reweight, ybst_til_tolab(soft_counterevent))
          if (passcuts_nbody) then
            pass_cuts_check = .true.
            call set_alphas( &
              event_momenta(0:3, 1:nexternal, soft_counterevent))
            call include_multichannel_enhance(1)
            if (abrv(1:2) /= 'vi' .and. &
                (.not. has_nlo_contribution_bundle() .or. &
                 nbody_contribution == 1)) then
              call compute_born()
            end if
            if (abrv /= 'born') then
              if (abrv(1:2) /= 'vi' .and. &
                  has_nlo_contribution_bundle()) then
                call compute_decay_width_counterterm()
              end if
              if (has_nlo_contribution_bundle()) then
                call begin_bundle_virtual_tricks()
              end if
              call compute_nbody_noborn()
              if (has_nlo_contribution_bundle()) then
                call finish_bundle_virtual_tricks()
              end if
            end if
          end if
        else if (.not. has_nlo_contribution_bundle() .or. &
                 nbody_contribution == 1) then
          skip_nplusone = .true.
        end if

        ! Pair the one sampled resolved sector with its reduced sector while
        ! the latter's underlying Born decay tree is still resident.  The
        ! production contribution remains first, which is required when the
        ! process-specific decay metadata are initialized for the first time.
        if (has_nlo_contribution_bundle() .and. &
            nbody_contribution == selected_contribution .and. &
            abrv(1:4) /= 'born' .and. abrv(1:4) /= 'bovi' .and. &
            abrv(1:2) /= 'vi' .and. .not. skip_nplusone) then
          nbody = .false.
          call evaluate_additive_n1body_sector( &
               xx, vegas_wgt, nndim, abrv, nfks_picked, volume, &
               event_momenta, p_born, calculated_born, &
               wgt_me_born, wgt_me_real, .true.)
          nbody = .true.
          nplusone_evaluated = .true.
        end if
      end do
    end if

    if (abrv(1:4) /= 'born' .and. abrv(1:4) /= 'bovi' .and. &
        abrv(1:2) /= 'vi' .and. &
        .not. skip_nplusone .and. .not. nplusone_evaluated) then
      nbody = .false.
      call evaluate_additive_n1body_sector( &
           xx, vegas_wgt, nndim, abrv, nfks_picked, volume, &
           event_momenta, p_born, calculated_born, &
           wgt_me_born, wgt_me_real, .false.)
    end if

    call include_pdf_and_alphas()
    if (doreweight) then
      if (do_rwgt_scale .or. do_rwgt_decay_scale) call reweight_scale()
      if (do_rwgt_pdf) call reweight_pdf()
    end if

    call get_wgt_no_nbody(sampled_weight)
    if (conditional_decay_integer) then
      call fill_mc_integer(conditional_dimension, conditional_integer, &
                           abs(sampled_weight)*volume)
    else
      call fill_mc_integer(max(ini_fin_fks(ichan), 1), picked_integer, &
                           abs(sampled_weight)*volume)
    end if

    call fill_plots()
    call fill_mint_function(f)
  end function sigint_impl


  subroutine evaluate_additive_n1body_sector( &
       xx, vegas_wgt, nndim, abrv, ifks, volume, event_momenta, p_born, &
       calculated_born, wgt_me_born, wgt_me_real, reuse_born)
    implicit none
    double precision, intent(in) :: xx(ndimmax), vegas_wgt, volume
    integer, intent(in) :: nndim, ifks
    character(len=4), intent(in) :: abrv
    double precision, intent(inout) :: &
         event_momenta(0:3, nexternal, soft_counterevent:real_event)
    double precision, intent(inout) :: p_born(0:3, nexternal - 1)
    logical, intent(inout) :: calculated_born
    double precision, intent(inout) :: wgt_me_born, wgt_me_real
    logical, intent(in) :: reuse_born

    double precision :: jacobian, momentum(0:3, nexternal), reweight
    double precision :: vegas_variables(99)
    integer :: radiation_block
    logical :: passcuts_nbody, passcuts_n1body, passcuts_coll

    calculated_born = .false.
    wgt_me_born = 0d0
    wgt_me_real = 0d0
    jacobian = 1d0/volume
    radiation_block = 1
    if (has_nlo_contribution_bundle()) &
         radiation_block = contribution_for_fks(ifks)
    call update_vegas_x_impl( &
         xx, vegas_variables, nndim, abrv, radiation_block)
    call update_fks_dir_impl(ifks)
    if (reuse_born) then
      call generate_momenta_reusing_born( &
           nndim, iconfig, jacobian, vegas_variables, momentum)
    else
      call generate_momenta( &
           nndim, iconfig, jacobian, vegas_variables, momentum)
    end if
    if (p_born(0, 1) < 0d0) return

    call compute_prefactors_n1body(vegas_wgt)
    passcuts_nbody = passcuts( &
         event_momenta(0, 1, soft_counterevent), reweight, &
         ybst_til_tolab(soft_counterevent))
    passcuts_coll = passcuts_nbody
    if (.not. passcuts_coll) then
      passcuts_coll = passcuts( &
           event_momenta(0, 1, collinear_counterevent), reweight, &
           ybst_til_tolab(collinear_counterevent))
    end if
    passcuts_n1body = passcuts( &
         momentum, reweight, ybst_til_tolab(real_event))

    if (passcuts_nbody .and. abrv /= 'real') then
      pass_cuts_check = .true.
      call set_alphas( &
           event_momenta(0:3, 1:nexternal, soft_counterevent))
      call include_multichannel_enhance(3)
      call compute_soft_counter_term()
      call compute_soft_collinear_ct_impl()
    end if
    if (passcuts_coll .and. abrv /= 'real') then
      call set_alphas( &
           event_momenta(0:3, 1:nexternal, collinear_counterevent))
      call compute_collinear_counter_term()
    end if
    if (passcuts_n1body) then
      pass_cuts_check = .true.
      call set_alphas(momentum)
      call include_multichannel_enhance(2)
      call compute_real_emission()
    end if
  end subroutine evaluate_additive_n1body_sector


  double precision function sigint_multiplicative_impl( &
       xx, vegas_wgt, ifl, f, ini_fin_fks, nndim, nbody, &
       event_momenta, p_born, virtual_over_born, calculated_born, &
       abrv, wgt_me_born, wgt_me_real)
    implicit none
    double precision, intent(in) :: xx(ndimmax), vegas_wgt
    integer, intent(in) :: ifl, ini_fin_fks(maxchannels), nndim
    double precision, intent(out) :: f(nintegrals)
    logical, intent(inout) :: nbody, calculated_born
    double precision, intent(inout) :: &
         event_momenta(0:3, nexternal, soft_counterevent:real_event)
    double precision, intent(inout) :: p_born(0:3, nexternal - 1)
    double precision, intent(inout) :: virtual_over_born
    character(len=4), intent(in) :: abrv
    double precision, intent(inout) :: wgt_me_born, wgt_me_real

    type(multiplicative_nlo_workspace), save :: workspace
    type(multiplicative_component_cache), save :: folded_production_cache
    complex(kind=8), allocatable, save :: contracted_weights(:)
    double precision, allocatable, save :: leaf_weights(:), total_weights(:)
    double precision, allocatable, save :: leaf_momenta(:, :)
    double precision, allocatable, save :: volumes(:)
    integer, allocatable, save :: picked_integers(:), configurations(:)
    integer, allocatable, save :: statuses(:), pdgs(:)
    logical, allocatable, save :: particle_from_decay(:)
    double precision :: jacobian, momentum(0:3, nexternal)
    double precision :: vegas_variables(99), reweight
    double precision :: production_boost(0:1), imaginary_tolerance
    integer :: amplitude_order, category, channel_count, component
    integer :: component_position, contribution, contribution_count
    integer :: ibody, ifks, leaf_capacity, particle_count
    integer :: production_position, radiation_block, weight
    integer(kind=8) :: leaf_mask, leaves_seen
    logical :: available, leaf_is_materialized, pass_leaf
    logical :: production_contribution, reuse_folded_production
    logical, save :: folded_production_integer_is_valid = .false.
    real :: plot_time_before, plot_time_after

    if (new_point .and. ifl /= 2) pass_cuts_check = .false.
    call print_fks_channel_map()
    if (ifl /= 0) then
      write (*, *) 'ERROR ifl not equal to zero in multiplicative sigint', ifl
      stop 1
    end if
    if (.not. has_nlo_contribution_bundle()) then
      call fail_driver( &
           'multiplicative NLO requires a complete contribution bundle')
    end if
    if (abrv == 'born' .or. abrv == 'real' .or. abrv == 'bovi' .or. &
        abrv(1:2) == 'vi') then
      call fail_driver( &
           'multiplicative NLO requires the full Born/real/virtual run mode')
    end if

    sigint_multiplicative_impl = 0d0
    f = 0d0
    icontr = 0
    iwgt = 0
    do amplitude_order = 0, n_ave_virt
      virt_wgt_mint(amplitude_order) = 0d0
      born_wgt_mint(amplitude_order) = 0d0
    end do
    virtual_over_born = 0d0
    wgt_me_born = 0d0
    wgt_me_real = 0d0
    call clear_spin_density_weight_lines()
    call reset_spin_density_fks_matrices()
    call set_spin_density_fks_collection(.true.)
    call prepare_generated_multiplicative_workspace(workspace)

    ! A production component can be shared by all decay replicas only when
    ! its scales are functions of the fixed production core.  DECAYED scale
    ! choices deliberately retain the complete evaluation path.
    reuse_folded_production = decay_fold_reuses_production() .and. &
         folded_production_cache%valid .and. &
         .not. use_decayed_production_ren_scale_momenta()
    if (.not. decay_fold_reuses_production() .or. &
        use_decayed_production_ren_scale_momenta()) then
      call invalidate_multiplicative_component_cache( &
           folded_production_cache)
    end if
    if (decay_fold_group_active()) &
         call record_decay_fold_production_cache(reuse_folded_production)

    contribution_count = nlo_contribution_count()
    if (workspace%corrected_count /= contribution_count) then
      call fail_driver( &
           'generated B/R branches do not match the NLO contributions')
    end if
    if (allocated(picked_integers)) then
      if (size(picked_integers) /= contribution_count) then
        deallocate(picked_integers, configurations, volumes)
      end if
    end if
    if (.not. allocated(picked_integers)) then
      allocate(picked_integers(contribution_count))
      allocate(configurations(contribution_count))
      allocate(volumes(contribution_count))
    end if
    production_position = 0
    production_boost = 0d0
    if (decay_fold_group_active() .and. &
        .not. decay_fold_reuses_production()) then
      folded_production_integer_is_valid = .false.
    else if (.not. decay_fold_group_active()) then
      folded_production_integer_is_valid = .false.
    end if

    ! Every corrected physical block samples its own local FKS sector.  The
    ! initial/final outer channel partitions production only; all decay
    ! sectors use their complete local category.
    do contribution = 1, contribution_count
      production_contribution = &
           .not. contribution_is_nlo_decay(contribution)
      if (production_contribution) then
        category = ini_fin_fks(ichan)
      else
        category = 0
      end if
      channel_count = contribution_channel_count(contribution, category)
      if (production_contribution .and. &
          decay_fold_reuses_production()) then
        if (.not. folded_production_integer_is_valid) then
          call fail_driver( &
               'a decay fold lost its production FKS-channel choice')
        end if
      else
        call get_mc_integer( &
             multiplicative_mc_dimension(contribution, category, &
                                         production_contribution), &
             channel_count, picked_integers(contribution), &
             volumes(contribution))
        if (production_contribution .and. &
            decay_fold_group_active()) then
          folded_production_integer_is_valid = .true.
        end if
      end if
      configurations(contribution) = &
           contribution_channel_configuration( &
           contribution, category, picked_integers(contribution))
      call set_multiplicative_real_configuration( &
           workspace, contribution, configurations(contribution))
    end do

    ! Reduce each physical block and immediately evaluate its resolved
    ! sector.  Both maps share all but the final three random variables, so
    ! the second pass retains the complete underlying-Born decay tree and
    ! regenerates only the local FKS radiation.  No cuts are applied to these
    ! internal slots: their weighted densities are combined first.
    nbody = .true.
    do contribution = 1, contribution_count
      calculated_born = .false.
      wgt_me_born = 0d0
      wgt_me_real = 0d0
      ifks = contribution_representative_fks(contribution)
      radiation_block = contribution
      call update_vegas_x_impl( &
           xx, vegas_variables, nndim, abrv, radiation_block)
      call update_fks_dir_impl(ifks)
      production_contribution = &
           .not. contribution_is_nlo_decay(contribution)
      jacobian = 1d0
      if (production_contribution .and. &
          ini_fin_fks(ichan) /= 0) jacobian = 0.5d0
      call generate_momenta( &
           nndim, iconfig, jacobian, vegas_variables, momentum)
      if (p_born(0, 1) < 0d0) then
        call set_spin_density_fks_collection(.false.)
        return
      end if

      component_position = &
           workspace%contribution_positions(contribution)
      if (production_contribution) production_position = component_position
      ! The second solution of a massive real-emission map has no soft/Born
      ! counterevent.  It contributes only to R, so leave its B density and
      ! snapshot absent here; they are completed as an exact zero below.
      if (event_momenta(0, 1, soft_counterevent) > 0d0) then
        call set_alphas( &
             event_momenta(0:3, 1:nexternal, soft_counterevent))
        if (production_contribution .and. reuse_folded_production) then
          ! Rebuild only the cheap Born densities of the newly sampled decay
          ! blocks.  The expensive production B/R component is restored after
          ! the current decay lines have been PDF/alpha_s dressed.
          call record_multiplicative_spectator_born_densities( &
               component_position)
        else
          call compute_prefactors_nbody(vegas_wgt)
          ! The multiplicative contraction contains every physical block, so
          ! attach the global Born SDE partition exactly once, on production.
          if (production_contribution) call include_multichannel_enhance(1)
          call compute_born()
          call begin_bundle_virtual_tricks()
          call compute_nbody_noborn()
          call finish_bundle_virtual_tricks()
        end if

        if (production_contribution) then
          production_boost(spin_density_bornlike_branch) = &
               ybst_til_tolab(soft_counterevent)
          ! Production materializes one coherent LO embedding for all blocks;
          ! decay increments must not replace those Born snapshots.
          do component = 1, workspace%component_count
            call capture_multiplicative_snapshot( &
                 workspace, component, spin_density_bornlike_branch, &
                 soft_counterevent)
          end do
        end if
      end if

      ! Evaluate this block's resolved sector and local counterterms using
      ! the Born state retained above.  Counterterms feed B; only the
      ! resolved event feeds R.
      nbody = .false.
      ifks = configurations(contribution)
      calculated_born = .false.
      wgt_me_born = 0d0
      wgt_me_real = 0d0
      jacobian = 1d0/volumes(contribution)
      call update_vegas_x_impl( &
           xx, vegas_variables, nndim, abrv, contribution)
      call update_fks_dir_impl(ifks)
      call generate_momenta_reusing_born( &
           nndim, iconfig, jacobian, vegas_variables, momentum)
      if (p_born(0, 1) >= 0d0) then
        component_position = &
             workspace%contribution_positions(contribution)
        call capture_multiplicative_snapshot( &
             workspace, component_position, spin_density_real_branch, &
             real_event)
        if (production_contribution) then
          production_boost(spin_density_real_branch) = &
               ybst_til_tolab(real_event)
        end if

        if (.not. (production_contribution .and. &
                   reuse_folded_production)) then
          call compute_prefactors_n1body(vegas_wgt)
          if (event_momenta(0, 1, soft_counterevent) > 0d0) then
            call set_alphas( &
                 event_momenta(0:3, 1:nexternal, soft_counterevent))
            if (production_contribution) call include_multichannel_enhance(3)
            call compute_soft_counter_term()
            if (event_momenta(0, 1, soft_collinear_counterevent) > 0d0) then
              call compute_soft_collinear_ct_impl()
            end if
          end if
          if (event_momenta(0, 1, collinear_counterevent) > 0d0) then
            call set_alphas( &
                 event_momenta(0:3, 1:nexternal, collinear_counterevent))
            call compute_collinear_counter_term()
          end if
          call set_alphas(momentum)
          if (production_contribution) call include_multichannel_enhance(2)
          call compute_real_emission()
        end if
      end if
      nbody = .true.
    end do
    nbody = .false.
    if (production_position == 0) then
      call fail_driver('multiplicative branches contain no production block')
    end if

    leaf_capacity = multiplicative_leaf_particle_capacity()
    if (allocated(statuses)) then
      if (size(statuses) /= leaf_capacity) then
        deallocate(leaf_momenta, statuses, pdgs, particle_from_decay)
      end if
    end if
    if (.not. allocated(statuses)) then
      allocate(leaf_momenta(0:3, leaf_capacity))
      allocate(statuses(leaf_capacity), pdgs(leaf_capacity))
      allocate(particle_from_decay(leaf_capacity))
    end if

    call include_pdf_and_alphas()
    if (doreweight) then
      if (do_rwgt_scale .or. do_rwgt_decay_scale) call reweight_scale()
      if (do_rwgt_pdf) call reweight_pdf()
    end if
    if (icontr < 1) then
      ! A sampled factorized map can reject every local event.  This is an
      ! ordinary zero integrand point, not an incomplete branch workspace.
      ! Still classify its cut acceptance from the materialized leaves so
      ! MINT can distinguish a zero channel from inefficient phase space.
      call reset_multiplicative_leaf_iterator(workspace)
      do
        call next_multiplicative_leaf(workspace, leaf_mask, available)
        if (.not. available) exit
        call test_multiplicative_leaf_cuts( &
             workspace, production_position, production_boost, &
             leaf_momenta, statuses, pdgs, particle_from_decay, &
             particle_count, reweight, leaf_is_materialized, pass_leaf)
        if (.not. leaf_is_materialized) cycle
        if (pass_leaf) then
          pass_cuts_check = .true.
          exit
        end if
      end do
      do contribution = 1, contribution_count
        production_contribution = &
             .not. contribution_is_nlo_decay(contribution)
        if (production_contribution) then
          category = ini_fin_fks(ichan)
        else
          category = 0
        end if
        call fill_mc_integer( &
             multiplicative_mc_dimension( &
             contribution, category, production_contribution), &
             picked_integers(contribution), 0d0)
      end do
      call set_spin_density_fks_collection(.false.)
      return
    end if
    if (iwgt < 1) then
      call fail_driver('multiplicative NLO produced no evaluated weights')
    end if
    call set_multiplicative_weight_count(workspace, iwgt)
    call aggregate_spin_density_weight_lines(workspace)
    call complete_multiplicative_zero_branches(workspace)
    if (reuse_folded_production) then
      call restore_multiplicative_component_cache( &
           workspace, production_position, vegas_wgt, &
           folded_production_cache)
    else if (decay_fold_group_active() .and. &
             .not. use_decayed_production_ren_scale_momenta()) then
      call store_multiplicative_component_cache( &
           workspace, production_position, vegas_wgt, &
           folded_production_cache)
    end if

    if (allocated(contracted_weights)) then
      if (size(contracted_weights) /= iwgt) then
        deallocate(contracted_weights, leaf_weights, total_weights)
      end if
    end if
    if (.not. allocated(contracted_weights)) then
      allocate(contracted_weights(iwgt), leaf_weights(iwgt))
      allocate(total_weights(iwgt))
    end if
    total_weights = 0d0
    leaves_seen = 0_8
    call reset_multiplicative_leaf_iterator(workspace)
    do
      call next_multiplicative_leaf(workspace, leaf_mask, available)
      if (.not. available) exit
      leaves_seen = leaves_seen + 1_8
      call test_multiplicative_leaf_cuts( &
           workspace, production_position, production_boost, &
           leaf_momenta, statuses, pdgs, particle_from_decay, &
           particle_count, reweight, leaf_is_materialized, pass_leaf)
      if (.not. leaf_is_materialized) cycle
      if (.not. pass_leaf) cycle
      pass_cuts_check = .true.

      ! Cuts are independent of the matrix-element value.  Applying them
      ! before the spin contraction avoids both a second materialization and
      ! contractions for rejected leaves while preserving IR-safe analysis.
      call contract_multiplicative_leaf(workspace, contracted_weights)
      if (all(abs(contracted_weights) == 0d0)) cycle
      do weight = 1, iwgt
        imaginary_tolerance = 1d-8*max( &
             1d0, abs(real(contracted_weights(weight), kind=8)))
        if (abs(aimag(contracted_weights(weight))) > &
            imaginary_tolerance) then
          call fail_driver( &
               'a multiplicative spin-density contraction is not real')
        end if
        leaf_weights(weight) = real(contracted_weights(weight), kind=8)
      end do
      leaf_weights = leaf_weights*reweight
      total_weights = total_weights + leaf_weights

      ibody = 2
      if (any(workspace%branch_by_component == &
              spin_density_real_branch)) ibody = 1
      ! A multiplicative product spans several perturbative orders and has no
      ! unique legacy split-order tag.  Inclusive histograms remain defined.
      orders_tag_plot = -1
      call cpu_time(plot_time_before)
      call outfun_multiplicative_impl( &
           leaf_momenta, particle_count, &
           production_boost( &
           workspace%branch_by_component(production_position)), &
           leaf_weights, pdgs, statuses, ibody)
      call cpu_time(plot_time_after)
      t_plot = t_plot + plot_time_after - plot_time_before
    end do
    if (leaves_seen /= workspace%leaf_count) then
      call fail_driver('the multiplicative B/R iterator missed a leaf')
    end if

    do contribution = 1, contribution_count
      production_contribution = &
           .not. contribution_is_nlo_decay(contribution)
      if (production_contribution) then
        category = ini_fin_fks(ichan)
      else
        category = 0
      end if
      call fill_mc_integer( &
           multiplicative_mc_dimension(contribution, category, &
                                       production_contribution), &
           picked_integers(contribution), &
           abs(total_weights(1))*volumes(contribution))
    end do

    call fill_multiplicative_mint_function( &
         f, total_weights(1), virtual_over_born)
    call set_spin_density_fks_collection(.false.)
  end function sigint_multiplicative_impl


  subroutine test_multiplicative_leaf_cuts( &
       workspace, production_position, production_boost, leaf_momenta, &
       statuses, pdgs, particle_from_decay, particle_count, reweight, &
       materialized, passes)
    implicit none
    type(multiplicative_nlo_workspace), intent(in) :: workspace
    integer, intent(in) :: production_position
    double precision, intent(in) :: production_boost(0:1)
    double precision, intent(out) :: leaf_momenta(0:, :)
    integer, intent(out) :: statuses(:), pdgs(:), particle_count
    logical, intent(out) :: particle_from_decay(:)
    double precision, intent(out) :: reweight
    logical, intent(out) :: materialized, passes
    integer :: production_branch

    materialized = multiplicative_leaf_has_snapshots(workspace)
    passes = .false.
    reweight = 0d0
    particle_count = 0
    if (.not. materialized) return

    call materialize_multiplicative_leaf( &
         workspace, leaf_momenta, statuses, pdgs, &
         particle_from_decay, particle_count)
    production_branch = &
         workspace%branch_by_component(production_position)
    passes = passcuts_multiplicative( &
         leaf_momenta, particle_count, statuses, pdgs, &
         particle_from_decay, reweight, &
         production_boost(production_branch))
  end subroutine test_multiplicative_leaf_cuts


  integer function multiplicative_mc_dimension( &
       contribution, category, production_contribution)
    implicit none
    integer, intent(in) :: contribution, category
    logical, intent(in) :: production_contribution

    if (production_contribution) then
      multiplicative_mc_dimension = max(category, 1)
    else
      ! Dimensions one and two retain the production initial/final grids.
      multiplicative_mc_dimension = contribution + 2
    end if
  end function multiplicative_mc_dimension


  integer function additive_mc_dimension(contribution, category)
    implicit none
    integer, intent(in) :: contribution, category

    ! Keep separate conditional grids for the two outer initial/final channel
    ! partitions: their local configuration counts need not be identical.
    ! Dimensions one and two remain the legacy global FKS grids.
    additive_mc_dimension = 2 + 2*(contribution - 1) + &
         min(2, max(category, 1))
  end function additive_mc_dimension


  double precision function additive_contribution_probability( &
       contribution, category)
    integer, intent(in) :: contribution, category
    integer :: configuration, dimension, position

    dimension = max(category, 1)
    additive_contribution_probability = 0d0
    do position = 1, fks_channel_count(category)
      configuration = fks_channel_configuration(category, position)
      if (contribution_for_fks(configuration) /= contribution) cycle
      additive_contribution_probability = &
           additive_contribution_probability + &
           mc_integer_interval_volume(dimension, position)
    end do
    if (additive_contribution_probability <= 0d0 .or. &
        additive_contribution_probability > 1d0 + 64d0*epsilon(1d0)) then
      call fail_driver( &
           'an additive folded contribution has invalid probability')
    end if
  end function additive_contribution_probability


  subroutine fill_multiplicative_mint_function( &
       f, central_weight, virtual_ratio)
    implicit none
    double precision, intent(out) :: f(nintegrals)
    double precision, intent(in) :: central_weight, virtual_ratio
    integer :: amplitude_order, born_index, component
    integer :: component_integral, contribution, virtual_index

    f = 0d0
    f(1) = abs(central_weight)
    f(2) = central_weight
    f(4) = virtual_ratio
    ! A multiplicative product has no additive production/decay component
    ! decomposition.  The corresponding MINT slots are nevertheless the
    ! grid-learning inputs for each block's private radiation variables.
    ! Feed the complete integrand to every corrected-block proxy so those
    ! dimensions adapt to the function being integrated instead of seeing an
    ! identically zero grid weight.
    do contribution = 1, nlo_contribution_count()
      component = bundle_nlo_component(contribution)
      component_integral = first_bundle_component_integral + component - 1
      if (component_integral > nintegrals) then
        call fail_driver('a multiplicative grid proxy is out of range')
      end if
      f(component_integral) = central_weight
    end do
    do amplitude_order = 0, n_ord_virt
      if (amplitude_order == 0) then
        f(3) = 0d0
        f(5) = 0d0
        f(6) = 0d0
        do virtual_index = 1, n_ord_virt
          f(3) = f(3) + virt_wgt_mint(virtual_index)
          f(6) = f(6) + born_wgt_mint(virtual_index)
        end do
        f(5) = abs(f(3))
      else
        virtual_index = 2*amplitude_order + 5
        born_index = 2*amplitude_order + 6
        f(virtual_index) = virt_wgt_mint(amplitude_order)
        f(born_index) = born_wgt_mint(amplitude_order)
      end if
    end do
  end subroutine fill_multiplicative_mint_function

  subroutine update_fks_dir_impl(nfks)
    implicit none
    integer, intent(in) :: nfks

    nfksprocess = nfks
    call fks_inc_chooser()
    call leshouche_inc_chooser()
    call setcuts()
    call setfksfactor()
  end subroutine update_fks_dir_impl


  subroutine update_vegas_x_impl(xx, x, nndim, abrv, radiation_block)
    implicit none
    double precision, intent(in) :: xx(ndimmax)
    double precision, intent(out) :: x(99)
    integer, intent(in) :: nndim, radiation_block
    character(len=4), intent(in) :: abrv
    integer :: i, shared_dimensions, source_start

    x = 0d0
    if (nndim < 3 .or. nndim > 99) then
      call fail_driver('invalid canonical phase-space dimension')
    end if
    shared_dimensions = nndim - 3
    do i = 1, shared_dimensions
      x(i) = xx(i)
    end do
    if (abrv == 'born' .or. abrv(1:2) == 'vi') then
      do i = nndim - 2, nndim
        x(i) = random_unit_interval(iconfig)
      end do
    else
      source_start = factorized_radiation_start(nndim, radiation_block)
      if (source_start + 2 > ndimmax) then
        call fail_driver('factorized radiative variables exceed ndimmax')
      end if
      do i = 0, 2
        x(nndim - 2 + i) = xx(source_start + i)
      end do
    end if
  end subroutine update_vegas_x_impl

  subroutine get_user_params_impl(ncall, nitmax, restart_mode, &
                                  ini_fin_fks, isum_hel, multi_channel, abrv, nbody, &
                                  mc_hel, random_offset_split)
    implicit none
    integer, intent(out) :: ncall, nitmax, restart_mode
    integer, intent(inout) :: ini_fin_fks(maxchannels), isum_hel
    logical, intent(inout) :: multi_channel
    integer, intent(inout) :: mc_hel
    character(len=4), intent(inout) :: abrv
    logical, intent(inout) :: nbody
    integer, intent(inout) :: random_offset_split
    integer :: i, kchan, parsed_integer
    double precision :: configurations(maxchannels)
    character(len=5) :: abrv_input
    character(len=140) :: buffer
    logical :: done

    if (.not. generated_data_initialized) then
      call fail_driver('generated driver data are not initialized')
    end if
    call validate_process_and_born_dimensions()
    call validate_fks_metadata()

    open (unit=83, file='input_app.txt', status='old')
    done = .false.
    nchans = 0
    do while (.not. done)
      read (83, '(a)', err=222, end=222) buffer
      if (buffer(1:7) == 'NPOINTS') then
        buffer = buffer(10:100)
        read (buffer, *) ncall
        write (*, *) 'Number of phase-space points per iteration:', ncall
      else if (buffer(1:11) == 'NITERATIONS') then
        read (buffer(14:), *) nitmax
        write (*, *) 'Maximum number of iterations is:', nitmax
      else if (buffer(1:8) == 'ACCURACY') then
        read (buffer(11:), *) accuracy
        write (*, *) 'Desired accuracy is:', accuracy
      else if (buffer(1:10) == 'ADAPT_GRID') then
        read (buffer(13:), *) parsed_integer
        write (*, *) 'Using adaptive grids:', parsed_integer
      else if (buffer(1:12) == 'MULTICHANNEL') then
        read (buffer(15:), *) parsed_integer
        if (parsed_integer == 1) then
          multi_channel = .true.
          write (*, *) 'Using Multi-channel integration'
        else
          multi_channel = .false.
          write (*, *) 'Not using Multi-channel integration'
        end if
      else if (buffer(1:12) == 'SUM_HELICITY') then
        read (buffer(15:), *) parsed_integer
        if (nincoming == 1) then
          write (*, *) 'Sum over helicities in the virtuals'// &
            ' for decay process'
          mc_hel = 0
        else if (parsed_integer == 0) then
          mc_hel = 0
          write (*, *) 'Explicitly summing over helicities'// &
            ' for the virtuals'
        else
          mc_hel = 1
          write (*, *) 'Do MC over helicities for the virtuals'
        end if
        isum_hel = 0
      else if (buffer(1:6) == 'NCHANS') then
        read (buffer(9:), *) nchans
        write (*, *) 'Number of channels to integrate together:', nchans
        if (nchans > maxchannels) then
          write (*, *) 'Too many integration channels to be '// &
            'integrated together. Increase maxchannels', nchans, &
            maxchannels
          stop 1
        end if
      else if (buffer(1:7) == 'CHANNEL') then
        if (nchans <= 0) then
          write (*, *) '"NCHANS" missing in input files'// &
            ' (still zero)', nchans
          stop
        end if
        read (buffer(10:), *) (configurations(kchan), kchan=1, nchans)
        do kchan = 1, nchans
          iconfigs(kchan) = int(configurations(kchan))
          parsed_integer = nint(configurations(kchan)*10d0) - &
                           iconfigs(kchan)*10
          if (parsed_integer == 0) then
            ini_fin_fks(kchan) = 0
          else if (parsed_integer == 1) then
            ini_fin_fks(kchan) = 1
          else if (parsed_integer == 2) then
            ini_fin_fks(kchan) = 2
          else
            write (*, *) 'ERROR: invalid configuration number', &
              configurations
            stop 1
          end if
          do i = 1, generated_mapconfig(0, 0)
            if (iconfigs(kchan) == generated_mapconfig(i, 0)) then
              iconfigs(kchan) = i
              exit
            end if
          end do
        end do
        write (*, *) 'Running Configuration Number(s): ', &
          (iconfigs(kchan), kchan=1, nchans)
        write (*, *) 'initial-or-final', &
          (ini_fin_fks(kchan), kchan=1, nchans)
      else if (buffer(1:5) == 'SPLIT') then
        read (buffer(8:), *) random_offset_split
        write (*, *) 'Splitting channel:', random_offset_split
      else if (buffer(1:8) == 'WGT_MULT') then
        read (buffer(11:), *) wgt_mult
        write (*, *) 'Weight multiplier:', wgt_mult
      else if (buffer(1:8) == 'RUN_MODE') then
        read (buffer(11:), *) abrv_input
        if (abrv_input(5:5) == '0') then
          nbody = .true.
        else
          nbody = .false.
        end if
        abrv = abrv_input(1:4)
        write (*, *) 'doing the ', abrv, ' of this channel'
        if (nbody) then
          write (*, *) 'integration Born/virtual with Sfunction=1'
        else
          write (*, *) 'Normal integration (Sfunction != 1)'
        end if
      else if (buffer(1:7) == 'RESTART') then
        read (buffer(10:), *) restart_mode
        if (restart_mode == 0) then
          write (*, *) 'RESTART: Fresh run'
        else if (restart_mode == -1) then
          write (*, *) 'RESTART: Use old grids, but refil plots'
        else if (restart_mode == 1) then
          write (*, *) 'RESTART: continue with existing run'
        else
          write (*, *) 'RESTART:', restart_mode
        end if
      end if
      cycle
222   done = .true.
    end do
    close (83)

    if (fks_configs == 1) then
      if (pdg_type_d(1, fks_i_d(1)) == -21 .and. abrv /= 'born') then
        if (amp_split_size == 1) then
          write (*, *) 'Process generated with [LOonly=QCD]. '// &
            'Setting abrv to "born".'
          abrv = 'born'
        else
          write (*, *) 'Process only with virtual corrections'// &
            'Setting abrv to "bovi".'
          abrv = 'bovi'
        end if
      end if
    end if
  end subroutine get_user_params_impl

  subroutine fail_driver(message)
    implicit none
    character(len=*), intent(in) :: message

    write (*, *) 'ERROR in driver_mintFO: ', trim(message)
    stop 1
  end subroutine fail_driver

end module driver_mintfo_module
