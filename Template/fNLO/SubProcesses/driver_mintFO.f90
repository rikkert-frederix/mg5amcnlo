module driver_mintfo_module
  use run_printout_module, only: write_run_summary
  use extra_weights, only: doreweight, dyn_scale, scalevarR, scalevarF, &
       lhaPDFid, nmemPDF
  use mint_module, only: maxchannels, n_ave_virt, average_virtual, &
                         virtual_fraction, min_virt_fraction_mint, n_ord_virt, ncalls0, &
                         itmax, imode, ndim, ndimmax, nintegrals, nchans, &
                         iconfig, ichan, &
                         iconfigs, accuracy, wgt_mult, new_point, pass_cuts_check, &
                         virt_wgt_mint, born_wgt_mint, mint, &
                         first_bundle_component_integral, &
                         first_multiplicative_linear_integral, &
                         first_multiplicative_group_integral, &
                         multiplicative_lo_integral, &
                         multiplicative_additive_integral
  use mint_module, only: ans_result => ans, unc_result => unc
  use FKSParams, only: min_virt_fraction, virt_fraction, &
       PrecisionVirtualAtRunTime, FKSParamReader
  use weight_lines, only: icontr, deallocate_weight_lines
  use process_dimensions, only: nexternal, nincoming, fks_configs, &
                                amp_split_size, lmaxconfigs, &
                                validate_process_and_born_dimensions, &
                                configure_event_capacity
  use fks_metadata, only: fks_i_d, pdg_type_d, validate_fks_metadata
  use fks_channel_map, only: fks_channel_count, &
                             fks_channel_configuration, &
                             print_fks_channel_map, get_born_fks_process
  use fks_random_module, only: random_unit_interval
  use run_state, only: lpp, fixed_fac_scale, muf1_over_ref, &
                       muf2_over_ref, muf1_ref_fixed, muf2_ref_fixed, &
                       do_rwgt_scale, do_rwgt_decay_scale, do_rwgt_pdf, &
                       q2fact
  use genps_fks, only: generate_momenta
  use decay_chain_metadata, only: real_phase_space_dimension, &
       decay_node_count, node_pdg
  use fnlo_scale_variations, only: configure_fnlo_scale_variations, &
       fnlo_scale_point_count, decode_fnlo_scale_point
  use nlo_contribution_bundle, only: has_nlo_contribution_bundle, &
       nlo_contribution_count, contribution_representative_fks, &
       contribution_for_fks, factorized_integration_dimension, &
       factorized_radiation_start, &
       multiplicative_event_capacity, &
       contribution_fks_channel_count, &
       contribution_fks_channel_configuration, &
       multiplicative_mc_integer_dimension, &
       nlo_virtual_grid_count, bundle_component_count, &
       bundle_component_label, bundle_species_is_nlo, &
       bundle_nlo_component
  use decay_chain_parameters, only: uses_multiplicative_nlo_combination, &
       decay_width_denominator_rescaling, decay_scale_species_count, &
       decay_scale_species, decay_lo_width, decay_nlo_width
  use multiplicative_phase_space, only: &
       multiplicative_phase_space_assembly, &
       initialize_multiplicative_phase_space_assembly, &
       capture_multiplicative_contribution, &
       restore_multiplicative_phase_space_assembly
  use multiplicative_density_terms, only: block_nlo_distribution, &
       multiplicative_density_tuple, density_cartesian_tuple_count, &
       decode_density_cartesian_tuple
  use multiplicative_block_distribution, only: &
       build_multiplicative_block_nlo_distribution, &
       build_multiplicative_lo_only_distribution
  use multiplicative_kinematics, only: realize_factorized_event_tuple
  use multiplicative_scale_state, only: &
       initialize_multiplicative_scale_references, &
       activate_multiplicative_block_reference, &
       build_multiplicative_scale_tables
  use spin_density_matrix_results, only: reset_spin_density_caches
  use multiplicative_production_channels, only: &
       multiplicative_production_channel_partition
  use multiplicative_runtime, only: multiplicative_event_evaluation, &
       evaluate_multiplicative_event_tuple
  use multiplicative_lambda_validation, only: &
       multiplicative_lambda_accumulator, &
       initialize_multiplicative_lambda_accumulator, &
       accumulate_multiplicative_lambda_atom, &
       formal_lambda_lo_weight, formal_lambda_additive_weight, &
       formal_lambda_block_linear_correction, &
       require_formal_lambda_closure, &
       require_formal_lambda_linear_closure
  use setscales_module, only: set_alphas
  use split_orders, only: check_amp_split
  use cuts_module, only: passcuts, passcuts_multiplicative
  use mc_integer_module, only: get_mc_integer, fill_mc_integer
  use timing_state, only: reset_timing_state, tBorn, tIS, tReal, &
                          tCount, tf_nb, tf_all, t_as, tr_s, tr_pdf, t_plot, &
                          t_cuts, t_isum, tOLP, tGenPS, t_coupl
  use fks_contributions_module, only: compute_prefactors_nbody, &
       compute_prefactors_n1body, include_multichannel_enhance, &
       compute_born, compute_decay_width_counterterm, &
       compute_nbody_noborn, compute_soft_counter_term, &
       compute_soft_collinear_ct_impl, compute_collinear_counter_term, &
       compute_real_emission
  use fks_weights_module, only: include_pdf_and_alphas, reweight_scale, &
       reweight_pdf, get_wgt_no_nbody, fill_plots, fill_mint_function, &
       begin_bundle_virtual_tricks, finish_bundle_virtual_tricks
  use fks_singular_module, only: fill_configurations_common, setfksfactor, &
       evaluate_born_matrix
  use madfks_plot_module, only: topout_impl, initplot_impl, &
       outfun_multiplicative_impl
  use fnlo_process_common, only: nfksprocess, soft_counterevent, &
                                 collinear_counterevent, &
                                 soft_collinear_counterevent, &
                                 real_event, &
                                 ybst_til_tolab, force_polecheck, &
                                 amp2, config_map, diagramsymmetryfactor
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

    double precision function dlum(bjorken_x)
      double precision, intent(in) :: bjorken_x(2)
    end function dlum

    double precision function dlum_configuration(configuration, bjorken_x)
      integer, intent(in) :: configuration
      double precision, intent(in) :: bjorken_x(2)
    end function dlum_configuration

    subroutine InitPDFm(set_index, member_index)
      integer, intent(in) :: set_index, member_index
    end subroutine InitPDFm

    integer function sdm_multiplicative_block_count()
    end function sdm_multiplicative_block_count

    integer function sdm_multiplicative_physical_block(position)
      integer, intent(in) :: position
    end function sdm_multiplicative_physical_block

    integer function sdm_multiplicative_block_pdg(block)
      integer, intent(in) :: block
    end function sdm_multiplicative_block_pdg

    integer function sdm_contribution_component_position(contribution)
      integer, intent(in) :: contribution
    end function sdm_contribution_component_position

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
    ! ``nfksprocess`` lives in a generated COMMON block and is not guaranteed
    ! to have a defined value when a fresh MINT executable starts.  Several
    ! initialization routines query the context-sensitive decay metadata
    ! before the first integrand call.  Pin that context to the production
    ! representative so fresh-grid and refine/restart jobs see the same decay
    ! chain graph.
    if (has_nlo_contribution_bundle()) then
      nfksprocess = contribution_representative_fks(1)
    end if
    call configure_event_capacity(multiplicative_event_capacity())
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
      if (uses_multiplicative_nlo_combination()) call initplot_impl()
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
    if (uses_multiplicative_nlo_combination()) then
      call write_multiplicative_validation_results()
      return
    end if
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


  subroutine write_multiplicative_validation_results()
    integer :: unit_number, ios, contribution, group, integral, position
    double precision :: block_sum

    open(newunit=unit_number, file='multiplicative_validation_results.dat', &
         status='replace', action='write', iostat=ios)
    if (ios /= 0) then
      call fail_driver('cannot open multiplicative_validation_results.dat')
    end if
    write(unit_number, '(a)') 'FORMAT 2'
    write(unit_number, '(a,1x,es24.16,1x,es24.16)') &
         'EXACT', ans_result(2, 0), unc_result(2, 0)
    write(unit_number, '(a,1x,es24.16,1x,es24.16)') &
         'LAMBDA_ONE', ans_result(2, 0), unc_result(2, 0)
    write(unit_number, '(a,1x,es24.16)') 'LAMBDA_ONE_CLOSURE', &
         0d0
    write(unit_number, '(a,1x,es24.16,1x,es24.16)') &
         'LO', ans_result(multiplicative_lo_integral, 0), &
         unc_result(multiplicative_lo_integral, 0)
    write(unit_number, '(a,1x,es24.16,1x,es24.16)') &
         'ADDITIVE', ans_result(multiplicative_additive_integral, 0), &
         unc_result(multiplicative_additive_integral, 0)
    ! The derivative is reported as a correlated central-value difference.
    ! Its separate uncertainty would require retaining the LO/additive
    ! covariance, so do not manufacture one from the marginal errors.
    write(unit_number, '(a,1x,es24.16)') 'LINEAR_CORRECTION', &
         ans_result(multiplicative_additive_integral, 0) - &
         ans_result(multiplicative_lo_integral, 0)
    write(unit_number, '(a,1x,i0)') &
         'BLOCK_LINEAR_COUNT', nlo_contribution_count()
    block_sum = 0d0
    do contribution = 1, nlo_contribution_count()
      integral = first_multiplicative_linear_integral + contribution - 1
      position = sdm_contribution_component_position(contribution)
      write(unit_number, &
           '(a,1x,i0,1x,i0,1x,es24.16,1x,es24.16)') &
           'BLOCK_LINEAR', contribution, &
           sdm_multiplicative_physical_block(position), &
           ans_result(integral, 0), unc_result(integral, 0)
      block_sum = block_sum + ans_result(integral, 0)
    end do
    write(unit_number, '(a,1x,es24.16)') &
         'BLOCK_LINEAR_SUM', block_sum
    write(unit_number, '(a,1x,es24.16)') &
         'BLOCK_LINEAR_CLOSURE', block_sum - &
         (ans_result(multiplicative_additive_integral, 0) - &
          ans_result(multiplicative_lo_integral, 0))
    write(unit_number, '(a,1x,i0)') &
         'BLOCK_LINEAR_GROUP_COUNT', 3*nlo_contribution_count()
    do contribution = 1, nlo_contribution_count()
      position = sdm_contribution_component_position(contribution)
      do group = 1, 3
        integral = first_multiplicative_group_integral + &
             3*(contribution - 1) + group - 1
        write(unit_number, &
             '(a,1x,i0,1x,i0,1x,i0,1x,es24.16,1x,es24.16)') &
             'BLOCK_LINEAR_GROUP', contribution, &
             sdm_multiplicative_physical_block(position), group, &
             ans_result(integral, 0), unc_result(integral, 0)
      end do
    end do
    write(unit_number, '(a)') 'END'
    close(unit_number)
  end subroutine write_multiplicative_validation_results

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
    double precision :: vegas_variables(99), mc_integer_weight
    integer :: nfks_born, nfks_picked, ifks, nfks_min, nfks_max
    integer :: amplitude_order, picked_integer, position, radiation_block
    integer :: nbody_contribution, nbody_contribution_max
    logical :: passcuts_nbody, passcuts_n1body, passcuts_coll
    logical :: skip_nplusone

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
    skip_nplusone = .false.

    if (uses_multiplicative_nlo_combination()) then
      call sigint_multiplicative_impl( &
           xx, vegas_wgt, f, ini_fin_fks, nndim, nbody, &
           event_momenta, p_born, abrv)
      return
    end if

    call get_mc_integer(max(ini_fin_fks(ichan), 1), &
                        fks_channel_count(ini_fin_fks(ichan)), picked_integer, &
                        volume)
    nfks_picked = fks_channel_configuration(ini_fin_fks(ichan), &
                                            picked_integer)

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
      end do
    end if

    if (abrv(1:4) /= 'born' .and. abrv(1:4) /= 'bovi' .and. &
        abrv(1:2) /= 'vi' .and. &
        .not. skip_nplusone) then
      nbody = .false.
      nfks_min = picked_integer
      nfks_max = picked_integer
      mc_integer_weight = 1d0/volume

      do position = nfks_min, nfks_max
        ifks = fks_channel_configuration(ini_fin_fks(ichan), position)
        calculated_born = .false.
        wgt_me_born = 0d0
        wgt_me_real = 0d0
        jacobian = mc_integer_weight
        radiation_block = 1
        if (has_nlo_contribution_bundle()) &
             radiation_block = contribution_for_fks(ifks)
        call update_vegas_x_impl(xx, vegas_variables, nndim, abrv, &
                                 radiation_block)
        call update_fks_dir_impl(ifks)
        call generate_momenta(nndim, iconfig, jacobian, &
                              vegas_variables, momentum)
        if (p_born(0, 1) < 0d0) cycle

        call compute_prefactors_n1body(vegas_wgt)
        passcuts_nbody = passcuts( &
                         event_momenta(0, 1, soft_counterevent), reweight, &
                         ybst_til_tolab(soft_counterevent))
        passcuts_coll = passcuts_nbody
        if (.not. passcuts_coll) then
          passcuts_coll = passcuts( &
                            event_momenta(0, 1, collinear_counterevent), &
                            reweight, &
                            ybst_til_tolab(collinear_counterevent))
        end if
        passcuts_n1body = passcuts( &
                             momentum, reweight, &
                             ybst_til_tolab(real_event))

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
      end do
    end if

    call include_pdf_and_alphas()
    if (doreweight) then
      if (do_rwgt_scale .or. do_rwgt_decay_scale) call reweight_scale()
      if (do_rwgt_pdf) call reweight_pdf()
    end if

    call get_wgt_no_nbody(sampled_weight)
    call fill_mc_integer(max(ini_fin_fks(ichan), 1), picked_integer, &
                         abs(sampled_weight)*volume)

    call fill_plots()
    call fill_mint_function(f)
  end function sigint_impl


  subroutine sigint_multiplicative_impl( &
       xx, vegas_wgt, f, ini_fin_fks, nndim, nbody, event_momenta, &
       p_born, abrv)
    double precision, intent(in) :: xx(ndimmax), vegas_wgt
    double precision, intent(out) :: f(nintegrals)
    integer, intent(in) :: ini_fin_fks(maxchannels), nndim
    logical, intent(inout) :: nbody
    double precision, intent(inout) :: &
         event_momenta(0:3,nexternal,soft_counterevent:real_event)
    double precision, intent(inout) :: p_born(0:3,nexternal-1)
    character(len=4), intent(in) :: abrv
    type(multiplicative_phase_space_assembly) :: assembly
    type(block_nlo_distribution), allocatable :: distributions(:)
    type(block_nlo_distribution) :: distribution
    type(multiplicative_density_tuple) :: tuple
    type(multiplicative_event_evaluation) :: evaluation
    type(multiplicative_event_evaluation) :: variation_evaluation
    type(multiplicative_lambda_accumulator) :: lambda_validation
    double precision :: vegas_variables(99), momentum(0:3,nexternal)
    double precision :: jacobian, volume
    double precision :: logarithmic_mu2_r(0:nexternal)
    double precision :: logarithmic_mu2_f(0:nexternal)
    double precision :: coupling_rescaling(0:nexternal,0:1)
    double precision :: production_mu2_r, production_mu2_f
    double precision :: luminosity, reweight, tuple_weight, total_weight
    double precision :: width_rescaling, production_channel_partition
    double precision :: channel_born_weight
    double precision :: nbody_channel_partition
    double precision :: raw_lo_weight, formal_lo_value, formal_block_value
    double precision :: formal_lo_rescaling
    double precision, allocatable :: validation_lo_widths(:)
    double precision, allocatable :: validation_nlo_widths(:)
    double precision, allocatable :: plotted_weight(:)
    double precision, allocatable :: radiation_grid_weights(:)
    double precision, allocatable :: radiation_grid_groups(:, :)
    double precision, allocatable :: linear_radiation_groups(:, :)
    integer, allocatable :: sampled_fks(:), sampled_integer(:)
    integer, allocatable :: sampled_dimension(:)
    integer, allocatable :: factor_indices(:)
    integer, allocatable :: contribution_positions(:)
    integer, allocatable :: validation_block_orders(:)
    integer, allocatable :: validation_width_blocks(:)
    double precision, allocatable :: sampled_volume(:)
    integer :: decay_block_factor_indices(0:nexternal)
    integer :: channel_event_slots(0:nexternal)
    logical, allocatable :: component_owned(:)
    integer :: contribution_count, contribution, category, channel_count
    integer :: picked, configuration, component_count, component_position
    integer :: position, block, tuple_index, tuple_count, sampled_count
    integer :: production_position, production_term_index
    integer :: radiation_term_index, radiation_event_slot
    integer :: radiation_grid_group
    integer :: luminosity_configuration
    integer :: weight_count, weight_index, dd, point, kr, kf
    integer :: pdf_set, pdf_member, species_count
    logical :: available, pass

    f = 0d0
    nbody = .false.
    if (trim(abrv) /= 'all') then
      call fail_driver( &
           'MULTIPLICATIVE mode requires one unsplit RUN_MODE=all job')
    end if
    contribution_count = nlo_contribution_count()
    ! The previous integrand call normally leaves the generated COMMON state
    ! on its last decay contribution.  Width metadata below describe the full
    ! chain and therefore require the production context, independently of
    ! which MINT channel this point will subsequently sample.
    call update_fks_dir_impl(contribution_representative_fks(1))
    component_count = sdm_multiplicative_block_count()
    if (contribution_count < 1 .or. component_count < contribution_count) then
      call fail_driver('the multiplicative block graph is incomplete')
    end if
    if (allocated(validation_block_orders)) then
      deallocate(validation_block_orders)
    end if
    if (allocated(validation_width_blocks)) then
      deallocate(validation_width_blocks)
    end if
    if (allocated(validation_lo_widths)) deallocate(validation_lo_widths)
    if (allocated(validation_nlo_widths)) deallocate(validation_nlo_widths)
    call initialize_multiplicative_lambda_accumulator( &
         lambda_validation, component_count)
    allocate(validation_block_orders(component_count))
    call initialize_validation_widths( &
         validation_width_blocks, validation_lo_widths, &
         validation_nlo_widths)
    if (allocated(sampled_fks)) deallocate(sampled_fks)
    if (allocated(sampled_integer)) deallocate(sampled_integer)
    if (allocated(sampled_dimension)) deallocate(sampled_dimension)
    if (allocated(sampled_volume)) deallocate(sampled_volume)
    if (allocated(distributions)) deallocate(distributions)
    if (allocated(component_owned)) deallocate(component_owned)
    if (allocated(factor_indices)) deallocate(factor_indices)
    if (allocated(plotted_weight)) deallocate(plotted_weight)
    if (allocated(radiation_grid_weights)) &
         deallocate(radiation_grid_weights)
    if (allocated(radiation_grid_groups)) &
         deallocate(radiation_grid_groups)
    if (allocated(linear_radiation_groups)) &
         deallocate(linear_radiation_groups)
    if (allocated(contribution_positions)) &
         deallocate(contribution_positions)
    allocate(sampled_fks(contribution_count))
    allocate(sampled_integer(contribution_count))
    allocate(sampled_dimension(contribution_count))
    allocate(sampled_volume(contribution_count))
    allocate(radiation_grid_weights(contribution_count))
    allocate(radiation_grid_groups(contribution_count, 3))
    allocate(linear_radiation_groups(contribution_count, 3))
    allocate(contribution_positions(contribution_count))
    allocate(distributions(component_count))
    allocate(component_owned(component_count))
    species_count = 0
    if (do_rwgt_decay_scale) species_count = decay_scale_species_count()
    allocate(factor_indices(species_count))
    weight_count = multiplicative_plot_weight_count()
    allocate(plotted_weight(weight_count))
    component_owned = .false.
    radiation_grid_weights = 0d0
    radiation_grid_groups = 0d0
    linear_radiation_groups = 0d0
    raw_lo_weight = 0d0
    contribution_positions = 0
    sampled_count = 0

    do contribution = 1, contribution_count
      category = 0
      if (contribution == 1) category = ini_fin_fks(ichan)
      channel_count = contribution_fks_channel_count( &
           contribution, category)
      if (channel_count < 1) return
      sampled_dimension(contribution) = &
           multiplicative_mc_integer_dimension(contribution, category)
      call get_mc_integer( &
           sampled_dimension(contribution), channel_count, picked, volume)
      sampled_integer(contribution) = picked
      sampled_volume(contribution) = volume
      sampled_count = contribution
      configuration = contribution_fks_channel_configuration( &
           contribution, category, picked)
      sampled_fks(contribution) = configuration

      call update_fks_dir_impl(configuration)
      call update_vegas_x_impl( &
           xx, vegas_variables, nndim, abrv, contribution)
      jacobian = 1d0/volume
      call generate_momenta( &
           nndim, iconfig, jacobian, vegas_variables, momentum)
      if (p_born(0,1) < 0d0) then
        call fill_multiplicative_discrete_grids()
        return
      end if
      if (contribution == 1) then
        call initialize_multiplicative_phase_space_assembly(assembly)
      end if
      call capture_multiplicative_contribution(assembly, contribution)
    end do
    call restore_multiplicative_phase_space_assembly(assembly)
    call initialize_multiplicative_scale_references()

    ! All terms in one production FKS family use one Born-level channel
    ! partition.  Materialize the all-Born tuple first so that every matrix
    ! element sees boosted block-local momenta, then let SBORN export the
    ! production provider's positive per-diagram amplitude weights.
    channel_event_slots = soft_counterevent
    call update_fks_dir_impl(sampled_fks(1))
    call realize_factorized_event_tuple(channel_event_slots, pass)
    if (.not. pass) then
      call fill_multiplicative_discrete_grids()
      return
    end if
    call activate_multiplicative_block_reference(0)
    call evaluate_born_matrix(soft_counterevent, channel_born_weight)
    production_channel_partition = &
         multiplicative_production_channel_partition( &
         ini_fin_fks, iconfigs, nchans, ichan, amp2, &
         config_map(:,0), diagramsymmetryfactor)
    ! The channel-weight contraction evaluates every spectator block while
    ! the production reference coupling is active.  Its AMP2 values have
    ! now been consumed; discard the accompanying density caches so the
    ! Cartesian contraction rebuilds each block at that block's reference.
    call reset_spin_density_caches()

    do contribution = 1, contribution_count
      call update_fks_dir_impl(sampled_fks(contribution))
      ! An n-body Born/virtual atom is common to every sampled FKS sector.
      ! Its block measure inherited the inverse discrete-sector probability,
      ! so cancel that probability here instead of summing several copies of
      ! the same atom.  Real and local counterevent atoms remain genuine
      ! sector sums and are deliberately left untouched.
      nbody_channel_partition = sampled_volume(contribution)
      ! Production is additionally present in both the initial- and
      ! final-state outer channel families.
      if (contribution == 1 .and. ini_fin_fks(ichan) /= 0) &
           nbody_channel_partition = 0.5d0*nbody_channel_partition
      call build_multiplicative_block_nlo_distribution( &
           contribution, distribution, available, &
           nbody_channel_partition)
      if (.not. available) then
        call fill_multiplicative_discrete_grids()
        return
      end if
      component_position = &
           sdm_contribution_component_position(contribution)
      if (component_position < 1 .or. &
          component_position > component_count .or. &
          component_owned(component_position)) then
        call fail_driver( &
             'the NLO-contribution/component map is not one-to-one')
      end if
      distributions(component_position) = distribution
      component_owned(component_position) = .true.
      contribution_positions(contribution) = component_position
    end do
    do position = 1, component_count
      if (component_owned(position)) cycle
      block = sdm_multiplicative_physical_block(position)
      call build_multiplicative_lo_only_distribution( &
           block, distributions(position))
      component_owned(position) = .true.
    end do

    production_position = contribution_positions(1)
    if (production_position < 1 .or. &
        production_position > component_count .or. &
        distributions(production_position)%block /= 0) then
      call fail_driver( &
           'the production contribution has no density component')
    end if

    tuple_count = density_cartesian_tuple_count(distributions)
    total_weight = 0d0
    do tuple_index = 1, tuple_count
      call decode_density_cartesian_tuple( &
           distributions, tuple_index, tuple)
      production_term_index = tuple%term_indices(production_position)
      luminosity_configuration = distributions(production_position)% &
           terms(production_term_index)%luminosity_configuration
      if (luminosity_configuration < 1 .or. &
          luminosity_configuration > fks_configs) then
        call fail_driver( &
             'a production density term has no luminosity configuration')
      end if
      call realize_factorized_event_tuple(tuple%event_slots, pass)
      if (.not. pass) cycle
      call build_multiplicative_scale_tables( &
           tuple%event_slots, logarithmic_mu2_r, logarithmic_mu2_f, &
           coupling_rescaling, production_mu2_r, production_mu2_f)
      call evaluate_multiplicative_event_tuple( &
           distributions, tuple_index, logarithmic_mu2_r, &
           logarithmic_mu2_f, coupling_rescaling, vegas_wgt, &
           PrecisionVirtualAtRunTime, &
           evaluation, .true.)
      if (.not. evaluation%available) cycle
      if (abs(aimag(evaluation%partonic_weight)) > &
          1d-8*max(1d0,abs(dble(evaluation%partonic_weight)))) then
        call fail_driver( &
             'a multiplicative density contraction is not real')
      end if

      q2fact = production_mu2_f
      luminosity = dlum_configuration( &
           luminosity_configuration, evaluation%bjorken_x)
      pass = passcuts_multiplicative( &
           evaluation%momenta, evaluation%pdgs, &
           evaluation%origin_blocks, evaluation%y_to_lab, reweight)
      if (.not. pass) cycle
      pass_cuts_check = .true.
      ! Generated DLUM providers already apply the GeV^-2-to-pb conversion
      ! for two incoming beams, exactly as in the additive weight path.
      tuple_weight = dble(evaluation%partonic_weight)*luminosity* &
           production_channel_partition
      tuple_weight = tuple_weight*reweight
      total_weight = total_weight + tuple_weight
      plotted_weight = 0d0
      plotted_weight(1) = tuple_weight

      weight_index = 1
      if (do_rwgt_scale .or. do_rwgt_decay_scale) then
        do dd = 1, dyn_scale(0)
          do point = 1, fnlo_scale_point_count(dd)
            call decode_fnlo_scale_point( &
                 dd, point, kr, kf, factor_indices)
            call map_decay_factor_indices( &
                 factor_indices, decay_block_factor_indices)
            call build_multiplicative_scale_tables( &
                 tuple%event_slots, logarithmic_mu2_r, logarithmic_mu2_f, &
                 coupling_rescaling, production_mu2_r, production_mu2_f, &
                 scalevarR(kr), scalevarF(kf), &
                 decay_block_factor_indices, dyn_scale(dd))
            call evaluate_multiplicative_event_tuple( &
                 distributions, tuple_index, logarithmic_mu2_r, &
                 logarithmic_mu2_f, coupling_rescaling, vegas_wgt, &
                 PrecisionVirtualAtRunTime, &
                 variation_evaluation, .true.)
            weight_index = weight_index + 1
            if (.not. variation_evaluation%available) cycle
            call require_real_multiplicative_weight( &
                 variation_evaluation%partonic_weight)
            q2fact = production_mu2_f
            luminosity = dlum_configuration( &
                 luminosity_configuration, &
                 variation_evaluation%bjorken_x)
            width_rescaling = &
                 decay_width_denominator_rescaling(factor_indices)
            plotted_weight(weight_index) = &
                 dble(variation_evaluation%partonic_weight)*luminosity* &
                 production_channel_partition*width_rescaling*reweight
          end do
        end do
      end if

      if (do_rwgt_pdf) then
        do pdf_set = 1, lhaPDFid(0)
          do pdf_member = 0, nmemPDF(pdf_set)
            call InitPDFm(pdf_set, pdf_member)
            call build_multiplicative_scale_tables( &
                 tuple%event_slots, logarithmic_mu2_r, logarithmic_mu2_f, &
                 coupling_rescaling, production_mu2_r, production_mu2_f)
            call evaluate_multiplicative_event_tuple( &
                 distributions, tuple_index, logarithmic_mu2_r, &
                 logarithmic_mu2_f, coupling_rescaling, vegas_wgt, &
                 PrecisionVirtualAtRunTime, &
                 variation_evaluation, .true.)
            weight_index = weight_index + 1
            if (.not. variation_evaluation%available) cycle
            call require_real_multiplicative_weight( &
                 variation_evaluation%partonic_weight)
            q2fact = production_mu2_f
            luminosity = dlum_configuration( &
                 luminosity_configuration, &
                 variation_evaluation%bjorken_x)
            plotted_weight(weight_index) = &
                 dble(variation_evaluation%partonic_weight)*luminosity* &
                 production_channel_partition*reweight
          end do
        end do
        call InitPDFm(1, 0)
      end if
      if (weight_index /= weight_count) then
        call fail_driver('the multiplicative plot-weight layout is invalid')
      end if
      call outfun_multiplicative_impl( &
           evaluation%momenta, evaluation%y_to_lab, plotted_weight, &
           evaluation%pdgs, evaluation%origin_blocks)
      call extract_validation_block_orders( &
           tuple, distributions, evaluation%nlo_order, &
           validation_block_orders)
      if (evaluation%nlo_order == 0) &
           raw_lo_weight = raw_lo_weight + tuple_weight
      do contribution = 1, contribution_count
        position = contribution_positions(contribution)
        if (position < 1 .or. position > component_count) then
          call fail_driver( &
               'an NLO contribution has no radiation-grid component')
        end if
        ! A block-local radiation grid learns the complete physical NLO
        ! distribution of that block with every spectator at LO.  The common
        ! all-LO tuple is therefore included in every block grid, together
        ! with the lambda-linear correction owned by that block.  R/S and
        ! C/SC generally live at different mapped momenta, so every atom is
        ! first evaluated and combined with the LO spectator blocks on its
        ! own momenta.  Only those final measured weights are paired here.
        ! This keeps the local FKS cancellation inside each positive grid
        ! proxy without allowing unrelated Born/virtual cancellations to
        ! make the proxy accidentally tiny.  Taking ABS atom by atom would
        ! instead train on the non-integrable 1/xi pieces that FKS
        ! subtraction is specifically meant to cancel.
        ! Products in which another block is NLO remain in TOTAL_WEIGHT but
        ! must not distort this one block's importance map with unrelated
        ! higher-order tails.
        if (evaluation%nlo_order == 0 .or. &
            (validation_block_orders(position) == 1 .and. &
             evaluation%nlo_order == 1)) then
          radiation_term_index = tuple%term_indices(position)
          radiation_event_slot = distributions(position)% &
               terms(radiation_term_index)%event_slot
          radiation_grid_group = 1
          if (radiation_event_slot == real_event .or. &
              (radiation_event_slot == soft_counterevent .and. &
               distributions(position)%terms(radiation_term_index)% &
               sign < 0)) then
            radiation_grid_group = 2
          else if (radiation_event_slot == collinear_counterevent .or. &
                   radiation_event_slot == &
                   soft_collinear_counterevent) then
            radiation_grid_group = 3
          end if
          radiation_grid_groups(contribution, radiation_grid_group) = &
               radiation_grid_groups(contribution, &
                                      radiation_grid_group) + tuple_weight
          if (evaluation%nlo_order == 1) then
            linear_radiation_groups(contribution, radiation_grid_group) = &
                 linear_radiation_groups(contribution, &
                                         radiation_grid_group) + tuple_weight
          end if
        end if
      end do
      call accumulate_multiplicative_lambda_atom( &
           lambda_validation, validation_block_orders, tuple_weight)
    end do

    call require_formal_lambda_closure( &
         lambda_validation, total_weight, validation_lo_widths, &
         validation_nlo_widths)
    call require_formal_lambda_linear_closure( &
         lambda_validation, validation_width_blocks, &
         validation_lo_widths, validation_nlo_widths)
    f(1) = abs(total_weight)
    f(2) = total_weight
    formal_lo_value = formal_lambda_lo_weight( &
         lambda_validation, validation_lo_widths, validation_nlo_widths)
    f(multiplicative_lo_integral) = formal_lo_value
    f(multiplicative_additive_integral) = formal_lambda_additive_weight( &
         lambda_validation, validation_lo_widths, validation_nlo_widths)
    do contribution = 1, contribution_count
      position = contribution_positions(contribution)
      f(first_multiplicative_linear_integral + contribution - 1) = &
           formal_lambda_block_linear_correction( &
           lambda_validation, position, validation_width_blocks, &
           validation_lo_widths, validation_nlo_widths)
      formal_block_value = &
           f(first_multiplicative_linear_integral + contribution - 1)
      formal_lo_rescaling = 0d0
      if (raw_lo_weight /= 0d0) &
           formal_lo_rescaling = formal_lo_value/raw_lo_weight
      linear_radiation_groups(contribution, 2:3) = &
           formal_lo_rescaling*linear_radiation_groups(contribution, 2:3)
      ! The Born-like group owns the signed n-body NLO atom and the formal
      ! derivative of this block's width denominator.  Defining it by
      ! closure keeps all three diagnostics finite at Born zeros.
      linear_radiation_groups(contribution, 1) = formal_block_value - &
           sum(linear_radiation_groups(contribution, 2:3))
      do radiation_grid_group = 1, 3
        f(first_multiplicative_group_integral + &
          3*(contribution - 1) + radiation_grid_group - 1) = &
             linear_radiation_groups(contribution, radiation_grid_group)
      end do
    end do
    do contribution = 1, contribution_count
      radiation_grid_weights(contribution) = &
           sum(abs(radiation_grid_groups(contribution, :)))
      component_position = bundle_nlo_component(contribution)
      f(first_bundle_component_integral + component_position - 1) = &
           radiation_grid_weights(contribution)
    end do
    call fill_multiplicative_discrete_grids(radiation_grid_weights)

  contains

    subroutine fill_multiplicative_discrete_grids(weights)
      double precision, intent(in), optional :: weights(:)
      double precision :: grid_weight
      integer :: selected

      ! A rejected block phase-space point can return before the later
      ! blocks have been sampled.  Only those grids whose dimension,
      ! interval and volume are already defined may be filled.  Each block's
      ! discrete selector learns from the same positive, locally subtracted
      ! training weight as that block's three continuous radiation grids,
      ! rather than from unrelated large corrections in another block.
      if (present(weights)) then
        if (size(weights) < sampled_count) then
          call fail_driver( &
               'a radiation-grid weight vector has the wrong size')
        end if
      end if
      do selected = 1, sampled_count
        grid_weight = 0d0
        if (present(weights)) grid_weight = weights(selected)
        call fill_mc_integer( &
             sampled_dimension(selected), sampled_integer(selected), &
             grid_weight*sampled_volume(selected))
      end do
    end subroutine fill_multiplicative_discrete_grids


    subroutine initialize_validation_widths( &
         width_blocks, lo_widths, nlo_widths)
      integer, allocatable, intent(out) :: width_blocks(:)
      double precision, allocatable, intent(out) :: lo_widths(:)
      double precision, allocatable, intent(out) :: nlo_widths(:)
      integer :: corrected_count, node, occurrence, owner, candidate
      integer :: pdg

      corrected_count = 0
      do node = 1, decay_node_count()
        if (bundle_species_is_nlo(node_pdg(node))) then
          corrected_count = corrected_count + 1
        end if
      end do
      allocate(width_blocks(corrected_count))
      allocate(lo_widths(corrected_count))
      allocate(nlo_widths(corrected_count))
      occurrence = 0
      do node = 1, decay_node_count()
        pdg = node_pdg(node)
        if (.not. bundle_species_is_nlo(pdg)) cycle
        owner = 0
        do candidate = 1, component_count
          if (sdm_multiplicative_physical_block(candidate) /= node) cycle
          owner = candidate
          exit
        end do
        if (owner == 0) then
          call fail_driver( &
               'a corrected width has no multiplicative decay block')
        end if
        occurrence = occurrence + 1
        width_blocks(occurrence) = owner
        lo_widths(occurrence) = decay_lo_width(pdg)
        nlo_widths(occurrence) = decay_nlo_width(pdg)
      end do
    end subroutine initialize_validation_widths


    subroutine extract_validation_block_orders( &
         selected_tuple, selected_distributions, evaluated_order, &
         block_orders)
      type(multiplicative_density_tuple), intent(in) :: selected_tuple
      type(block_nlo_distribution), intent(in) :: selected_distributions(:)
      integer, intent(in) :: evaluated_order
      integer, intent(out) :: block_orders(:)
      integer :: component, term

      if (size(block_orders) /= size(selected_distributions) .or. &
          selected_tuple%distribution_count /= &
          size(selected_distributions)) then
        call fail_driver('the validation block layout is inconsistent')
      end if
      do component = 1, size(selected_distributions)
        term = selected_tuple%term_indices(component)
        if (term < 1 .or. &
            term > selected_distributions(component)%term_count) then
          call fail_driver('a validation tuple term is out of range')
        end if
        block_orders(component) = &
             selected_distributions(component)%terms(term)%nlo_order
      end do
      if (sum(block_orders) /= selected_tuple%nlo_order .or. &
          sum(block_orders) /= evaluated_order) then
        call fail_driver('the validation NLO order does not close')
      end if
    end subroutine extract_validation_block_orders


    integer function multiplicative_plot_weight_count()
      integer :: dynamic_index, set_index

      multiplicative_plot_weight_count = 1
      if (do_rwgt_scale .or. do_rwgt_decay_scale) then
        do dynamic_index = 1, dyn_scale(0)
          multiplicative_plot_weight_count = &
               multiplicative_plot_weight_count + &
               fnlo_scale_point_count(dynamic_index)
        end do
      end if
      if (do_rwgt_pdf) then
        do set_index = 1, lhaPDFid(0)
          multiplicative_plot_weight_count = &
               multiplicative_plot_weight_count + nmemPDF(set_index) + 1
        end do
      end if
    end function multiplicative_plot_weight_count


    subroutine map_decay_factor_indices(species_factors, block_factors)
      integer, intent(in) :: species_factors(:)
      integer, intent(out) :: block_factors(0:nexternal)
      integer :: component, physical_block, pdg, species

      block_factors = 1
      if (size(species_factors) == 0) return
      if (size(species_factors) /= decay_scale_species_count()) then
        call fail_driver('the decay scale-factor vector has the wrong size')
      end if
      do component = 1, component_count
        physical_block = sdm_multiplicative_physical_block(component)
        if (physical_block == 0) cycle
        pdg = sdm_multiplicative_block_pdg(physical_block)
        do species = 1, size(species_factors)
          if (abs(pdg) /= abs(decay_scale_species(species))) cycle
          block_factors(physical_block) = species_factors(species)
          exit
        end do
      end do
    end subroutine map_decay_factor_indices


    subroutine require_real_multiplicative_weight(weight)
      complex(kind=8), intent(in) :: weight

      if (abs(aimag(weight)) > 1d-8*max(1d0,abs(dble(weight)))) then
        call fail_driver( &
             'a varied multiplicative density contraction is not real')
      end if
    end subroutine require_real_multiplicative_weight

  end subroutine sigint_multiplicative_impl

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
