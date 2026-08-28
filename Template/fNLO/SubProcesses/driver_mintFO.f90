module driver_mintfo_module
  use run_printout_module, only: write_run_summary
  use extra_weights, only: doreweight
  use mint_module, only: maxchannels, n_ave_virt, average_virtual, &
                         virtual_fraction, min_virt_fraction_mint, n_ord_virt, ncalls0, &
                         itmax, imode, ndim, ndimmax, nintegrals, nchans, &
                         iconfig, ichan, &
                         iconfigs, accuracy, wgt_mult, new_point, pass_cuts_check, &
                         virt_wgt_mint, born_wgt_mint, mint
  use mint_module, only: ans_result => ans, unc_result => unc
  use FKSParams, only: min_virt_fraction, virt_fraction, FKSParamReader
  use weight_lines, only: icontr, deallocate_weight_lines
  use process_dimensions, only: nexternal, nincoming, fks_configs, &
                                amp_split_size, lmaxconfigs, &
                                validate_process_and_born_dimensions
  use fks_metadata, only: fks_i_d, pdg_type_d, validate_fks_metadata
  use fks_channel_map, only: fks_channel_count, &
                             fks_channel_configuration, &
                             print_fks_channel_map, get_born_fks_process
  use fks_random_module, only: random_unit_interval
  use run_state, only: lpp, fixed_fac_scale, muf1_over_ref, &
                       muf2_over_ref, muf1_ref_fixed, muf2_ref_fixed, &
                       do_rwgt_scale, do_rwgt_pdf
  use genps_fks, only: generate_momenta
  use decay_chain_metadata, only: real_phase_space_dimension
  use nlo_contribution_bundle, only: has_nlo_contribution_bundle, &
       nlo_contribution_count, contribution_representative_fks, &
       nlo_virtual_grid_count
  use setscales_module, only: set_alphas
  use split_orders, only: check_amp_split
  use cuts_module, only: passcuts
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
  use fks_singular_module, only: fill_configurations_common, setfksfactor
  use madfks_plot_module, only: topout_impl
  use fnlo_process_common, only: nfksprocess, soft_counterevent, &
                                 collinear_counterevent, &
                                 real_event, &
                                 ybst_til_tolab, force_polecheck
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

    ndim = real_phase_space_dimension()
    if (abs(lpp(1)) >= 1) ndim = ndim + 1
    if (abs(lpp(2)) >= 1) ndim = ndim + 1
    nndim = ndim

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
        do_rwgt_pdf = .false.
      else
        doreweight = do_rwgt_scale .or. do_rwgt_pdf
      end if
      write (*, *) 'imode is ', imode
      call mint(sigint)
      call topout_impl()
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
    integer :: amplitude_order, picked_integer, position
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

    call update_vegas_x_impl(xx, vegas_variables, nndim, abrv)
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
        call update_fks_dir_impl(ifks)
        call generate_momenta(nndim, iconfig, jacobian, &
                              vegas_variables, momentum)
        if (p_born(0, 1) < 0d0) cycle

        call compute_prefactors_n1body(vegas_wgt, jacobian)
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
          call compute_real_emission(momentum)
        end if
      end do
    end if

    call include_pdf_and_alphas()
    if (doreweight) then
      if (do_rwgt_scale) call reweight_scale()
      if (do_rwgt_pdf) call reweight_pdf()
    end if

    call get_wgt_no_nbody(sampled_weight)
    call fill_mc_integer(max(ini_fin_fks(ichan), 1), picked_integer, &
                         abs(sampled_weight)*volume)

    call fill_plots()
    call fill_mint_function(f)
  end function sigint_impl

  subroutine update_fks_dir_impl(nfks)
    implicit none
    integer, intent(in) :: nfks

    nfksprocess = nfks
    call fks_inc_chooser()
    call leshouche_inc_chooser()
    call setcuts()
    call setfksfactor()
  end subroutine update_fks_dir_impl


  subroutine update_vegas_x_impl(xx, x, nndim, abrv)
    implicit none
    double precision, intent(in) :: xx(ndimmax)
    double precision, intent(out) :: x(99)
    integer, intent(in) :: nndim
    character(len=4), intent(in) :: abrv
    integer :: i, copied_dimensions

    x = 0d0
    copied_dimensions = min(nndim, ndimmax, 99)
    if (abrv == 'born' .or. abrv(1:2) == 'vi') then
      do i = 1, min(copied_dimensions, nndim - 3)
        x(i) = xx(i)
      end do
      do i = max(1, nndim - 2), min(nndim, 99)
        x(i) = random_unit_interval(iconfig)
      end do
    else
      do i = 1, copied_dimensions
        x(i) = xx(i)
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
