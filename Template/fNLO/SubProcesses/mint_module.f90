!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! MINT Integrator Package
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Original version by Paolo Nason (for POWHEG (BOX))
! Modified by Rikkert Frederix (for MadGraph5_aMC@NLO)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!      subroutine mint(fun,ndim,ncalls0,itmax,imode,
! ndim=number of dimensions
! ncalls0=# of calls per iteration
! itmax =# of iterations
! fun(xx,www,ifirst): returns the function to be integrated multiplied by www;
!                     xx(1:ndim) are the variables of integration
!                     ifirst=0: normal behaviour
! imode: integer flag
!
! imode=-1:
! same as imode=0 as far as this routine is concerned, except for the
! fact that a grid is read at the beginning (rather than initialized).
! The return value of imode will be zero.
!
! imode=0:
! When called with imode=0 the routine integrates the absolute value of
! the function and sets up a grid xgrid(0:50,ndim) such that in each
! ndim-1 dimensional slice (i.e. xgrid(m-1,n)<xx(n)<xgrid(m,n)) the
! contribution of the integral is the same the array xgrid is setup at
! this stage; ans and err are the integral and its error
!
! Added the posibility to keep track of more than one integral:
!
! nintegrals=1 : the function that is used to update the grids. This is
! the ABS cross section.
! nintegrals=2 : the actual cross section. This includes virtual corrections.
! nintegrals=3 : the cross section from the M_Virt/M_Born ratio alone:
! this defines the average virtual that is added to each phase-space
! point
! nintegrals=4 : the cross section of the actual virtual minus the
! average virtual. This is used to determine the fraction of phase-space
! points for which we include the virtual.
! nintegrals=5 : abs of 3
! nintegrals=6 : born
! nintegrals=7..6+2*n_ave_virt : virtual and born order by order
! remaining integrals: resolved Born/production/decay/width components
!

module mint_module
  use fnlo_process_common, only: fnlo_maxchannels, fks_configs, &
                                 amp_split_size
  use FKSParams ! contains use_poly_virtual
  use mc_integer_module, only: regrid_MC_integer, empty_MC_integer, &
                               reset_MC_grid
  use nlo_contribution_bundle, only: factorized_shared_dimension
  use factorized_mint_policy, only: factorized_mint_grid_weight, &
       factorized_mint_uses_uniform_channels, &
       factorized_mint_shows_multiplicative_validation
  use polynomial_fit, only: init_polyfit, add_point_polyfit, &
                            do_polyfit, get_polyfit, save_polyfit, restore_polyfit
  implicit none
  private
  integer, parameter, private :: nintervals = 32    ! max number of intervals in the integration grids
  ! The canonical phase-space vector remains bounded by 99 entries.  MINT's
  ! vector can be larger because every factorized NLO block owns an
  ! additional copy of the three FKS radiation variables.
  integer, parameter, public  :: ndimmax = 3*99     ! max number of integration-grid dimensions
  ! A contribution bundle can assign one virtual grid to every
  ! contribution/split-order pair.  Size the workspace from the generated
  ! process instead of imposing the historical ten-grid limit.  Keep ten as
  ! a floor so ordinary-process checkpoint layouts remain unchanged.
  integer, parameter, public  :: n_ave_virt = &
       max(10, fks_configs*amp_split_size)
  integer, parameter, public  :: max_bundle_components = fks_configs + 2
  integer, parameter, public  :: first_bundle_component_integral = &
       7 + 2*n_ave_virt
  integer, parameter, public :: first_multiplicative_linear_integral = &
       first_bundle_component_integral + max_bundle_components
  integer, parameter, public :: first_multiplicative_group_integral = &
       first_multiplicative_linear_integral + fks_configs
  ! The multiplicative path does not use the two legacy virtual-ratio
  ! diagnostics.  Reuse only those inert slots for formal-lambda validation;
  ! slots 3 and 6 participate in virtual-grid control and must remain zero.
  ! Existing MINT checkpoint dimensions remain byte-for-byte compatible.
  integer, parameter, public :: multiplicative_lo_integral = &
       4
  integer, parameter, public :: multiplicative_additive_integral = &
       5
  integer, parameter, public  :: nintegrals = &
       6 + 2*n_ave_virt + max_bundle_components + 4*fks_configs
  integer, parameter, private :: nintervals_virt = 8! max number of intervals in the grids for the approx virtual
  integer, parameter, private :: min_inter = 4      ! minimal number of intervals
  integer, parameter, private :: min_it0 = 4        ! minimal number of iterations in the mint step 0 phase
  integer, parameter, private :: mint_grid_format_version = 2
  character(len=*), parameter, private :: mint_grid_magic = 'MINT_GRID'
! Maximum points to try per iteration when too few non-zero points are found.
  integer, parameter, private :: max_points = 100000
  integer, parameter, public  :: maxchannels = fnlo_maxchannels
  ! Note that the number of intervals in the integration grids, 'nintervals', cannot be arbitrarily large.
  ! It should be equal to
  !     nintervals = min_inter * 2^n,
  ! where 'n' is an integer smaller than or equal to min_it0.
  !
  ! The number of intergrals should be equal to
  !     nintegrals=6+2*n_ave_virt+max_bundle_components+4*fks_configs
  !

! public variables
  integer, public :: ncalls0, ndim, itmax, imode, n_ord_virt, nchans, iconfig, ichan
  integer, dimension(maxchannels), public :: iconfigs
  double precision, public :: accuracy, min_virt_fraction_mint, wgt_mult
  double precision, dimension(0:n_ave_virt, maxchannels), public :: average_virtual
  double precision, dimension(0:n_ave_virt), public :: virt_wgt_mint, born_wgt_mint, polyfit
  double precision, dimension(maxchannels), public :: virtual_fraction
  double precision, dimension(nintegrals, 0:maxchannels), public :: ans, unc
  logical, public :: new_point, pass_cuts_check

! private variables
  integer, private :: nit, nit_included, kpoint_iter, nint_used, nint_used_virt, min_it, ncalls, pass_cuts_point, ng, npg, k
  integer, private :: born_dimensions = 0
  integer, allocatable, private :: icell(:), ncell(:)
  integer, dimension(nintegrals), private :: non_zero_point, ntotcalls
  integer, allocatable, private :: nhits(:, :, :)
  integer, allocatable, private :: nhits_in_grids(:)
  integer, allocatable, private :: nvirt(:, :, :, :), nvirt_acc(:, :, :, :)
  logical, private :: double_points, reset, firsttime
  logical, allocatable, private :: regridded(:)
  double precision, allocatable, private :: xgrid(:, :, :), xacc(:, :, :)
  double precision, allocatable, private :: vtot(:, :), etot(:, :), chi2(:, :)
  double precision, dimension(nintegrals, 3), private :: ans3, unc3
  double precision, dimension(nintegrals), private :: ans_l3, unc_l3, chi2_l3, f
  double precision, allocatable, private :: ans_chan(:)
  double precision, dimension(2), private :: HwU_values
  double precision, allocatable, private :: ave_virt(:, :, :, :), &
    ave_virt_acc(:, :, :, :), ave_born_acc(:, :, :, :)
  double precision, private :: vol_chan
  double precision, allocatable, private :: rand(:)
  double precision, allocatable, private :: xgrid_new(:, :)
  logical, private :: bad_iteration = .false.
  logical, private :: mint_state_initialized = .false.
  double precision, private :: even_dng = 0d0
  integer, private :: even_current_dim = 0
  integer, allocatable, private :: even_iii(:), even_kkk(:)

! functions and subroutines:
  public :: mint
  private :: initialise_mint, setup_basic_mint &
       &, update_accumulated_results, prepare_next_iteration &
       &, check_desired_accuracy, update_integration_grids &
       &, combine_final_three_iterations &
       &, print_accumulated_last_three &
       &, update_virtual_fraction, combine_iterations &
       &, print_results_accumulated, check_fractional_uncertainty &
       &, print_results_current_iteration &
       &, compute_fractional_uncertainty, combine_results_channels &
       &, check_for_special_channels_loop &
       &, combine_special_channels, get_amount_of_points &
       &, add_point_to_grids &
       &, accumulate_the_point, compute_integrand, get_random_x &
       &, start_iteration, reset_accumulated_grids &
       &, check_evenly_random_numbers, finalise_mint, write_results &
       &, write_channel_info, setup_imode_m1, setup_imode_0 &
       &, reset_mint_grids, setup_common, write_grids_to_file &
       &, double_grid, regrid, smooth_xacc, nextlexi, init_ave_virt&
       &, get_ave_virt, fill_ave_virt, regrid_ave_virt, double_ave_virt&
       &, get_channel, close_run_zero_res &
       &, initialize_even_random_numbers, get_ran
contains

  character(len=13) function integral_title(index)
    implicit none
    integer, intent(in) :: index

    select case (index)
    case (1)
      integral_title = 'ABS integral '
    case (2)
      integral_title = 'Integral     '
    case (3)
      integral_title = 'Virtual      '
    case (4)
      if (show_multiplicative_validation_titles()) then
        integral_title = 'Mult LO      '
      else
        integral_title = 'Virtual ratio'
      end if
    case (5)
      if (show_multiplicative_validation_titles()) then
        integral_title = 'Mult additive'
      else
        integral_title = 'ABS virtual  '
      end if
    case (6)
      integral_title = 'Born         '
    case default
      if (index >= first_multiplicative_group_integral) then
        write (integral_title, '(a1,i2,a1,i1,8x)') &
             'G', (index - first_multiplicative_group_integral)/3 + 1, &
             '.', mod(index - first_multiplicative_group_integral, 3) + 1
      else if (index >= first_multiplicative_linear_integral) then
        write (integral_title, '(a1,i3,9x)') &
             'L', index - first_multiplicative_linear_integral + 1
      else if (index >= first_bundle_component_integral) then
        write (integral_title, '(a1,i3,9x)') &
             'C', index - first_bundle_component_integral + 1
      else if (mod(index, 2) == 1) then
        write (integral_title, '(a1,i3,9x)') 'V', (index - 5)/2
      else
        write (integral_title, '(a1,i3,9x)') 'B', (index - 6)/2
      end if
    end select
  end function integral_title


  logical function show_multiplicative_validation_titles()
    show_multiplicative_validation_titles = &
         factorized_mint_shows_multiplicative_validation()
  end function show_multiplicative_validation_titles

  subroutine initialize_mint_state
    implicit none
    integer :: expected_born_dimensions

    if (ndim .lt. 1 .or. ndim .gt. ndimmax) then
      write (*, *) 'Invalid MINT dimension:', ndim
      stop 1
    end if
    if (nchans .lt. 1 .or. nchans .gt. maxchannels) then
      write (*, *) 'Invalid number of MINT channels:', nchans
      stop 1
    end if
    if (n_ord_virt .lt. 0 .or. n_ord_virt .gt. n_ave_virt) then
      write (*, *) 'Invalid number of virtual orders:', n_ord_virt
      stop 1
    end if

    expected_born_dimensions = factorized_shared_dimension(ndim)
    if (mint_state_initialized) then
      if (size(icell) .eq. ndim .and. size(nhits, 3) .eq. nchans .and. &
          size(nvirt, 2) .eq. expected_born_dimensions .and. &
          ubound(nvirt, 3) .eq. n_ord_virt) return
      call finalize_mint_state
    end if

    born_dimensions = expected_born_dimensions

    allocate (icell(ndim), ncell(ndim), rand(ndim))
    allocate (nhits(nintervals, ndim, nchans))
    allocate (nhits_in_grids(nchans), regridded(nchans))
    allocate (xgrid(0:nintervals, ndim, nchans))
    allocate (xacc(0:nintervals, ndim, nchans))
    allocate (vtot(nintegrals, 0:nchans))
    allocate (etot(nintegrals, 0:nchans))
    allocate (chi2(nintegrals, 0:nchans))
    allocate (ans_chan(0:nchans))
    allocate (xgrid_new(0:nintervals, ndim))
    allocate (nvirt(nintervals_virt, born_dimensions, &
                    0:n_ord_virt, nchans))
    allocate (nvirt_acc(nintervals_virt, born_dimensions, &
                        0:n_ord_virt, nchans))
    allocate (ave_virt(nintervals_virt, born_dimensions, &
                       0:n_ord_virt, nchans))
    allocate (ave_virt_acc(nintervals_virt, born_dimensions, &
                           0:n_ord_virt, nchans))
    allocate (ave_born_acc(nintervals_virt, born_dimensions, &
                           0:n_ord_virt, nchans))
    allocate (even_iii(ndim), even_kkk(ndim))

    icell = 0
    ncell = 0
    rand = 0d0
    nhits = 0
    nhits_in_grids = 0
    regridded = .false.
    xgrid = 0d0
    xacc = 0d0
    vtot = 0d0
    etot = 0d0
    chi2 = 0d0
    ans_chan = 0d0
    xgrid_new = 0d0
    nvirt = 0
    nvirt_acc = 0
    ave_virt = 0d0
    ave_virt_acc = 0d0
    ave_born_acc = 0d0
    even_iii = 1
    even_kkk = 1
    even_dng = 0d0
    even_current_dim = 0
    bad_iteration = .false.
    mint_state_initialized = .true.
  end subroutine initialize_mint_state

  subroutine finalize_mint_state
    implicit none

    if (allocated(icell)) deallocate (icell)
    if (allocated(ncell)) deallocate (ncell)
    if (allocated(rand)) deallocate (rand)
    if (allocated(nhits)) deallocate (nhits)
    if (allocated(nhits_in_grids)) deallocate (nhits_in_grids)
    if (allocated(regridded)) deallocate (regridded)
    if (allocated(xgrid)) deallocate (xgrid)
    if (allocated(xacc)) deallocate (xacc)
    if (allocated(vtot)) deallocate (vtot)
    if (allocated(etot)) deallocate (etot)
    if (allocated(chi2)) deallocate (chi2)
    if (allocated(ans_chan)) deallocate (ans_chan)
    if (allocated(xgrid_new)) deallocate (xgrid_new)
    if (allocated(nvirt)) deallocate (nvirt)
    if (allocated(nvirt_acc)) deallocate (nvirt_acc)
    if (allocated(ave_virt)) deallocate (ave_virt)
    if (allocated(ave_virt_acc)) deallocate (ave_virt_acc)
    if (allocated(ave_born_acc)) deallocate (ave_born_acc)
    if (allocated(even_iii)) deallocate (even_iii)
    if (allocated(even_kkk)) deallocate (even_kkk)
    even_dng = 0d0
    even_current_dim = 0
    born_dimensions = 0
    bad_iteration = .false.
    mint_state_initialized = .false.
  end subroutine finalize_mint_state

  subroutine mint(fun)
    implicit none
    integer kpoint
    double precision :: vol
    double precision, dimension(ndimmax) :: x
    double precision, external :: fun
    logical :: enough_points, channel_loop_done
    call initialise_mint
    do while (nit .lt. itmax)
      call start_iteration
2     kpoint_iter = kpoint_iter + 1
      do kpoint = 1, ncalls
        new_point = .true.
        call get_random_x(x, vol)
        call compute_integrand(fun, x, vol)
        call accumulate_the_point(x)
      end do
      call get_amount_of_points(enough_points)
      if (.not. enough_points) goto 2
      if (imode .eq. 0 .and. nit .eq. 1 .and. double_points) then
        call check_for_special_channels_loop(channel_loop_done)
        if (.not. channel_loop_done) goto 2
        call combine_special_channels
      else
        call combine_results_channels
      end if
      call update_accumulated_results
    end do
    call finalise_mint
  end subroutine mint

  subroutine initialise_mint
    implicit none
    call initialize_mint_state
    if (imode .ne. 0) call read_grids_from_file
    call setup_basic_mint
    if (imode .eq. 0) then
      call setup_imode_0
    elseif (imode .eq. -1) then
      call setup_imode_m1
    else
      write (*, *) 'Unsupported MINT mode in the fixed-order template:', imode
      stop 1
    end if
    call setup_common
  end subroutine initialise_mint

  subroutine setup_basic_mint
    implicit none
    ! if ncalls0 is greater than 0, use the default running, i.e. do not
    ! double the integration points after each iteration and use a fixed number
    ! of intervals in the grids.
    if (ncalls0 .gt. 0) then
      double_points = .false.
      nint_used = nintervals
      nint_used_virt = nintervals_virt
    else
      ! if ncalls0.le.0, reset it and double the integration points per iteration
      ncalls0 = 80*ndim*(nchans/3 + 1)
      double_points = .true.
      if (imode .eq. -1) then
        nint_used = nintervals
        nint_used_virt = nintervals_virt
      else
        nint_used = min_inter
        nint_used_virt = min_inter
      end if
    end if
    reset = .false.
    ncalls = 0  ! # PS points (updated below)
  end subroutine setup_basic_mint

  subroutine update_accumulated_results
    implicit none
    double precision, dimension(nintegrals) :: efrac
    logical :: iterations_done
    call compute_fractional_uncertainty(efrac)
    call print_results_current_iteration(efrac)
    call check_fractional_uncertainty(efrac)
    if (reset) return ! iteration was not accurate enough: do not include it
    call combine_iterations
    call combine_final_three_iterations
    call HwU_accum_iter(.true., ntotcalls(1), HwU_values)
    if (imode .eq. 0) then
      call update_virtual_fraction
      call update_integration_grids
    end if
    call check_desired_accuracy(iterations_done)
    if (.not. iterations_done) then
      call prepare_next_iteration
    else
      nit = itmax
    end if
  end subroutine update_accumulated_results

  subroutine prepare_next_iteration
    implicit none
    integer :: kchan, kdim, k_ord_virt
    if (double_points) then
! Double the number of intervals in the grids if not yet reach the maximum
      if (2*nint_used .le. nintervals) then
        do kchan = 1, nchans
          do kdim = 1, ndim
            call double_grid(kdim, kchan)
          end do
        end do
        nint_used = 2*nint_used
      end if
      if (2*nint_used_virt .le. nintervals_virt) then
        do k_ord_virt = 0, n_ord_virt
          call double_ave_virt(k_ord_virt)
        end do
        nint_used_virt = 2*nint_used_virt
      end if
! double the number of points for the next iteration
      ncalls0 = ncalls0*2
    end if
  end subroutine prepare_next_iteration

  subroutine check_desired_accuracy(iterations_done)
    implicit none
    logical :: iterations_done
    integer :: i
! Quit if the desired accuracy has been reached
    iterations_done = .false.
    if (nit_included .ge. min_it .and. accuracy .gt. 0d0) then
      if (unc(1, 0)/ans(1, 0)*max(1d0, chi2(1, 0)/dble(nit_included - 1)) .lt. accuracy) then
        write (*, *) 'Found desired accuracy'
        iterations_done = .true.
      elseif (unc_l3(1)/ans_l3(1)*max(1d0, chi2_l3(1)) .lt. accuracy) then
        write (*, *) 'Found desired accuracy in last 3 iterations'
        iterations_done = .true.
        ! overwrite results with the results from the last three iterations
        do i = 1, nintegrals
          ans(i, 0) = ans_l3(i)
          unc(i, 0) = unc_l3(i)
          chi2(i, 0) = chi2_l3(i)*dble(nit_included - 1)
        end do
      end if
    end if
  end subroutine check_desired_accuracy

  subroutine update_integration_grids
    implicit none
    integer :: kchan, kdim, k_ord_virt
    do kchan = 1, nchans
      do kdim = 1, ndim
        call regrid(kdim, kchan)
      end do
      ! overwrite xgrid with the new xgrid
      if (regridded(kchan)) xgrid(1:nint_used, 1:ndim, kchan) = xgrid_new(1:nint_used, 1:ndim)
    end do
    if (use_poly_virtual) then
      call do_polyfit()
    else
      do k_ord_virt = 0, n_ord_virt
        call regrid_ave_virt(k_ord_virt)
      end do
    end if
! Regrid the MC over integers (used for the MC over FKS dirs)
    call regrid_MC_integer
  end subroutine update_integration_grids

  subroutine combine_final_three_iterations
    implicit none
    integer :: i, j
! Update the results of the last tree iterations
    do j = 1, 2
      ans3(1:nintegrals, j) = ans3(1:nintegrals, j + 1)
      unc3(1:nintegrals, j) = unc3(1:nintegrals, j + 1)
    end do
    ans3(1:nintegrals, 3) = vtot(1:nintegrals, 0)
    unc3(1:nintegrals, 3) = etot(1:nintegrals, 0)
! Compute the results of the last three iterations
    if (nit_included .ge. 4) then
      do i = 1, nintegrals
        ans_l3(i) = 0d0
        unc_l3(i) = ans3(i, 1)*1d99
        chi2_l3(i) = 0d0
        do j = 1, 3 ! the three final iterations
          if (i .ne. 1 .and. (unc_l3(i) .eq. 0d0 .or. unc3(i, j) .eq. 0d0)) then
            continue ! do not do anything
          else
            ans_l3(i) = (ans_l3(i)/unc_l3(i) + ans3(i, j)/unc3(i, j))/(1d0/unc_l3(i) + 1d0/unc3(i, j))
            unc_l3(i) = 1d0/sqrt(1d0/unc_l3(i)**2 + 1/unc3(i, j)**2)
            chi2_l3(i) = chi2_l3(i) + (ans3(i, j) - ans_l3(i))**2/unc3(i, j)**2
          end if
        end do
        chi2_l3(i) = chi2_l3(i)/2d0 ! three iterations, so 2 degrees of freedom
      end do
      call print_accumulated_last_three
    end if
  end subroutine combine_final_three_iterations

  subroutine print_accumulated_last_three
    implicit none
    integer :: i
    double precision, dimension(nintegrals) :: efrac
    do i = 1, 2
      if (ans_l3(i) .ne. 0d0) then
        efrac(i) = abs(unc_l3(i)/ans_l3(i))
      else
        efrac(i) = 0d0
      end if
      if (ans_l3(i) .ne. 0d0 .and. unc_l3(i) .ne. 0d0) then
        write (*, '(a,1x,e10.4,1x,a,1x,e10.4,1x,a,1x,f7.3,1x,a)') &
          'accumulated results last 3 iterations '//integral_title(i)//' =', &
          ans_l3(i), ' +/- ', unc_l3(i), ' (', efrac(i)*100d0, '%)'
      end if
    end do
    write (*, '(a,1x,e10.4)') 'accumulated result last 3 iterrations Chi^2 per DoF =' &
      , chi2_l3(1)
  end subroutine print_accumulated_last_three

  subroutine update_virtual_fraction
! Update the fraction of phase-space points for which we include the virtual corrections
! in the calculation
    implicit none
    integer kchan, k_ord_virt
    double precision :: error_virt
    do kchan = 1, nchans
      error_virt = 0d0
      do k_ord_virt = 1, n_ord_virt
        error_virt = error_virt + etot(2*k_ord_virt + 5, kchan)**2
      end do
      error_virt = sqrt(error_virt)
      virtual_fraction(kchan) = max(min(virtual_fraction(kchan) &
                                        *max(min(2d0*error_virt/etot(1, kchan), 2d0), 0.25d0), 1d0) &
                                    , Min_virt_fraction_mint)
    end do
  end subroutine update_virtual_fraction

  subroutine combine_iterations
    implicit none
    integer i, kchan
    HwU_values(1) = etot(1, 0)
    HwU_values(2) = unc(1, 0)
    if (nit .eq. 1) then ! first iteration
      ans(1:nintegrals, 0:nchans) = vtot(1:nintegrals, 0:nchans)
      unc(1:nintegrals, 0:nchans) = etot(1:nintegrals, 0:nchans)
      ans_chan(0:nchans) = ans(1, 0:nchans)
      write (*, '(a,1x,e10.4)') 'Chi^2 per d.o.f.', 0d0
    else
      do kchan = nchans, 0, -1 ! go backwards so that kchan=0 goes last
        ! (this makes sure central value is correctly updated).
        do i = 1, nintegrals
          if (i .ne. 1 .and. (etot(i, 0) .eq. 0d0 .or. unc(i, 0) .eq. 0d0)) then
            continue ! do not do anything
          else
            ! Use one common inverse-variance weight for every channel so
            ! their accumulated central values still add to the inclusive
            ! result.  The previous 1/sigma weighting was inconsistent with
            ! the inverse-variance uncertainty immediately below.
            ans(i, kchan) = &
                 (ans(i, kchan)/unc(i, 0)**2 + &
                  vtot(i, kchan)/etot(i, 0)**2)/ &
                 (1d0/unc(i, 0)**2 + 1d0/etot(i, 0)**2)
            unc(i, kchan) = 1d0/sqrt(1d0/unc(i, kchan)**2 + 1d0/etot(i, kchan)**2)
            chi2(i, kchan) = chi2(i, kchan) + (vtot(i, kchan) - ans(i, kchan))**2/etot(i, kchan)**2
          end if
        end do
        ans_chan(kchan) = ans(1, kchan)
      end do
      write (*, '(a,1x,e10.4)') 'Chi^2=', (vtot(1, 0) - ans(1, 0))**2/etot(1, 0)**2
    end if
    nit_included = nit_included + 1
    call print_results_accumulated
  end subroutine combine_iterations

  subroutine print_results_accumulated
    implicit none
    integer i
    double precision, dimension(nintegrals) :: efrac
    do i = 1, nintegrals
      if (ans(i, 0) .ne. 0d0) then
        efrac(i) = abs(unc(i, 0)/ans(i, 0))
      else
        efrac(i) = 0d0
      end if
      if (ans(i, 0) .ne. 0d0 .and. unc(i, 0) .ne. 0d0) then
        write (*, '(a,1x,e10.4,1x,a,1x,e10.4,1x,a,1x,f7.3,1x,a)') &
          'accumulated results '//integral_title(i)//' =', ans(i, 0), &
          ' +/- ', unc(i, 0), ' (', efrac(i)*100d0, '%)'
      end if
    end do
    if (nit_included .le. 1) then
      write (*, '(a,1x,e10.4)') 'accumulated result Chi^2 per DoF =', 0d0
    else
      write (*, '(a,1x,e10.4)') 'accumulated result Chi^2 per DoF =', chi2(1, 0)/dble(nit_included - 1)
    end if
  end subroutine print_results_accumulated

  subroutine check_fractional_uncertainty(efrac)
    implicit none
    double precision, dimension(nintegrals) :: efrac
! If there was a large fluctation in this iteration, be careful with
! including it in the accumalated results and plots.
    if (efrac(1) .gt. 0.3d0 .and. nit .gt. 3) then
! Do not include the results in the plots
      call HwU_accum_iter(.false., ntotcalls(1), HwU_values)
! Do not include the results in the updating of the grids.
      write (*, *) 'Large fluctuation ( >30 % ). Not including iteration in results.'
! empty the accumulated results in the MC over integers
      call empty_MC_integer
! empty the accumulated results for the MINT grids (Cannot really
! skip the increase of the upper bounding envelope. So, simply
! continue here. Note that no matter how large the integrand for the
! PS point, the upper bounding envelope is at most increased by a
! factor 2, so this should be fine).
      reset = .true.
! double the number of points for the next iteration
      if (double_points) ncalls0 = ncalls0*2
      if (bad_iteration .and. imode .eq. 0 .and. double_points) then
! 2nd bad iteration is a row. Reset grids
        write (*, *) '2nd bad iteration in a row. Resetting grids and starting from scratch...'
        if (double_points) then
          if (imode .eq. 0) nint_used = min_inter ! reset number of intervals
          ncalls0 = ncalls0/8   ! Start with larger number
        end if
        call reset_mint_grids
        call reset_MC_grid  ! reset the grid for the integers
        call initplot  ! Also reset all the plots
        call setup_common
        bad_iteration = .false.
      else
        bad_iteration = .true.
      end if
    else
      bad_iteration = .false.
    end if
  end subroutine check_fractional_uncertainty

  subroutine print_results_current_iteration(efrac)
    implicit none
    integer :: i
    double precision, dimension(nintegrals) :: efrac
    do i = 1, nintegrals
      if (vtot(i, 0) .ne. 0d0 .and. etot(i, 0) .ne. 0d0) then
        write (*, '(a,1x,e10.4,1x,a,1x,e10.4,1x,a,1x,f7.3,1x,a)') &
          integral_title(i)//' =', vtot(i, 0), ' +/- ', etot(i, 0), &
          ' (', efrac(i)*100d0, '%)'
      end if
    end do
  end subroutine print_results_current_iteration

  subroutine compute_fractional_uncertainty(efrac)
    implicit none
    integer :: i
    double precision, dimension(nintegrals) :: efrac
    do i = 1, nintegrals
      if (vtot(i, 0) .ne. 0d0) then
        efrac(i) = abs(etot(i, 0)/vtot(i, 0))
      else
        efrac(i) = 0d0
      end if
    end do
  end subroutine compute_fractional_uncertainty

  subroutine combine_results_channels
    implicit none
    integer :: kchan
    vtot(1:nintegrals, 0) = sum(vtot(1:nintegrals, 1:nchans), dim=2)
    etot(1:nintegrals, 0) = sum(etot(1:nintegrals, 1:nchans), dim=2)
    do kchan = 0, nchans
      vtot(1:nintegrals, kchan) = vtot(1:nintegrals, kchan)/dble(ntotcalls(1:nintegrals))
      etot(1:nintegrals, kchan) = etot(1:nintegrals, kchan)/dble(ntotcalls(1:nintegrals))
      etot(1:nintegrals, kchan) = sqrt(abs(etot(1:nintegrals, kchan) - vtot(1:nintegrals, kchan)**2) &
                                       /dble(ntotcalls(1:nintegrals)))
    end do
  end subroutine combine_results_channels

  subroutine check_for_special_channels_loop(channel_loop_done)
    implicit none
    logical :: channel_loop_done
    integer :: kchan
    do kchan = nchans, 1, -1
      if (ans_chan(kchan) .eq. 1d0) then
! results of the current channel
        vtot(1:nintegrals, kchan) = vtot(1:nintegrals, kchan)/dble(ntotcalls(1:nintegrals))
        etot(1:nintegrals, kchan) = etot(1:nintegrals, kchan)/dble(ntotcalls(1:nintegrals))
        etot(1:nintegrals, kchan) = sqrt(abs(etot(1:nintegrals, kchan) - vtot(1:nintegrals, kchan)**2) &
                                         /dble(ntotcalls(1:nintegrals)))
        if (kchan .eq. nchans) then
! done all channels
          channel_loop_done = .true.
          return
        end if
! prepare for the next channel
        ans_chan(kchan) = 0d0
        ans_chan(kchan + 1) = 1d0
        ntotcalls(1:nintegrals) = 0
        non_zero_point(1:nintegrals) = 0
        pass_cuts_point = 0
        kpoint_iter = 0
        channel_loop_done = .false.
        return
      end if
    end do
  end subroutine check_for_special_channels_loop

  subroutine combine_special_channels
    implicit none
! set the total result for the first iteration to the sum over all the channels
    vtot(1:nintegrals, 0) = sum(vtot(1:nintegrals, 1:nchans), dim=2)
    etot(1:nintegrals, 0) = sum(etot(1:nintegrals, 1:nchans)**2, dim=2)
    etot(1:nintegrals, 0) = sqrt(etot(1:nintegrals, 0))
    ncalls0 = ncalls0*nchans
  end subroutine combine_special_channels

  subroutine get_amount_of_points(enough_points)
    ! fill the ntotcalls() array with the total number of calls used
    ! and check if this is enough for this iteration.
    implicit none
    logical :: enough_points
    integer :: i
    do i = 1, nintegrals
! Number of phase-space points used
      ntotcalls(i) = ncalls*kpoint_iter
! Special for the computation of the 'computed virtual'
      if (i .eq. 4 .and. non_zero_point(i) .ne. 0) &
        ntotcalls(i) = non_zero_point(i)
    end do

    if (.not. double_points) then
! If not doubling the number of integration points for each iteration, nothing
! needs to be done here.
      enough_points = .true.
      return
    end if
    if (pass_cuts_point .lt. 25) then
! Not enough points have passed to cuts to get a reliable estimate
      if (ntotcalls(1) .gt. max_points) then
! tried many points already. Need to crash.
        write (*, *) 'ERROR: NOT ENOUGH POINTS PASS THE CUTS. '// &
          'RESULTS CANNOT BE TRUSTED. '// &
          'LOOSEN THE GENERATION CUTS, OR ADAPT SET_TAU_MIN()'// &
          ' IN SETCUTS.F ACCORDINGLY.'
        stop 1
      else
        enough_points = .false.
        return
      end if
    end if
    if (non_zero_point(1) .lt. int(0.99*ncalls)) then
! Not enough (non-zero) points have been generated
      if (pass_cuts_point .gt. ncalls .and. &
          non_zero_point(1) .lt. 2) then
! Many points passed the cuts, but less than 2 non-zero integrand
! values: must be that the PDFs or the matrix elements (e.g. coupling
! constants) are numerically zero. End the run gracefully
        if (nit .gt. 1 .or. imode .ne. 0) then
          write (*, *) 'THE INTEGRAL APPEARS TO BE ZERO: END THE RUN GRACEFULLY.'
          write (*, *) 'TRIED', ntotcalls(1), 'PS POINTS AND ONLY ' &
            , non_zero_point(1), ' GAVE A NON-ZERO INTEGRAND.'
          call close_run_zero_res
          stop 0
        else
! This is for the special channels loop. Simply assume that the result
! for this channel is zero, and go to the next channel. If all
! channels give a zero result, end the run gracefully.
          vtot(1, ichan) = 0d0
          if (ichan .eq. nchans .and. all(vtot(1, 1:nchans) .eq. 0d0)) then
            write (*, *) 'THE INTEGRAL APPEARS TO BE ZERO: END THE RUN GRACEFULLY.'
            write (*, *) 'TRIED', ntotcalls(1), 'PS POINTS AND ONLY ' &
              , non_zero_point(1), ' GAVE A NON-ZERO INTEGRAND.'
            call close_run_zero_res
            stop 0
          end if
          enough_points = .true.
          return
        end if
      else
        if (ntotcalls(1) .lt. max_points) then
          enough_points = .false.
          return
        end if
      end if
    end if
    enough_points = .true.
  end subroutine get_amount_of_points

  subroutine add_point_to_grids(x)
    implicit none
    integer :: kdim, k_ord_virt, ithree, isix
    double precision, dimension(ndimmax) :: x
    double precision :: virtual, born, grid_weight
! accumulate the function in xacc(icell(kdim),kdim) to adjust the grid later
    do kdim = 1, ndim
      call factorized_mint_grid_weight( &
           ndim, kdim, f, multiplicative_lo_integral, &
           first_bundle_component_integral, grid_weight)
      xacc(icell(kdim), kdim, ichan) = &
           xacc(icell(kdim), kdim, ichan) + grid_weight
    end do
! Set the Born contribution (to compute the average_virtual) to zero if
! the virtual was not computed for this phase-space point. Compensate by
! including the virtual_fraction.
    do k_ord_virt = 0, n_ord_virt
      if (k_ord_virt .eq. 0) then
        ithree = 3
        isix = 6
      else
        ithree = 2*k_ord_virt + 5
        isix = 2*k_ord_virt + 6
      end if
      if (f(ithree) .ne. 0d0) then
        born = f(isix)
        ! virt_wgt_mint=(virtual-average_virtual*born)/virtual_fraction. Compensate:
        if (use_poly_virtual) then
          virtual = f(ithree)*virtual_fraction(ichan) + &
                    polyfit(k_ord_virt)*f(isix)
          call add_point_polyfit(ichan, k_ord_virt, &
                                 x(1:born_dimensions), &
                                 virtual/born, born/wgt_mult)
        else
          virtual = f(ithree)*virtual_fraction(ichan) + &
                    average_virtual(k_ord_virt, ichan)*f(isix)
          call fill_ave_virt(x, k_ord_virt, virtual, born)
        end if
      else
        f(isix) = 0d0
      end if
    end do
  end subroutine add_point_to_grids

  subroutine accumulate_the_point(x)
    implicit none
    integer :: i
    double precision, dimension(ndimmax) :: x
    call add_point_to_grids(x)
    do i = 1, nintegrals
      if (f(i) .ne. 0d0) non_zero_point(i) = non_zero_point(i) + 1
    end do
    if (pass_cuts_check) pass_cuts_point = pass_cuts_point + 1
! Add the PS point to the result of this iteration
    vtot(1:nintegrals, ichan) = vtot(1:nintegrals, ichan) + f(1:nintegrals)
    etot(1:nintegrals, ichan) = etot(1:nintegrals, ichan) + f(1:nintegrals)**2
! Accumulate the points in the HwU histograms
    if (f(1) .ne. 0d0) call HwU_add_points
  end subroutine accumulate_the_point

  subroutine compute_integrand(fun, x, vol)
    implicit none
    integer :: ifirst
    double precision :: dummy, vol
    double precision, dimension(nintegrals) :: f1
    double precision, dimension(ndimmax) :: x
    double precision, external :: fun
    ! contribution to integral
    ifirst = 0
    dummy = fun(x, vol, ifirst, f1)
    f(1:nintegrals) = f1(1:nintegrals)
  end subroutine compute_integrand

  subroutine get_random_x(x, vol)
    implicit none
    integer :: kdim, k_ord_virt
    double precision :: vol, dx
    double precision, dimension(ndimmax) :: x
    call get_channel
! find random x, and its random cell
    do kdim = 1, ndim
      rand(kdim) = get_ran()
      ncell(kdim) = min(int(rand(kdim)*nint_used) + 1, nint_used)
      rand(kdim) = rand(kdim)*nint_used - (ncell(kdim) - 1)
    end do
    vol = 1d0/vol_chan*wgt_mult
! convert 'flat x' ('rand') to 'vegas x' ('x') and include jacobian ('vol')
    do kdim = 1, ndim
      icell(kdim) = ncell(kdim)
      dx = xgrid(icell(kdim), kdim, ichan) - xgrid(icell(kdim) - 1, kdim, ichan)
      vol = vol*dx*nint_used
      x(kdim) = xgrid(icell(kdim) - 1, kdim, ichan) + rand(kdim)*dx
      nhits(icell(kdim), kdim, ichan) = nhits(icell(kdim), kdim, ichan) + 1
    end do
    do k_ord_virt = 0, n_ord_virt
      if (use_poly_virtual) then
        call get_polyfit(ichan, k_ord_virt, &
                         x(1:born_dimensions), polyfit(k_ord_virt))
      else
        call get_ave_virt(x, k_ord_virt)
      end if
    end do
  end subroutine get_random_x

  subroutine start_iteration
    implicit none
    call write_channel_info
    nit = nit + 1
    write (*, *) '------- iteration', nit
    call check_evenly_random_numbers
    if (imode .eq. 0) then
      call reset_accumulated_grids
    end if
    vtot(1:nintegrals, 0:nchans) = 0d0
    etot(1:nintegrals, 0:nchans) = 0d0
    kpoint_iter = 0
    non_zero_point(1:nintegrals) = 0
    pass_cuts_point = 0
  end subroutine start_iteration

  subroutine reset_accumulated_grids
    implicit none
    integer :: kchan
    do kchan = 1, nchans
      ! only reset if grids were updated (or there is a forced reset)
      if (regridded(kchan) .or. reset) then
        if (regridded(kchan) .and. .not. reset) then
          ! set nhits_in_grids equal to the number of points used for the last update
          nhits_in_grids(kchan) = sum(nhits(1:nint_used, 1, kchan), dim=1)
        elseif (regridded(kchan) .and. reset) then
          nhits_in_grids(kchan) = 0
        end if
        xacc(0:nint_used, 1:ndim, kchan) = 0d0
        nhits(1:nint_used, 1:ndim, kchan) = 0
      end if
    end do
    reset = .false.
  end subroutine reset_accumulated_grids

  subroutine check_evenly_random_numbers
    implicit none
    if (ncalls .ne. ncalls0) then
      ! Uses more evenly distributed random numbers. This overwrites
      ! the number of calls
      call initialize_even_random_numbers
      write (*, *) 'Update # PS points: ', ncalls0, ' --> ', ncalls
    end if
  end subroutine check_evenly_random_numbers

  subroutine finalise_mint
    implicit none
    integer :: kchan
    call write_channel_info
    if (nit_included .ge. 2) then
      chi2(1, 0:nchans) = chi2(1, 0:nchans)/dble(nit_included - 1)
    else
      chi2(1, 0:nchans) = 0d0
    end if
    write (*, *) '-------'
    ncalls0 = ncalls*kpoint_iter ! return number of points used
    if (double_points) then
      itmax = 2
    else
      itmax = nit_included
    end if
    do kchan = 1, nchans
      if (regridded(kchan)) then
        ! set equal to number of points used for the last update
        nhits_in_grids(kchan) = sum(nhits(1:nint_used, 1, kchan), dim=1)
      end if
    end do
    call write_grids_to_file
    call write_results
  end subroutine finalise_mint

  subroutine write_results
    implicit none
    integer :: kchan
    write (*, *) 'Final result [ABS]:', ans(1, 0), ' +/-', unc(1, 0)
    write (*, *) 'Final result:', ans(2, 0), ' +/-', unc(2, 0)
    write (*, *) 'chi**2 per D.o.F.:', chi2(1, 0)
    open (unit=58, file='results.dat', status='unknown')
    do kchan = 0, nchans
      write (58, *) ans(1, kchan), unc(2, kchan), 0d0, 0, 0, 0, 0, 0d0, 0d0, ans(2, kchan)
    end do
    close (58)
  end subroutine write_results

  subroutine write_channel_info
    implicit none
    integer :: kchan, np
    do kchan = 1, nchans
      np = sum(nhits(1:nint_used, 1, kchan))
      write (*, 250) 'channel', kchan, ':', iconfigs(kchan) &
        , regridded(kchan), np, nhits_in_grids(kchan) &
        , ans_chan(kchan), ans(2, kchan), virtual_fraction(kchan)
    end do
    return
250 format(a7, i5, 1x, a1, 1x, i5, 1x, l1, 1x, i8, 1x, i8, 2x, e10.4, 2x, e10.4, 2x, e10.4)
  end subroutine write_channel_info

  subroutine setup_imode_m1
    implicit none
    imode = 0
    min_it = min_it0
    ans_chan(1:nchans) = ans(1, 1:nchans)
    ans_chan(0) = sum(ans(1, 1:nchans))
  end subroutine setup_imode_m1

  subroutine setup_imode_0
    implicit none
    min_it = min_it0
    call reset_mint_grids
  end subroutine setup_imode_0

  subroutine reset_mint_grids
    implicit none
    integer :: kint
    do kint = 0, nint_used
      xgrid(kint, 1:ndim, 1:nchans) = dble(kint)/nint_used
    end do
    nhits(1:nint_used, 1:ndim, 1:nchans) = 0
    regridded(1:nchans) = .true.
    nhits_in_grids(1:nchans) = 0
    if (use_poly_virtual) then
      call init_polyfit(born_dimensions, nchans, n_ord_virt, 1000)
    else
      call init_ave_virt
    end if
    virtual_fraction(1:nchans) = max(virt_fraction, min_virt_fraction)
    average_virtual(0:n_ave_virt, 1:nchans) = 0d0
    ans_chan(0:nchans) = 0d0
    if (double_points) then
      ! when doubling points, start with the very first channel only. For the
      ! first iteration, we compute each channel separately.
      ans_chan(0) = 1d0
      ans_chan(1) = 1d0
      ncalls0 = ncalls0/nchans
    end if
  end subroutine reset_mint_grids

  subroutine setup_common
    implicit none
    nit = 0
    nit_included = 0
    ans(1:nintegrals, 0:nchans) = 0d0
    unc(1:nintegrals, 0:nchans) = 0d0
    ans3(1:nintegrals, 1:3) = 0d0
    unc3(1:nintegrals, 1:3) = 0d0
    HwU_values(1:2) = 0d0
  end subroutine setup_common

  subroutine write_grids_to_file
! Write the MINT integration grids to file
    implicit none
    integer :: i, j, k, kchan, poly_virtual_flag
    poly_virtual_flag = merge(1, 0, use_poly_virtual)
    open (unit=12, file='mint_grids', status='replace', action='write')
    write (12, *) mint_grid_magic, mint_grid_format_version
    write (12, *) 'META', ndim, nchans, nintegrals, nintervals, &
         nintervals_virt, born_dimensions, n_ord_virt, poly_virtual_flag
    do kchan = 1, nchans
      do j = 0, nintervals
        write (12, *) 'AVE', (xgrid(j, i, kchan), i=1, ndim)
      end do
      if (.not. use_poly_virtual) then
        do j = 1, nintervals_virt
          do k = 0, n_ord_virt
            write (12, *) 'AVE', &
                 (ave_virt(j, i, k, kchan), i=1, born_dimensions)
          end do
        end do
      end if
      write (12, *) 'SUM', (ans(i, kchan), i=1, nintegrals)
      write (12, *) 'QSM', (unc(i, kchan), i=1, nintegrals)
      write (12, *) 'SPE', ncalls0, itmax, nhits_in_grids(kchan)
      write (12, *) 'AVE', virtual_fraction(kchan), average_virtual(0, kchan)
    end do
    if (use_poly_virtual) call save_polyfit(12)
    close (12)
  end subroutine write_grids_to_file

  subroutine read_grids_from_file
! Read the MINT integration grids from file
    implicit none
    integer :: i, j, k, kchan, idum, ios, version
    integer :: stored_ndim, stored_nchans, stored_nintegrals
    integer :: stored_intervals, stored_virtual_intervals
    integer :: stored_born_dimensions, stored_virtual_orders
    integer :: stored_poly_virtual
    integer, dimension(maxchannels) :: points
    character(len=3) :: dummy
    character(len=16) :: record_name
    double precision, allocatable :: stored_average(:)
    if (.not. mint_state_initialized) call initialize_mint_state
    open (unit=12, file='mint_grids', status='old', action='read')
    read (12, *, iostat=ios) record_name, version
    if (ios /= 0 .or. trim(record_name) /= mint_grid_magic) then
      write (*, *) 'Unsupported legacy MINT grid file; regenerate the grid'
      stop 1
    end if
    if (version /= mint_grid_format_version) then
      write (*, *) 'Unsupported MINT grid format version:', version
      stop 1
    end if
    read (12, *, iostat=ios) record_name, stored_ndim, stored_nchans, &
         stored_nintegrals, stored_intervals, stored_virtual_intervals, &
         stored_born_dimensions, stored_virtual_orders, &
         stored_poly_virtual
    if (ios /= 0 .or. trim(record_name) /= 'META') then
      write (*, *) 'Malformed MINT grid metadata record'
      stop 1
    end if
    if (stored_ndim /= ndim .or. stored_nchans /= nchans .or. &
        stored_nintegrals /= nintegrals .or. &
        stored_intervals /= nintervals .or. &
        stored_virtual_intervals /= nintervals_virt .or. &
        stored_born_dimensions /= born_dimensions) then
      write (*, *) 'MINT grid metadata is incompatible with this process'
      stop 1
    end if
    if (stored_virtual_orders < 0 .or. &
        stored_virtual_orders > n_ave_virt) then
      write (*, *) 'Invalid virtual-order count in MINT grid:', &
           stored_virtual_orders
      stop 1
    end if
    if (stored_poly_virtual /= merge(1, 0, use_poly_virtual)) then
      write (*, *) 'MINT grid virtual-approximation policy is incompatible'
      stop 1
    end if
    ans(1, 0) = 0d0
    unc(1, 0) = 0d0
    ave_virt = 0d0
    average_virtual = 0d0
    if (.not. use_poly_virtual) allocate(stored_average(born_dimensions))
    do kchan = 1, nchans
      do j = 0, nintervals
        read (12, *) dummy, (xgrid(j, i, kchan), i=1, ndim)
      end do
      if (.not. use_poly_virtual) then
        do j = 1, nintervals_virt
          do k = 0, stored_virtual_orders
            read (12, *) dummy, stored_average
            if (k <= n_ord_virt) &
                 ave_virt(j, 1:born_dimensions, k, kchan) = &
                 stored_average
          end do
        end do
      end if
      read (12, *) dummy, (ans(i, kchan), i=1, nintegrals)
      read (12, *) dummy, (unc(i, kchan), i=1, nintegrals)
      read (12, *) dummy, idum, idum, nhits_in_grids(kchan)
      read (12, *) dummy, virtual_fraction(kchan), average_virtual(0, kchan)
      ans(1, 0) = ans(1, 0) + ans(1, kchan)
      unc(1, 0) = unc(1, 0) + unc(1, kchan)**2
    end do
    unc(1, 0) = sqrt(unc(1, 0))
    if (allocated(stored_average)) deallocate(stored_average)
    ! polyfit stuff:
    if (use_poly_virtual) then
      if (stored_virtual_orders == n_ord_virt) then
        do kchan = 1, nchans
          read (12, *) dummy, points(kchan)
        end do
        do kchan = 1, nchans
          backspace (12)
        end do
        call init_polyfit(born_dimensions, nchans, n_ord_virt, &
                          maxval(points(1:nchans)))
        call restore_polyfit(12)
        call do_polyfit()
      else
        ! A cheap no-loop grid still contains a valid phase-space map.  Its
        ! polynomial payload has fewer coefficient columns, so retain the
        ! map and start only the virtual residual model from an empty state.
        call init_polyfit(born_dimensions, nchans, n_ord_virt, 1000)
      end if
    end if
    close (12)
! check for zero cross-section: if restoring grids corresponding to
! sigma=0, just terminate the run
    if (imode .ne. 0 .and. ans(1, 0) .eq. 0d0 .and. unc(1, 0) .eq. 0d0) then
      call initplot()
      call close_run_zero_res
      stop 0
    end if
  end subroutine read_grids_from_file

  subroutine double_grid(kdim, kchan)
    implicit none
    integer :: kchan, kdim, i
    do i = nint_used, 1, -1
      xgrid(i*2, kdim, kchan) = xgrid(i, kdim, kchan)
      xgrid(i*2 - 1, kdim, kchan) = (xgrid(i, kdim, kchan) + xgrid(i - 1, kdim, kchan))/2d0
      if ((.not. regridded(kchan)) .and. (.not. reset)) then
        nhits(i*2, kdim, kchan) = nhits(i, kdim, kchan)/2
        nhits(i*2 - 1, kdim, kchan) = nhits(i, kdim, kchan) - nhits(i*2, kdim, kchan)
        if (nhits(i, kdim, kchan) .ne. 0) then
          xacc(i*2, kdim, kchan) = xacc(i, kdim, kchan)*nhits(i*2, kdim, kchan)/dble(nhits(i, kdim, kchan))
          xacc(i*2 - 1, kdim, kchan) = xacc(i, kdim, kchan)*nhits(i*2 - 1, kdim, kchan)/dble(nhits(i, kdim, kchan))
        else
          xacc(i*2, kdim, kchan) = 0d0
          xacc(i*2 - 1, kdim, kchan) = 0d0
        end if
      end if
    end do
  end subroutine double_grid

  subroutine regrid(kdim, kchan)
    implicit none
    integer :: kdim, kchan, kint, jint
    double precision :: r, total
    double precision, parameter :: tiny = 1d-8
! compute total number of points and update grids if large
    regridded(kchan) = .false.
    if (sum(nhits(1:nint_used, kdim, kchan), dim=1) .lt. nint(0.9*nhits_in_grids(kchan))) return
    regridded(kchan) = .true.
! Use the same smoothing as in VEGAS uses for the grids, i.e. use the
! average of the central and the two neighbouring grid points: (Only do
! this if we are already at the maximum intervals, because the doubling
! of the grids also includes a smoothing).
    if (nint_used .eq. nintervals) then
      call smooth_xacc(kdim, kchan)
    end if
    do kint = 1, nint_used
      if (nhits(kint, kdim, kchan) .ne. 0) then
        xacc(kint, kdim, kchan) = abs(xacc(kint, kdim, kchan))/nhits(kint, kdim, kchan)
      else
        xacc(kint, kdim, kchan) = 0d0
      end if
    end do
! Overwrite xacc so that it accumulates the cross section with each
! successive interval.  It already contains a factor equal to the
! interval size. Thus the integral of rho is performed by summing up
    total = sum(xacc(1:nint_used, kdim, kchan), dim=1)
    do kint = 1, nint_used
      if (nhits(kint, kdim, kchan) .ne. 0) then
!     take logarithm to help convergence (taken from LO dsample.f)
        if (xacc(kint, kdim, kchan) .ne. total) then
          xacc(kint, kdim, kchan) = ((xacc(kint, kdim, kchan)/total - 1d0)/log(xacc(kint, kdim, kchan)/total))**1.5
        else
          xacc(kint, kdim, kchan) = 1d0
        end if
        xacc(kint, kdim, kchan) = xacc(kint - 1, kdim, kchan) + abs(xacc(kint, kdim, kchan))
      else
        xacc(kint, kdim, kchan) = xacc(kint - 1, kdim, kchan)
      end if
    end do
! No valid points. Simply return
    if (xacc(nint_used, kdim, kchan) .eq. 0d0) return
! normalise xacc so that it goes from 0 to 1.
    xacc(1:nint_used, kdim, kchan) = xacc(1:nint_used, kdim, kchan)/xacc(nint_used, kdim, kchan)
! Check that we have a reasonable result and update the accumulated results if need be
    do kint = 1, nint_used
      if (xacc(kint, kdim, kchan) .lt. (xacc(kint - 1, kdim, kchan) + tiny)) then
        xacc(kint, kdim, kchan) = xacc(kint - 1, kdim, kchan) + tiny
      end if
    end do
! it could happen that the change above yielded xacc() values greater than 1: one more update needed
    xacc(nint_used, kdim, kchan) = 1d0
    do kint = 1, nint_used
      if (xacc(nint_used - kint, kdim, kchan) .gt. (xacc(nint_used - kint + 1, kdim, kchan) - tiny)) then
        xacc(nint_used - kint, kdim, kchan) = 1d0 - dble(kint)*tiny
      else
        exit
      end if
    end do
! adjust 'xgrid_new' (temporary grid) so that each element contains identical cross section
    xgrid_new(0, kdim) = 0d0
    do kint = 1, nint_used
      r = dble(kint)/dble(nint_used)
      do jint = 1, nint_used
        if (r .lt. xacc(jint, kdim, kchan)) then
          xgrid_new(kint, kdim) = xgrid(jint - 1, kdim, kchan) + (r - xacc(jint - 1, kdim, kchan))/ &
                                  (xacc(jint, kdim, kchan) - xacc(jint - 1, kdim, kchan))* &
                                  (xgrid(jint, kdim, kchan) - xgrid(jint - 1, kdim, kchan))
          goto 11
        end if
      end do
      if (jint .ne. nint_used + 1 .and. kint .ne. nint_used) then
        write (*, *) 'ERROR', jint, nint_used
        stop 1
      end if
11    continue
    end do
    xgrid_new(nint_used, kdim) = 1d0
  end subroutine regrid

  subroutine smooth_xacc(kdim, kchan)
    implicit none
    integer :: kdim, kchan, kint, kk, itot, kkint
    double precision :: tot
    integer, parameter :: isize = 1
    integer, dimension(1:nintervals) :: local_nhits
    double precision, dimension(1:nintervals) :: local_xacc
    do kint = 1, nint_used
      tot = 0d0
      itot = 0
      do kk = -isize, isize
        kkint = kint + kk
        if (kkint .le. 0) kkint = 1
        if (kkint .ge. nint_used + 1) kkint = nint_used
        tot = tot + xacc(kkint, kdim, kchan)
        itot = itot + nhits(kkint, kdim, kchan)
      end do
      local_xacc(kint) = tot/dble(2*isize + 1)
      local_nhits(kint) = nint(itot/dble(2*isize + 1))
    end do
    xacc(1:nint_used, kdim, kchan) = local_xacc(1:nint_used)
    nhits(1:nint_used, kdim, kchan) = local_nhits(1:nint_used)
  end subroutine smooth_xacc

  subroutine nextlexi(iii, kkk, iret)
! kkk: array of integers 1 <= kkk(j) <= iii(j), j=1,ndim
! at each call iii is increased lexicographycally.
! for example, starting from ndim=3, kkk=(1,1,1), iii=(2,3,2)
! subsequent calls to nextlexi return
!         kkk(1)      kkk(2)      kkk(3)    iret
! 0 calls   1           1           1       0
! 1         1           1           2       0
! 2         1           2           1       0
! 3         1           2           2       0
! 4         1           3           1       0
! 5         1           3           2       0
! 6         2           1           1       0
! 7         2           1           2       0
! 8         2           2           1       0
! 9         2           2           2       0
! 10        2           3           1       0
! 11        2           3           2       0
! 12        2           3           2       1
    implicit none
    integer :: k, iret
    integer, dimension(ndimmax) :: kkk, iii
    k = ndim
1   continue
    if (kkk(k) .lt. iii(k)) then
      kkk(k) = kkk(k) + 1
      iret = 0
      return
    else
      kkk(k) = 1
      k = k - 1
      if (k .eq. 0) then
        iret = 1
        return
      end if
      goto 1
    end if
  end subroutine nextlexi

  subroutine init_ave_virt
    implicit none
    if (n_ord_virt .gt. n_ave_virt) then
      write (*, *) 'Too many grids to keep track off', n_ord_virt, n_ave_virt
      stop 1
    end if
    nvirt(1:nint_used_virt, 1:born_dimensions, &
          0:n_ord_virt, 1:nchans) = 0
    ave_virt(1:nint_used_virt, 1:born_dimensions, &
             0:n_ord_virt, 1:nchans) = 0d0
    nvirt_acc(1:nint_used_virt, 1:born_dimensions, &
              0:n_ord_virt, 1:nchans) = 0
    ave_virt_acc(1:nint_used_virt, 1:born_dimensions, &
                 0:n_ord_virt, 1:nchans) = 0d0
    ave_born_acc(1:nint_used_virt, 1:born_dimensions, &
                 0:n_ord_virt, 1:nchans) = 0d0
  end subroutine init_ave_virt

  subroutine get_ave_virt(x, k_ord_virt)
    implicit none
    integer :: kdim, ncell, k_ord_virt
    double precision, dimension(ndimmax) :: x
    average_virtual(k_ord_virt, ichan) = 0d0
    do kdim = 1, born_dimensions
      ncell = min(int(x(kdim)*nint_used_virt) + 1, nint_used_virt)
      average_virtual(k_ord_virt, ichan) = average_virtual(k_ord_virt, ichan) &
                                           + ave_virt(ncell, kdim, k_ord_virt, ichan)
    end do
    average_virtual(k_ord_virt, ichan) = &
         average_virtual(k_ord_virt, ichan)/born_dimensions
  end subroutine get_ave_virt

  subroutine fill_ave_virt(x, k_ord_virt, virtual, born)
    implicit none
    integer :: kdim, ncell, k_ord_virt
    double precision, dimension(ndimmax) :: x
    double precision :: virtual, born
    do kdim = 1, born_dimensions
      ncell = min(int(x(kdim)*nint_used_virt) + 1, nint_used_virt)
      nvirt_acc(ncell, kdim, k_ord_virt, ichan) = nvirt_acc(ncell, kdim, k_ord_virt, ichan) + 1
      ave_virt_acc(ncell, kdim, k_ord_virt, ichan) = ave_virt_acc(ncell, kdim, k_ord_virt, ichan) + virtual
      ave_born_acc(ncell, kdim, k_ord_virt, ichan) = ave_born_acc(ncell, kdim, k_ord_virt, ichan) + born
    end do
  end subroutine fill_ave_virt

  subroutine regrid_ave_virt(k_ord_virt)
    implicit none
    integer kchan, kdim, i, k_ord_virt
! need to solve for k_new = (virt+k_old*born)/born = virt/born + k_old
    do kchan = 1, nchans
      do kdim = 1, born_dimensions
        do i = 1, nint_used_virt
          if (ave_born_acc(i, kdim, k_ord_virt, kchan) .eq. 0d0) cycle
          if (ave_virt(i, kdim, k_ord_virt, kchan) .eq. 0d0) then ! i.e. first iteration
            ave_virt(i, kdim, k_ord_virt, kchan) = ave_virt_acc(i, kdim, k_ord_virt, kchan) &
                                                   /ave_born_acc(i, kdim, k_ord_virt, kchan) + ave_virt(i, kdim, k_ord_virt, kchan)
          else  ! give some importance to the iterations already done
            ave_virt(i, kdim, k_ord_virt, kchan) = (nvirt_acc(i, kdim, k_ord_virt, kchan)* &
                                                    ave_virt_acc(i, kdim, k_ord_virt, kchan) &
                                                    /ave_born_acc(i, kdim, k_ord_virt, kchan) &
                                                    + nvirt(i, kdim, k_ord_virt, kchan) &
                                                    *ave_virt(i, kdim, k_ord_virt, kchan)) &
                                                   /dble(nvirt_acc(i, kdim, k_ord_virt, kchan) &
                                                         + nvirt(i, kdim, k_ord_virt, kchan))
          end if
        end do
      end do
    end do
! reset the acc values
    nvirt(1:nint_used_virt, 1:born_dimensions, &
          k_ord_virt, 1:nchans) = &
      nvirt(1:nint_used_virt, 1:born_dimensions, &
            k_ord_virt, 1:nchans) &
      + nvirt_acc(1:nint_used_virt, 1:born_dimensions, &
                  k_ord_virt, 1:nchans)
    nvirt_acc(1:nint_used_virt, 1:born_dimensions, &
              k_ord_virt, 1:nchans) = 0
    ave_born_acc(1:nint_used_virt, 1:born_dimensions, &
                 k_ord_virt, 1:nchans) = 0d0
    ave_virt_acc(1:nint_used_virt, 1:born_dimensions, &
                 k_ord_virt, 1:nchans) = 0d0
  end subroutine regrid_ave_virt

  subroutine double_ave_virt(k_ord_virt)
    implicit none
    integer :: kdim, i, k_ord_virt, kchan
    do kchan = 1, nchans
      do kdim = 1, born_dimensions
        do i = nint_used_virt, 1, -1
          ave_virt(i*2, kdim, k_ord_virt, kchan) = ave_virt(i, kdim, k_ord_virt, kchan)
          if (nvirt(i, kdim, k_ord_virt, kchan) .ne. 0) then
            nvirt(i*2, kdim, k_ord_virt, kchan) = max(nvirt(i, kdim, k_ord_virt, kchan)/2, 1)
          else
            nvirt(i*2, kdim, k_ord_virt, kchan) = 0
          end if
          if (i .ne. 1) then
            ave_virt(i*2 - 1, kdim, k_ord_virt, kchan) = (ave_virt(i, kdim, k_ord_virt, kchan) &
                                                          + ave_virt(i - 1, kdim, k_ord_virt, kchan))/2d0
            if (nvirt(i, kdim, k_ord_virt, kchan) + nvirt(i - 1, kdim, k_ord_virt, kchan) .ne. 0) then
              nvirt(i*2 - 1, kdim, k_ord_virt, kchan) = &
                max((nvirt(i, kdim, k_ord_virt, kchan) + nvirt(i - 1, kdim, k_ord_virt, kchan))/4, 1)
            else
              nvirt(i*2 - 1, kdim, k_ord_virt, kchan) = 0
            end if
          else
            if (nvirt(1, kdim, k_ord_virt, kchan) .ne. 0) then
              nvirt(1, kdim, k_ord_virt, kchan) = max(nvirt(1, kdim, k_ord_virt, kchan)/2, 1)
            else
              nvirt(1, kdim, k_ord_virt, kchan) = 0
            end if
          end if
        end do
      end do
    end do
  end subroutine double_ave_virt

  subroutine get_channel
! Picks one random 'ichan' among the 'nchans' integration channels and
! updates the selected integration-channel state.
    implicit none
    double precision :: trgt, total
    double precision, external :: ran2
    if (nchans .eq. 1) then
      ichan = 1
      iconfig = iconfigs(ichan)
      vol_chan = 1d0
    elseif (nchans .gt. 1) then
      if (factorized_mint_uses_uniform_channels()) then
        ! The absolute multiplicative integrand can be dominated by rare
        ! loop--real products.  Probabilities learned from that integral
        ! consequently starve otherwise important production maps and make
        ! the lower-order projections extremely noisy.  Keep sampling the
        ! independently adapted channel grids uniformly; the AMP-based
        ! pointwise partition still performs the multi-channel enhancement.
        ichan = int(ran2()*nchans) + 1
        iconfig = iconfigs(ichan)
        vol_chan = 1d0/dble(nchans)
      elseif (ans_chan(0) .le. 0d0) then
!     pick one at random (flat)
        ichan = int(ran2()*nchans) + 1
        iconfig = iconfigs(ichan)
        vol_chan = 1d0/dble(nchans)
      else
!     pick one at random (weighted by cross section)
        total = sum(ans_chan(1:nchans))
        if (abs(total - ans_chan(0))/(total + ans_chan(0)) .gt. 1d-8) then
          write (*, *) 'ERROR: total should be equal to ans', total, ans_chan(0)
          stop 1
        end if
        trgt = ans_chan(0)*ran2()
        total = 0d0
        ichan = 0
        do while (total .lt. trgt)
          ichan = ichan + 1
          total = total + ans_chan(ichan)
        end do
        if (ichan .eq. 0 .or. ichan .gt. nchans) then
          write (*, *) 'ERROR: ichan cannot be zero or larger than nchans', ichan, nchans
          stop 1
        end if
        iconfig = iconfigs(ichan)
        vol_chan = ans_chan(ichan)/ans_chan(0)
      end if
    end if
  end subroutine get_channel

  subroutine close_run_zero_res
    implicit none
    integer :: kchan
    xgrid(0:nintervals, 1:ndim, 1:nchans) = 0d0
    ave_virt(1:nintervals_virt, 1:born_dimensions, &
             0:n_ord_virt, 1:nchans) = 0d0
    ans(1:nintegrals, 1:nchans) = 0d0
    unc(1:nintegrals, 1:nchans) = 0d0
    nhits_in_grids(1:nchans) = 0
    virtual_fraction(1:nchans) = 1d0
    average_virtual(0:n_ord_virt, 1:nchans) = 0d0
    call write_grids_to_file
    call write_results
    call regrid_MC_integer
    open (unit=12, file='res.dat', status='unknown')
    do kchan = 0, nchans
      write (12, *) ans(1, kchan), unc(1, kchan), ans(2, kchan), unc(2, kchan) &
        , itmax, ncalls0, 0d0
    end do
    close (12)
  end subroutine close_run_zero_res

  subroutine initialize_even_random_numbers
! Recompute the number of calls. Uses the algorithm from VEGAS
    implicit none
! Make sure that hypercubes are newly initialized
    firsttime = .true.
! Number of elements in which we can split one dimension
    ng = int((dble(ncalls0)/2d0)**(1d0/dble(ndim)))
! Total number of hypercubes
    k = ng**ndim
! Number of PS points in each hypercube (at least 2)
    npg = max(ncalls0/k, 2)
! Number of PS points for this iteration
    ncalls = npg*k
  end subroutine initialize_even_random_numbers

  function get_ran()
    implicit none
    double precision :: get_ran
    double precision, external ::  ran2
    integer :: i, iret
    if (firsttime) then
! initialise the hypercubes
      even_dng = 1d0/dble(ng)
      even_current_dim = 0
      do i = 1, ndim
        even_iii(i) = ng
        even_kkk(i) = 1
      end do
      firsttime = .false.
    end if
    even_current_dim = mod(even_current_dim, ndim) + 1
! This is the random number in the hypercube 'k' for current_dim
    get_ran = even_dng*(ran2() + dble(even_kkk(even_current_dim) - 1))
! Got random numbers for all dimensions, update kkk() for the next call
    if (even_current_dim .eq. ndim) then
      call nextlexi(even_iii, even_kkk, iret)
      if (iret .eq. 1) then
        call nextlexi(even_iii, even_kkk, iret)
      end if
    end if
  end function get_ran

end module mint_module
