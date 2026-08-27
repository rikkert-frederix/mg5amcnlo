module test_soft_col_limits_module
  use decay_chain_metadata, only: real_phase_space_dimension, &
                                  has_decay_chains
  use decay_chain_kinematics, only: fks_leg_mass, &
                                    minimum_core_final_mass
  use process_dimensions, only: nexternal, nincoming, fks_configs, &
                                nsplitorders, amp_split_size, order_names, &
                                validate_process_and_born_dimensions
  use fks_metadata, only: fks_i_d, pdg_type_d, validate_fks_metadata
  use run_state, only: lpp, ebeam, ptj, ptl, &
                       mll, mll_sf, ptgmin
  use fnlo_process_common, only: is_a_j, is_a_lp, is_a_lm, is_a_ph, &
                                 soft_counterevent, collinear_counterevent, &
                                 real_event
  use split_orders, only: amp_split_pos_to_orders
  use genps_fks, only: generate_momenta
  use fks_diagnostics, only: xprintout, checkres2
  use mint_module, only: ndim, iconfig, ichan, iconfigs, new_point
  use fks_random_module, only: random_unit_interval
  use fks_singular_module, only: fill_configurations_common, setfksfactor
  implicit none
  private

  integer, parameter :: maximum_limit_points = 15
  integer, parameter :: used_limit_points = 10
  integer, parameter :: random_vector_size = 99
  double precision, parameter :: zero = 0d0
  double precision, parameter :: one = 1d0
  double precision, parameter :: max_fail = 0.3d0

  integer, allocatable, save :: born_mapconfig(:)
  double precision, allocatable, save :: generated_masses(:)
  logical, save :: generated_data_initialized = .false.

  double precision, allocatable, save :: fxl(:), wfxl(:)
  double precision, allocatable, save :: limit(:), wlimit(:)
  double precision, allocatable, save :: lxp(:, :), xp(:, :, :)
  double precision, allocatable, save :: p(:, :), x(:)
  double precision, allocatable, save :: fxl_split(:, :)
  double precision, allocatable, save :: wfxl_split(:, :)
  double precision, allocatable, save :: limit_split(:, :)
  double precision, allocatable, save :: wlimit_split(:, :)
  integer, allocatable, save :: orders(:), nerr(:)
  double precision, allocatable, save :: fail_frac(:)

  double precision, allocatable, save :: event_momenta(:, :, :)
  double precision, allocatable, save :: event_jacobian(:)
  double precision, allocatable, save :: p_born(:, :)
  double precision, allocatable, save :: event_fks_momentum(:, :)
  double precision, allocatable, save :: event_xi(:), event_y(:)
  double precision, allocatable, save :: etmin(:), etmax(:)
  double precision, allocatable, save :: mxxmin(:, :)
  double precision, allocatable, save :: amp_split(:)
  logical, save :: work_state_initialized = .false.

  public :: run_test_soft_col_limits
  public :: init_test_limits_data
  public :: finalize_test_soft_col_limits

  interface
    subroutine init_test_limits_data_bridge()
      implicit none
    end subroutine init_test_limits_data_bridge

    subroutine test_limits_select_fks_bridge(configuration, i_fks, j_fks)
      implicit none
      integer, intent(in) :: configuration
      integer, intent(out) :: i_fks, j_fks
    end subroutine test_limits_select_fks_bridge

    subroutine test_limits_set_nndim_bridge(nndim_value)
      implicit none
      integer, intent(in) :: nndim_value
    end subroutine test_limits_set_nndim_bridge

    subroutine test_limits_set_controls_bridge(xi_fixed, y_fixed, &
                                               calculated_born, soft_test, collinear_test)
      implicit none
      double precision, intent(in) :: xi_fixed, y_fixed
      logical, intent(in) :: calculated_born, soft_test, collinear_test
    end subroutine test_limits_set_controls_bridge

    subroutine test_limits_sync_state_bridge(event_momenta_out, &
                                             event_jacobian_out, born_momenta, &
                                             event_xi_out, event_y_out, &
                                             event_fks_momentum_out)
      use process_dimensions, only: nexternal
      use fnlo_process_common, only: real_event
      implicit none
      double precision, intent(out) :: &
        event_momenta_out(0:3, nexternal, 0:real_event)
      double precision, intent(out) :: event_jacobian_out(0:real_event)
      double precision, intent(out) :: born_momenta(0:3, nexternal - 1)
      double precision, intent(out) :: event_xi_out(0:real_event)
      double precision, intent(out) :: event_y_out(0:real_event)
      double precision, intent(out) :: &
        event_fks_momentum_out(0:3, 0:real_event)
    end subroutine test_limits_sync_state_bridge

    subroutine test_limits_setcuts_bridge(etmin_out, etmax_out, mxxmin_out)
      use process_dimensions, only: nexternal, nincoming
      implicit none
      double precision, intent(out) :: &
        etmin_out(nincoming + 1:nexternal - 1)
      double precision, intent(out) :: &
        etmax_out(nincoming + 1:nexternal - 1)
      double precision, intent(out) :: &
        mxxmin_out(nincoming + 1:nexternal - 1, &
                   nincoming + 1:nexternal - 1)
    end subroutine test_limits_setcuts_bridge

    subroutine test_limits_sreal_bridge(event_slot, momentum, xi_fks, y_fks, &
                                        weight, split_weights)
      use process_dimensions, only: nexternal, amp_split_size
      implicit none
      integer, intent(in) :: event_slot
      double precision, intent(in) :: momentum(0:3, nexternal)
      double precision, intent(in) :: xi_fks, y_fks
      double precision, intent(out) :: weight
      double precision, intent(out) :: split_weights(amp_split_size)
    end subroutine test_limits_sreal_bridge

    subroutine setrun()
      implicit none
    end subroutine setrun

    subroutine setpara(card_name)
      implicit none
      character(len=*), intent(in) :: card_name
    end subroutine setpara

    subroutine init_fks_singular_bridge()
      implicit none
    end subroutine init_fks_singular_bridge

    subroutine sborn(momentum, weight)
      use process_dimensions, only: nexternal
      implicit none
      double precision, intent(in) :: momentum(0:3, nexternal - 1)
      complex(kind=kind(0d0)), intent(out) :: weight(2)
    end subroutine sborn

  end interface

contains

  subroutine run_test_soft_col_limits()
    implicit none
    integer :: i, j, k, l, jj, bs_min, bs_max, iconfig_in, nsofttests
    integer :: ncolltests, imax, iflag, iret, ntry, fks_conf_number
    integer :: fks_loop_min, fks_loop_max, fks_loop, iamp
    integer :: i_fks, j_fks, nfksprocess
    double precision :: wgt, fx, totmass, sister_mass
    double precision :: xi_i_fks_fix_save, y_ij_fks_fix_save
    double precision :: xi_i_fks_fix, y_ij_fks_fix
    logical :: calculated_born, softtest, colltest
    complex(kind=kind(0d0)) :: wgt1(2)

    call validate_process_and_born_dimensions()
    call validate_fks_metadata()
    call initialize_work_state()

    if (fks_configs .eq. 1) then
      if (pdg_type_d(1, fks_i_d(1)) .eq. -21) then
        write (*, *) 'Process generated with [LOonly=QCD]. No tests to do.'
        return
      end if
    end if
    write (*, *) 'Enter xi_i, y_ij to be used in coll/soft tests'
    write (*, *) ' Enter -2 to generate them randomly'
    read (*, *) xi_i_fks_fix_save, y_ij_fks_fix_save

    write (*, *) 'Enter number of tests for soft and collinear limits'
    read (*, *) nsofttests, ncolltests

    call setrun               !Sets up run parameters
    call setpara('param_card.dat') !Sets up couplings and masses
    call init_fks_singular_bridge()
    call fill_configurations_common
    call init_test_limits_data_bridge()
    call validate_generated_data()

!
    write (*, *) 'Give FKS configuration number ("0" loops over all)'
    read (*, *) fks_conf_number

    if (fks_conf_number .eq. 0) then
      fks_loop_min = 1
      fks_loop_max = fks_configs
    else
      fks_loop_min = fks_conf_number
      fks_loop_max = fks_conf_number
    end if

    fks_configuration_loop: do fks_loop = fks_loop_min, fks_loop_max
      nFKSprocess = fks_loop
      call test_limits_select_fks_bridge(nfksprocess, i_fks, j_fks)
      write (*, *) ''
      write (*, *) '================================================='
      write (*, *) ''
      write (*, *) 'NEW FKS CONFIGURATION:'
      write (*, *) 'FKS configuration number is ', nFKSprocess
      write (*, *) 'FKS partons are: i=', i_fks, '  j=', j_fks
      write (*, *) 'with PDGs:       i=', pdg_type_d(nfksprocess, i_fks), &
        '  j=', pdg_type_d(nfksprocess, j_fks)
!
      ndim = real_phase_space_dimension()
      if (abs(lpp(1)) .ge. 1) ndim = ndim + 1
      if (abs(lpp(2)) .ge. 1) ndim = ndim + 1
      if (ndim < 1 .or. ndim > random_vector_size) then
        call fail_test_limits('invalid integration dimension')
      end if
      call test_limits_set_nndim_bridge(ndim)

      call test_limits_setcuts_bridge(etmin, etmax, mxxmin)
! When doing hadron-hadron collision reduce the effect collision energy.
! Note that tests are always performed at fixed energy with Bjorken x=1.
      totmass = 0.0d0
      do i = nincoming + 1, nexternal
        if (is_a_j(i) .and. i .ne. nexternal) then
          totmass = totmass + max(ptj, generated_masses(i))
        elseif ((is_a_lp(i) .or. is_a_lm(i)) .and. i .ne. nexternal) then
          totmass = totmass + max(mll/2d0, mll_sf/2d0, ptl, &
                                  generated_masses(i))
        elseif (is_a_ph(i)) then
          totmass = totmass + ptgmin
        else
          if (any(mxxmin(i, i + 1:nexternal - 1) .gt. 0d0)) then
            do k = i + 1, nexternal - 1
              if (mxxmin(i, k) .gt. 0d0) then
                totmass = totmass + mxxmin(i, k)
              end if
            end do
          elseif (etmin(i) .gt. 0d0) then
            totmass = totmass + max(etmin(i), generated_masses(i))
          else
            totmass = totmass + generated_masses(i)
          end if
        end if
      end do
      if (has_decay_chains()) then
        totmass = max(totmass, minimum_core_final_mass())
      end if
      if (lpp(1) .ne. 0) ebeam(1) = max(ebeam(1) &
        &      /20d0, totmass)
      if (lpp(2) .ne. 0) ebeam(2) = max(ebeam(2) &
        &      /20d0, totmass)

      write (*, *) '  '
      write (*, *) '  '
      write (*, *) "Enter graph number (iconfig), " &
        &      //"'0' loops over all graphs"
      read (*, *) iconfig_in

      if (iconfig_in .eq. 0) then
        bs_min = 1
        bs_max = born_mapconfig(0)
      elseif (iconfig_in .eq. -1) then
        bs_min = 1
        bs_max = 1
      else
        bs_min = iconfig_in
        bs_max = iconfig_in
      end if

      born_configuration_loop: do iconfig = bs_min, bs_max
        ichan = 1
        iconfigs(1) = iconfig
        call setfksfactor
        call test_limits_setcuts_bridge(etmin, etmax, mxxmin)
        ntry = 1

        softtest = .false.
        colltest = .false.
        calculated_born = .false.
        xi_i_fks_fix = 0d0
        y_ij_fks_fix = 0d0
        call test_limits_set_controls_bridge(xi_i_fks_fix, y_ij_fks_fix, &
                                             calculated_born, softtest, colltest)

        do jj = 1, ndim
          x(jj) = random_unit_interval(iconfig)
        end do
        new_point = .true.
        wgt = 1d0
        call generate_momenta(ndim, iconfig, wgt, x, p)
        call sync_momentum_state()
        calculated_born = .false.
        call test_limits_set_controls_bridge(xi_i_fks_fix, y_ij_fks_fix, &
                                             calculated_born, softtest, colltest)
        do while ((wgt .lt. 0 .or. p(0, 1) .le. 0d0 .or. p_born(0, 1) .le. 0d0 &
       &         ) .and. ntry .lt. 1000)
          do jj = 1, ndim
            x(jj) = random_unit_interval(iconfig)
          end do
          new_point = .true.
          wgt = 1d0
          call generate_momenta(ndim, iconfig, wgt, x, p)
          call sync_momentum_state()
          calculated_born = .false.
          call test_limits_set_controls_bridge(xi_i_fks_fix, y_ij_fks_fix, &
                                               calculated_born, softtest, colltest)
          ntry = ntry + 1
        end do
        if (ntry .ge. 1000) then
          write (*, *) 'No points passed cuts...'
          write (12, *) 'ERROR: no points passed cuts... Cannot perform ', &
            'ME tests properly for config', iconfig
          cycle born_configuration_loop
        end if
        call sborn(p_born, wgt1)

        write (*, *) ''
        write (*, *) ''
        write (*, *) ''

        softtest = .true.
        colltest = .false.
        calculated_born = .false.
        xi_i_fks_fix = zero
        y_ij_fks_fix = zero
        call test_limits_set_controls_bridge(xi_i_fks_fix, y_ij_fks_fix, &
                                             calculated_born, softtest, colltest)
        nerr(:) = 0
        imax = used_limit_points
        do j = 1, nsofttests
          do iamp = 1, amp_split_size
            do i = 1, imax
              fxl_split(i, iamp) = 0d0
              wfxl_split(i, iamp) = 0d0
              limit_split(i, iamp) = 0d0
              wlimit_split(i, iamp) = 0d0
            end do
          end do
          if (nsofttests .le. 10) then
            write (*, *) ' '
            write (*, *) ' '
          end if
          y_ij_fks_fix = y_ij_fks_fix_save
          xi_i_fks_fix = 0.1d0
          calculated_born = .false.
          call test_limits_set_controls_bridge(xi_i_fks_fix, y_ij_fks_fix, &
                                               calculated_born, softtest, colltest)
          ntry = 1
          wgt = 1d0
          do jj = 1, ndim
            x(jj) = random_unit_interval(iconfig)
          end do
          new_point = .true.
          call generate_momenta(ndim, iconfig, wgt, x, p)
          call sync_momentum_state()
          do while ((wgt .lt. 0 .or. p(0, 1) .le. 0d0) .and. ntry .lt. 1000)
            wgt = 1d0
            do jj = 1, ndim
              x(jj) = random_unit_interval(iconfig)
            end do
            new_point = .true.
            call generate_momenta(ndim, iconfig, wgt, x, p)
            call sync_momentum_state()
            ntry = ntry + 1
          end do
          if (nsofttests .le. 10) write (*, *) 'ntry', ntry
          calculated_born = .false.
          call test_limits_set_controls_bridge(xi_i_fks_fix, y_ij_fks_fix, &
                                               calculated_born, softtest, colltest)
          call test_limits_sreal_bridge( &
            soft_counterevent, event_momenta(:, :, soft_counterevent), zero, &
            event_y(real_event), fx, amp_split)
          fxl(1) = fx*wgt
          wfxl(1) = event_jacobian(soft_counterevent)
          do iamp = 1, amp_split_size
            fxl_split(1, iamp) = &
              amp_split(iamp)*event_jacobian(soft_counterevent)
            wfxl_split(1, iamp) = event_jacobian(soft_counterevent)
          end do
          call test_limits_sreal_bridge( &
            real_event, p, event_xi(real_event), event_y(real_event), &
            fx, amp_split)
          limit(1) = fx*wgt
          wlimit(1) = wgt
          do iamp = 1, amp_split_size
            limit_split(1, iamp) = amp_split(iamp)*wgt
            wlimit_split(1, iamp) = wgt
          end do

          do k = 1, nexternal
            do l = 0, 3
              lxp(l, k) = event_momenta(l, k, soft_counterevent)
              xp(1, l, k) = p(l, k)
            end do
          end do
          do l = 0, 3
            lxp(l, nexternal + 1) = &
              event_fks_momentum(l, soft_counterevent)
            xp(1, l, nexternal + 1) = &
              event_fks_momentum(l, real_event)
          end do

          do i = 2, imax
            xi_i_fks_fix = xi_i_fks_fix/10d0
            calculated_born = .false.
            call test_limits_set_controls_bridge(xi_i_fks_fix, &
              y_ij_fks_fix, calculated_born, softtest, colltest)
            wgt = 1d0
            call generate_momenta(ndim, iconfig, wgt, x, p)
            call sync_momentum_state()
            calculated_born = .false.
            call test_limits_set_controls_bridge(xi_i_fks_fix, &
              y_ij_fks_fix, calculated_born, softtest, colltest)
            call test_limits_sreal_bridge( &
              soft_counterevent, event_momenta(:, :, soft_counterevent), zero, &
              event_y(real_event), fx, amp_split)
            fxl(i) = fx*wgt
            wfxl(i) = event_jacobian(soft_counterevent)
            do iamp = 1, amp_split_size
              fxl_split(i, iamp) = &
                amp_split(iamp)*event_jacobian(soft_counterevent)
              wfxl_split(i, iamp) = event_jacobian(soft_counterevent)
            end do
            calculated_born = .false.
            call test_limits_set_controls_bridge(xi_i_fks_fix, &
              y_ij_fks_fix, calculated_born, softtest, colltest)
            call test_limits_sreal_bridge( &
              real_event, p, event_xi(real_event), event_y(real_event), &
              fx, amp_split)
            limit(i) = fx*wgt
            wlimit(i) = wgt
            do iamp = 1, amp_split_size
              limit_split(i, iamp) = amp_split(iamp)*wgt
              wlimit_split(i, iamp) = wgt
            end do
            do k = 1, nexternal
              do l = 0, 3
                xp(i, l, k) = p(l, k)
              end do
            end do
            do l = 0, 3
              xp(i, l, nexternal + 1) = &
                event_fks_momentum(l, real_event)
            end do
          end do

          if (nsofttests .le. 10) then
            write (*, *) 'Soft limit:'
            do i = 1, imax
              call xprintout(6, limit(i), fxl(i))
            end do
            do iamp = 1, amp_split_size
              if (limit_split(1, iamp) .ne. 0d0 .or. fxl_split(1 &
    &                  , iamp) .ne. 0d0) then
                write (*, *) '   Split-order', iamp
                call amp_split_pos_to_orders(iamp, orders)
                do i = 1, nsplitorders
                  write (*, *) '      ', trim(order_names(i)), ':', orders(i)
                end do
                do i = 1, imax
                  call xprintout(6, limit_split(i, iamp), fxl_split(i &
  &                        , iamp))
                end do
                iflag = 0
                call checkres2(limit_split(1, iamp), fxl_split(1 &
   &                     , iamp), wlimit_split(1, iamp), wfxl_split(1, iamp) &
   &                     , xp, lxp, iflag, imax, j, i_fks, j_fks &
   &                     , iret)
                write (*, *) 'RETURN CODE', iret
              end if
            end do
!
            write (80, *) '  '
            write (80, *) '****************************'
            write (80, *) '  '
            do k = 1, nexternal + 1
              write (80, *) ''
              write (80, *) 'part:', k
              do l = 0, 3
                write (80, *) 'comp:', l
                do i = 1, 10
                  call xprintout(80, xp(i, l, k), lxp(l, k))
                end do
              end do
            end do
          else
            iflag = 0
            call checkres2(limit, fxl, wlimit, wfxl, xp, lxp, &
     &               iflag, imax, j, i_fks, j_fks, iret)
            nerr(0) = nerr(0) + iret
            ! check the contributions coming from each splitorders
            ! only look at the non vanishing ones
            do iamp = 1, amp_split_size
              if (limit_split(1, iamp) .ne. 0d0 .or. fxl_split(1 &
    &                  , iamp) .ne. 0d0) then
                call checkres2(limit_split(1, iamp), fxl_split(1 &
   &                     , iamp), wlimit_split(1, iamp), wfxl_split(1, iamp) &
   &                     , xp, lxp, iflag, imax, j, i_fks, j_fks &
   &                     , iret)
                nerr(iamp) = nerr(iamp) + iret
              end if
            end do
          end if
        end do
        if (nsofttests .gt. 10) then
          write (*, *) 'Soft tests done for (Born) config', iconfig
          write (*, *) 'Failures:', nerr
          do iamp = 0, amp_split_size
            if (iamp == 1) cycle
            fail_frac(iamp) = nerr(iamp)/dble(nsofttests)
            if (iamp .ne. 0) then
              write (*, fmt="(a,i3,a)", advance="no") 'Split-order', iamp, ': '
              call amp_split_pos_to_orders(iamp, orders)
              do i = 1, nsplitorders
                write (*, fmt="(a,a,i3,a)", advance="no") &
                  trim(order_names(i)), ':', orders(i), '; '
              end do
            else
              write (*, fmt="(a)", advance="no") 'Sum of all orders: '
            end if
            if (fail_frac(iamp) .lt. max_fail) then
              write (*, 401) nFKSprocess, fail_frac(iamp)
            else
              write (*, 402) nFKSprocess, fail_frac(iamp)
            end if
          end do
        end if

        write (*, *) ''
        write (*, *) ''
        write (*, *) ''

        if (has_decay_chains()) then
          sister_mass = fks_leg_mass(nfksprocess, j_fks)
        else
          sister_mass = generated_masses(j_fks)
        end if
        if (sister_mass .ne. 0d0) then
          write (*, *) 'No collinear test for massive j_fks'
          cycle born_configuration_loop
        end if

        softtest = .false.
        colltest = .true.
        calculated_born = .false.
        xi_i_fks_fix = zero
        y_ij_fks_fix = zero
        call test_limits_set_controls_bridge(xi_i_fks_fix, y_ij_fks_fix, &
                                             calculated_born, softtest, colltest)

        nerr(:) = 0
        imax = used_limit_points
        do j = 1, ncolltests
          do iamp = 1, amp_split_size
            do i = 1, imax
              fxl_split(i, iamp) = 0d0
              wfxl_split(i, iamp) = 0d0
              limit_split(i, iamp) = 0d0
              wlimit_split(i, iamp) = 0d0
            end do
          end do
          if (ncolltests .le. 10) then
            write (*, *) ' '
            write (*, *) ' '
          end if

          y_ij_fks_fix = 0.9d0
          xi_i_fks_fix = xi_i_fks_fix_save
          calculated_born = .false.
          call test_limits_set_controls_bridge(xi_i_fks_fix, y_ij_fks_fix, &
                                               calculated_born, softtest, colltest)
          ntry = 1
          wgt = 1d0
          do jj = 1, ndim
            x(jj) = random_unit_interval(iconfig)
          end do
          new_point = .true.
          call generate_momenta(ndim, iconfig, wgt, x, p)
          call sync_momentum_state()
          do while ((wgt .lt. 0 .or. p(0, 1) .le. 0d0) .and. ntry .lt. 1000)
            wgt = 1d0
            do jj = 1, ndim
              x(jj) = random_unit_interval(iconfig)
            end do
            new_point = .true.
            call generate_momenta(ndim, iconfig, wgt, x, p)
            call sync_momentum_state()
            ntry = ntry + 1
          end do
          if (ncolltests .le. 10) write (*, *) 'ntry', ntry
          calculated_born = .false.
          call test_limits_set_controls_bridge(xi_i_fks_fix, y_ij_fks_fix, &
                                               calculated_born, softtest, colltest)
          call test_limits_sreal_bridge( &
            collinear_counterevent, &
            event_momenta(:, :, collinear_counterevent), &
            event_xi(collinear_counterevent), one, fx, amp_split)
          fxl(1) = fx*event_jacobian(collinear_counterevent)
          wfxl(1) = event_jacobian(collinear_counterevent)
          do iamp = 1, amp_split_size
            fxl_split(1, iamp) = &
              amp_split(iamp)*event_jacobian(collinear_counterevent)
            wfxl_split(1, iamp) = event_jacobian(collinear_counterevent)
          end do

          call test_limits_sreal_bridge( &
            real_event, p, event_xi(real_event), event_y(real_event), &
            fx, amp_split)
          limit(1) = fx*wgt
          wlimit(1) = wgt
          do iamp = 1, amp_split_size
            limit_split(1, iamp) = amp_split(iamp)*wgt
            wlimit_split(1, iamp) = wgt
          end do

          do k = 1, nexternal
            do l = 0, 3
              lxp(l, k) = event_momenta(l, k, collinear_counterevent)
              xp(1, l, k) = p(l, k)
            end do
          end do
          do l = 0, 3
            lxp(l, nexternal + 1) = &
              event_fks_momentum(l, collinear_counterevent)
            xp(1, l, nexternal + 1) = &
              event_fks_momentum(l, real_event)
          end do

          do i = 2, imax
            y_ij_fks_fix = 1 - 0.1d0**i
            calculated_born = .false.
            call test_limits_set_controls_bridge(xi_i_fks_fix, &
              y_ij_fks_fix, calculated_born, softtest, colltest)
            wgt = 1d0
            call generate_momenta(ndim, iconfig, wgt, x, p)
            call sync_momentum_state()
            calculated_born = .false.
            call test_limits_set_controls_bridge(xi_i_fks_fix, &
              y_ij_fks_fix, calculated_born, softtest, colltest)
            call test_limits_sreal_bridge( &
              collinear_counterevent, &
              event_momenta(:, :, collinear_counterevent), &
              event_xi(collinear_counterevent), one, fx, amp_split)
            fxl(i) = fx*event_jacobian(collinear_counterevent)
            wfxl(i) = event_jacobian(collinear_counterevent)
            do iamp = 1, amp_split_size
              fxl_split(i, iamp) = &
                amp_split(iamp)*event_jacobian(collinear_counterevent)
              wfxl_split(i, iamp) = event_jacobian(collinear_counterevent)
            end do
            calculated_born = .false.
            call test_limits_set_controls_bridge(xi_i_fks_fix, &
              y_ij_fks_fix, calculated_born, softtest, colltest)
            call test_limits_sreal_bridge( &
              real_event, p, event_xi(real_event), event_y(real_event), &
              fx, amp_split)
            limit(i) = fx*wgt
            wlimit(i) = wgt
            do iamp = 1, amp_split_size
              limit_split(i, iamp) = amp_split(iamp)*wgt
              wlimit_split(i, iamp) = wgt
            end do
            do k = 1, nexternal
              do l = 0, 3
                xp(i, l, k) = p(l, k)
              end do
            end do
            do l = 0, 3
              xp(i, l, nexternal + 1) = &
                event_fks_momentum(l, real_event)
            end do
          end do
          if (ncolltests .le. 10) then
            write (*, *) 'Collinear limit:'
            do i = 1, imax
              call xprintout(6, limit(i), fxl(i))
            end do
            do iamp = 1, amp_split_size
              if (limit_split(1, iamp) .ne. 0d0 .or. fxl_split(1 &
    &                  , iamp) .ne. 0d0) then
                write (*, *) '   Split-order', iamp
                call amp_split_pos_to_orders(iamp, orders)
                do i = 1, nsplitorders
                  write (*, *) '      ', trim(order_names(i)), ':', orders(i)
                end do
                do i = 1, imax
                  call xprintout(6, limit_split(i, iamp), fxl_split(i &
  &                        , iamp))
                end do
                iflag = 1
                call checkres2(limit_split(1, iamp), fxl_split(1 &
   &                     , iamp), wlimit_split(1, iamp), wfxl_split(1, iamp) &
   &                     , xp, lxp, iflag, imax, j, i_fks, j_fks &
   &                     , iret)
                write (*, *) 'RETURN CODE', iret
              end if
            end do
!
            write (80, *) '  '
            write (80, *) '****************************'
            write (80, *) '  '
            do k = 1, nexternal + 1
              write (80, *) ''
              write (80, *) 'part:', k
              do l = 0, 3
                write (80, *) 'comp:', l
                do i = 1, 10
                  call xprintout(80, xp(i, l, k), lxp(l, k))
                end do
              end do
            end do
          else
            iflag = 1
            call checkres2(limit, fxl, wlimit, wfxl, xp, lxp, &
     &               iflag, imax, j, i_fks, j_fks, iret)
            nerr(0) = nerr(0) + iret
            ! check the contributions coming from each splitorders
            ! only look at the non vanishing ones
            do iamp = 1, amp_split_size
              if (limit_split(1, iamp) .ne. 0d0 .or. fxl_split(1 &
    &                  , iamp) .ne. 0d0) then
                call checkres2(limit_split(1, iamp), fxl_split(1, iamp), &
   &                     wlimit_split(1, iamp), wfxl_split(1, iamp), xp, lxp, &
   &                     iflag, imax, j, i_fks, j_fks, iret)
                nerr(iamp) = nerr(iamp) + iret
              end if
            end do
          end if
        end do
        if (ncolltests .gt. 10) then
          write (*, *) 'Collinear tests done for (Born) config', iconfig
          write (*, *) 'Failures:', nerr
          do iamp = 0, amp_split_size
            if (iamp == 1) cycle
            fail_frac(iamp) = nerr(iamp)/dble(nsofttests)
            if (iamp .ne. 0) then
              write (*, fmt="(a,i3,a)", advance="no") 'Split-order', iamp, ': '
              call amp_split_pos_to_orders(iamp, orders)
              do i = 1, nsplitorders
                write (*, fmt="(a,a,i3,a)", advance="no") &
                  trim(order_names(i)), ':', orders(i), '; '
              end do
            else
              write (*, fmt="(a)", advance="no") 'Sum of all orders: '
            end if
            if (fail_frac(iamp) .lt. max_fail) then
              write (*, 501) nFKSprocess, fail_frac(iamp)
            else
              write (*, 502) nFKSprocess, fail_frac(iamp)
            end if
          end do
        end if

      end do born_configuration_loop
    end do fks_configuration_loop

    return
401 format('     Soft test ', i2, ' PASSED. Fraction of failures: ', f4.2)
402 format('     Soft test ', i2, ' FAILED. Fraction of failures: ', f4.2)
501 format('Collinear test ', i2, ' PASSED. Fraction of failures: ', f4.2)
502 format('Collinear test ', i2, ' FAILED. Fraction of failures: ', f4.2)
  end subroutine run_test_soft_col_limits

  subroutine init_test_limits_data(mapconfig_input, &
                                   masses_input)
    implicit none
    integer, intent(in) :: mapconfig_input(0:)
    double precision, intent(in) :: masses_input(:)
    integer :: allocation_status

    if (size(masses_input) /= nexternal) then
      call fail_test_limits('generated mass table has the wrong size')
    end if
    if (ubound(mapconfig_input, 1) < 0) then
      call fail_test_limits('generated Born configuration table is empty')
    end if
    if (mapconfig_input(0) < 0 .or. &
        mapconfig_input(0) > ubound(mapconfig_input, 1)) then
      call fail_test_limits('invalid generated Born configuration count')
    end if

    if (allocated(born_mapconfig)) deallocate (born_mapconfig)
    if (allocated(generated_masses)) deallocate (generated_masses)
    generated_data_initialized = .false.

    allocate (born_mapconfig(0:ubound(mapconfig_input, 1)), &
              generated_masses(nexternal), stat=allocation_status)
    if (allocation_status /= 0) then
      call fail_test_limits('could not allocate generated process data')
    end if
    born_mapconfig = mapconfig_input
    generated_masses = masses_input
    generated_data_initialized = .true.
  end subroutine init_test_limits_data

  subroutine initialize_work_state()
    implicit none
    integer :: allocation_status

    if (work_state_initialized) then
      call validate_work_state()
      return
    end if

    allocate (fxl(maximum_limit_points), wfxl(maximum_limit_points), &
              limit(maximum_limit_points), wlimit(maximum_limit_points), &
              lxp(0:3, nexternal + 1), &
              xp(maximum_limit_points, 0:3, nexternal + 1), &
              p(0:3, nexternal), x(random_vector_size), &
              fxl_split(maximum_limit_points, amp_split_size), &
              wfxl_split(maximum_limit_points, amp_split_size), &
              limit_split(maximum_limit_points, amp_split_size), &
              wlimit_split(maximum_limit_points, amp_split_size), &
              orders(nsplitorders), nerr(0:amp_split_size), &
              fail_frac(0:amp_split_size), &
              event_momenta(0:3, nexternal, 0:real_event), &
              event_jacobian(0:real_event), &
              p_born(0:3, nexternal - 1), &
              event_fks_momentum(0:3, 0:real_event), &
              event_xi(0:real_event), event_y(0:real_event), &
              etmin(nincoming + 1:nexternal - 1), &
              etmax(nincoming + 1:nexternal - 1), &
              mxxmin(nincoming + 1:nexternal - 1, nincoming + 1:nexternal - 1), &
              amp_split(amp_split_size), stat=allocation_status)
    if (allocation_status /= 0) then
      call fail_test_limits('could not allocate limit-test work state')
    end if

    fxl = zero
    wfxl = zero
    limit = zero
    wlimit = zero
    lxp = zero
    xp = zero
    p = zero
    x = zero
    fxl_split = zero
    wfxl_split = zero
    limit_split = zero
    wlimit_split = zero
    orders = 0
    nerr = 0
    fail_frac = zero
    event_momenta = zero
    event_jacobian = zero
    p_born = zero
    event_fks_momentum = zero
    event_xi = zero
    event_y = zero
    etmin = zero
    etmax = zero
    mxxmin = zero
    amp_split = zero
    work_state_initialized = .true.
  end subroutine initialize_work_state

  subroutine sync_momentum_state()
    implicit none

    call validate_work_state()
    call test_limits_sync_state_bridge( &
      event_momenta, event_jacobian, p_born, event_xi, event_y, &
      event_fks_momentum)
  end subroutine sync_momentum_state

  subroutine validate_generated_data()
    implicit none

    if (.not. generated_data_initialized) then
      call fail_test_limits('generated process data were not initialized')
    end if
    if (.not. allocated(born_mapconfig) .or. &
        .not. allocated(generated_masses)) then
      call fail_test_limits('generated process data are incomplete')
    end if
    if (lbound(born_mapconfig, 1) /= 0 .or. &
        size(generated_masses) /= nexternal) then
      call fail_test_limits('generated process data have invalid bounds')
    end if
    if (born_mapconfig(0) < 0 .or. &
        born_mapconfig(0) > ubound(born_mapconfig, 1)) then
      call fail_test_limits('generated Born configuration count is invalid')
    end if
  end subroutine validate_generated_data

  subroutine validate_work_state()
    implicit none

    if (.not. work_state_initialized) then
      call fail_test_limits('limit-test work state was not initialized')
    end if
    if (.not. allocated(p) .or. .not. allocated(x) .or. &
        .not. allocated(event_momenta) .or. &
        .not. allocated(event_jacobian) .or. &
        .not. allocated(amp_split)) then
      call fail_test_limits('limit-test work state is incomplete')
    end if
    if (size(p, 1) /= 4 .or. size(p, 2) /= nexternal .or. &
        size(x) /= random_vector_size .or. &
        size(event_momenta, 2) /= nexternal .or. &
        size(event_momenta, 3) /= real_event + 1 .or. &
        size(event_jacobian) /= real_event + 1 .or. &
        size(amp_split) /= amp_split_size) then
      call fail_test_limits('limit-test work state has invalid bounds')
    end if
  end subroutine validate_work_state

  subroutine finalize_test_soft_col_limits()
    implicit none

    if (allocated(born_mapconfig)) deallocate (born_mapconfig)
    if (allocated(generated_masses)) deallocate (generated_masses)
    generated_data_initialized = .false.

    if (allocated(fxl)) deallocate (fxl)
    if (allocated(wfxl)) deallocate (wfxl)
    if (allocated(limit)) deallocate (limit)
    if (allocated(wlimit)) deallocate (wlimit)
    if (allocated(lxp)) deallocate (lxp)
    if (allocated(xp)) deallocate (xp)
    if (allocated(p)) deallocate (p)
    if (allocated(x)) deallocate (x)
    if (allocated(fxl_split)) deallocate (fxl_split)
    if (allocated(wfxl_split)) deallocate (wfxl_split)
    if (allocated(limit_split)) deallocate (limit_split)
    if (allocated(wlimit_split)) deallocate (wlimit_split)
    if (allocated(orders)) deallocate (orders)
    if (allocated(nerr)) deallocate (nerr)
    if (allocated(fail_frac)) deallocate (fail_frac)
    if (allocated(event_momenta)) deallocate (event_momenta)
    if (allocated(event_jacobian)) deallocate (event_jacobian)
    if (allocated(p_born)) deallocate (p_born)
    if (allocated(event_fks_momentum)) deallocate (event_fks_momentum)
    if (allocated(event_xi)) deallocate (event_xi)
    if (allocated(event_y)) deallocate (event_y)
    if (allocated(etmin)) deallocate (etmin)
    if (allocated(etmax)) deallocate (etmax)
    if (allocated(mxxmin)) deallocate (mxxmin)
    if (allocated(amp_split)) deallocate (amp_split)
    work_state_initialized = .false.
  end subroutine finalize_test_soft_col_limits

  subroutine fail_test_limits(message)
    implicit none
    character(len=*), intent(in) :: message

    write (*, '(a)') 'test_soft_col_limits: '//trim(message)
    stop 1
  end subroutine fail_test_limits

end module test_soft_col_limits_module
