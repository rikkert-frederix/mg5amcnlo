module binoth_lha_madloop_backend
  use FKSParams
  use process_dimensions, only: nexternal, nincoming, max_bhel, &
       nsplitorders, amp_split_size, nlo_orders, order_names
  use split_orders, only: orders_to_amp_split_pos, &
       amp_split_pos_to_orders
  use fks_singular_module, only: getpoles
  implicit none
  private

  logical :: firsttime = .true.
  logical :: firsttime_conversion = .true.
  logical :: firsttime_run = .true.
  integer :: nbad = 0
  integer :: nsqso = 0
  integer :: ml_res_array_dim = 0
  logical, allocatable :: keep_order(:)
  double precision, allocatable :: virt_wgts(:,:)
  double precision, allocatable :: virt_wgts_hel(:,:)
  double precision, allocatable :: accuracies(:)
  double precision, allocatable :: amp_split_poles_ml(:,:)
  double precision, allocatable :: prec_found(:)
  character(len=1), allocatable :: include_hel(:)
  integer, allocatable :: goodhel(:), hel(:)

  double precision :: qes2
  common /coupl_es/ qes2

  integer :: ntot, nsun, nsps, nups, neps, n100
  integer :: nddp, nqdp, nini, n10, n1(0:9)
  common /ups_stats/ ntot, nsun, nsps, nups, neps, n100, nddp, nqdp, &
       nini, n10, n1

  double precision :: volh
  integer :: mc_hel, ihel
  logical :: fillh
  common /mc_int2/ volh, mc_hel, ihel, fillh

  logical :: force_polecheck, polecheck_passed
  common /to_polecheck/ force_polecheck, polecheck_passed

  integer :: ret_code_common
  common /to_ret_code/ ret_code_common

  logical :: cs_run
  common /to_cs_run/ cs_run

  interface
    subroutine binoth_lha_update_couplings(mu_r_value, alpha_s)
      implicit none
      double precision, intent(out) :: mu_r_value, alpha_s
    end subroutine binoth_lha_update_couplings

  end interface

  public :: binoth_lha_eval

contains

  subroutine binoth_lha_eval(p, born_wgt, virt_wgt, pmass, amp_split, &
       amp_split_finite_ml, amp_split_poles_fks)
    implicit none
    double precision, intent(in) :: p(0:3, nexternal-1)
    double precision, intent(inout) :: born_wgt
    double precision, intent(out) :: virt_wgt
    double precision, intent(in) :: pmass(nexternal)
    double precision, intent(inout) :: amp_split(amp_split_size)
    double precision, intent(inout) :: &
         amp_split_finite_ml(amp_split_size)
    double precision, intent(inout) :: &
         amp_split_poles_fks(amp_split_size,2)

    double precision, parameter :: pi = 3.1415926535897932385d0
    logical, parameter :: fksprefact = .true.
    integer, parameter :: nbadmax = 5
    double precision :: single_pole, double_pole
    double precision :: born_wgt_recomputed, born_wgt_recomp_direct
    double precision :: mu_r_value, ao2pi, conversion, alpha_s
    double precision :: madfks_single, madfks_double, tolerance
    double precision :: target, accum, hel_fact, born_hel_from_virt
    double precision :: avg_pole_res(2), pole_diff(2)
    integer :: ret_code, i, j, ioerr, ioerr_counter, dt(8)
    integer :: iamp
    integer :: order_name_width
    integer :: amp_orders(nsplitorders), split_amp_orders(nsplitorders)
    logical :: cpol
    integer, external :: getordpowfromindex_ml5
    double precision, external :: ran2
    ioerr_counter = 0
    order_name_width = maxval(len_trim(order_names))
    call binoth_lha_update_couplings(mu_r_value, alpha_s)
    ao2pi = alpha_s/(2d0*pi)

    virt_wgt = 0d0
    single_pole = 0d0
    double_pole = 0d0
    born_hel_from_virt = 0d0

    if (.not. allocated(amp_split_poles_ml)) then
      allocate(amp_split_poles_ml(amp_split_size,2))
      allocate(prec_found(amp_split_size))
    end if
    amp_split_finite_ml(1:amp_split_size) = 0d0
    amp_split_poles_ml(1:amp_split_size,1) = 0d0
    amp_split_poles_ml(1:amp_split_size,2) = 0d0
    prec_found(1:amp_split_size) = 0d0

    if (firsttime_run) then
      if (.not. force_polecheck) then
        call set_forbid_hel_doublecheck(.true.)
      end if
      call get_nsqso_loop(nsqso)
      call get_answer_dimension(ml_res_array_dim)
      allocate(accuracies(0:nsqso))
      allocate(virt_wgts(0:3,0:ml_res_array_dim))
      allocate(virt_wgts_hel(0:3,0:ml_res_array_dim))
      allocate(keep_order(nsqso))
      allocate(include_hel(max_bhel))
      allocate(goodhel(max_bhel))
      allocate(hel(0:max_bhel))
      call FORCE_STABILITY_CHECK(.true.)
      if (.not. force_polecheck) then
        call COLLIER_COMPUTE_UV_POLES(.false.)
        call COLLIER_COMPUTE_IR_POLES(.false.)
      else
        call COLLIER_COMPUTE_UV_POLES(.true.)
        call COLLIER_COMPUTE_IR_POLES(.true.)
      end if
      firsttime_run = .false.
    end if

    firsttime = firsttime .or. force_polecheck
    if (firsttime) then
      write (*,*) 'alpha_s value used for the virtuals' // &
           ' is (for the first PS point): ', alpha_s
      tolerance = IRPoleCheckThreshold/10d0
      call sloopmatrix_thres(p, virt_wgts, tolerance, accuracies, &
           ret_code)

      do i = 1, nsqso
        keep_order(i) = .true.
        do j = 1, nsplitorders
          if (getordpowfromindex_ml5(j,i) > nlo_orders(j)) then
            keep_order(i) = .false.
            exit
          end if
        end do
        if (keep_order(i)) then
          write (*,*) 'VIRT: keeping split order ', i
        else
          write (*,*) 'VIRT: not keeping split order ', i
        end if
      end do

      do i = 1, nsqso
        if (keep_order(i)) then
          virt_wgt = virt_wgt + virt_wgts(1,i)
          single_pole = single_pole + virt_wgts(2,i)
          double_pole = double_pole + virt_wgts(3,i)
          do j = 1, nsplitorders
            amp_orders(j) = getordpowfromindex_ml5(j,i)
          end do
          amp_split_finite_ml(orders_to_amp_split_pos(amp_orders)) = &
               virt_wgts(1,i)
          amp_split_poles_ml(orders_to_amp_split_pos(amp_orders),1) = &
               virt_wgts(2,i)
          amp_split_poles_ml(orders_to_amp_split_pos(amp_orders),2) = &
               virt_wgts(3,i)
          prec_found(orders_to_amp_split_pos(amp_orders)) = accuracies(i)
        end if
      end do
    else
      tolerance = PrecisionVirtualAtRunTime
      if (mc_hel == 0) then
        call sloopmatrix_thres(p, virt_wgts, tolerance, accuracies, &
             ret_code)
        do i = 1, nsqso
          if (keep_order(i)) then
            virt_wgt = virt_wgt + virt_wgts(1,i)
            single_pole = single_pole + virt_wgts(2,i)
            double_pole = double_pole + virt_wgts(3,i)
            do j = 1, nsplitorders
              amp_orders(j) = getordpowfromindex_ml5(j,i)
            end do
            amp_split_finite_ml(orders_to_amp_split_pos(amp_orders)) = &
                 virt_wgts(1,i)
            amp_split_poles_ml(orders_to_amp_split_pos(amp_orders),1) = &
                 virt_wgts(2,i)
            amp_split_poles_ml(orders_to_amp_split_pos(amp_orders),2) = &
                 virt_wgts(3,i)
            prec_found(orders_to_amp_split_pos(amp_orders)) = &
                 accuracies(i)
          end if
        end do
      else if (mc_hel == 1) then
        call PickHelicityMC(p, goodhel, hel, ihel, volh)
        fillh = .false.
        call sloopmatrixhel_thres(p, hel(ihel), virt_wgts_hel, &
             tolerance, accuracies, ret_code)
        hel_fact = dble(goodhel(ihel))/volh/4d0
        do i = 1, nsqso
          if (keep_order(i)) then
            born_hel_from_virt = born_hel_from_virt + virt_wgts_hel(0,i)
            virt_wgt = virt_wgt + virt_wgts_hel(1,i)*hel_fact
            single_pole = single_pole + virt_wgts_hel(2,i)*hel_fact
            double_pole = double_pole + virt_wgts_hel(3,i)*hel_fact
            do j = 1, nsplitorders
              amp_orders(j) = getordpowfromindex_ml5(j,i)
            end do
            amp_split_finite_ml(orders_to_amp_split_pos(amp_orders)) = &
                 virt_wgts_hel(1,i)*hel_fact
            amp_split_poles_ml(orders_to_amp_split_pos(amp_orders),1) = &
                 virt_wgts_hel(2,i)*hel_fact
            amp_split_poles_ml(orders_to_amp_split_pos(amp_orders),2) = &
                 virt_wgts_hel(3,i)*hel_fact
            prec_found(orders_to_amp_split_pos(amp_orders)) = &
                 accuracies(i)
          end if
        end do
        if (nincoming /= 2) then
          write (*,*) 'Cannot do MC over helicities for 1->N processes'
          stop
        end if
      else
        write (*,*) 'Can only do sum over helicities,' // &
             ' or pure MC over helicities', mc_hel
        stop
      end if
    end if

    if (cs_run) then
      print *, 'I am skipping checkpoles'
      return
    end if

    ! MadLoop already returns CDR virtuals, so no scheme conversion is applied.

    cpol = .false.
    ret_code_common = ret_code
    if ((firsttime .or. mc_hel == 0) .and. &
        mod(ret_code,100)/10 /= 3 .and. mod(ret_code,100)/10 /= 4) then
      call getpoles(p, qes2, madfks_double, madfks_single, fksprefact, &
           amp_split_poles_fks)
      polecheck_passed = .true.
      do iamp = 1, amp_split_size
        if (iamp /= 0) then
          if (amp_split_poles_fks(iamp,1) == 0d0 .and. &
              amp_split_poles_fks(iamp,1) == 0d0) cycle
        end if
        if (iamp == 0) then
          if (firsttime) then
            write (*,*) ''
            write (*,*) 'Sum of all split-orders'
          end if
        else
          if (firsttime) then
            write (*,*) ''
            write (*,*) 'Splitorders', iamp
            call amp_split_pos_to_orders(iamp, split_amp_orders)
            do i = 1, nsplitorders
              write (*,*) '      ', order_names(i)(1:order_name_width), ':', &
                   split_amp_orders(i)
            end do
          end if
          single_pole = amp_split_poles_ml(iamp,1)
          double_pole = amp_split_poles_ml(iamp,2)
          madfks_single = amp_split_poles_fks(iamp,1)
          madfks_double = amp_split_poles_fks(iamp,2)
        end if

        avg_pole_res(1) = (single_pole+madfks_single)/2d0
        avg_pole_res(2) = (double_pole+madfks_double)/2d0
        pole_diff(1) = dabs(single_pole-madfks_single)
        pole_diff(2) = dabs(double_pole-madfks_double)
        if (dabs(avg_pole_res(1))+dabs(avg_pole_res(2)) /= 0d0) then
          cpol = .not. ((pole_diff(1)+pole_diff(2))/ &
               (dabs(avg_pole_res(1))+dabs(avg_pole_res(2))) < &
               tolerance*10d0 .or. &
               (mod(ret_code,10) == 7 .and. .not. force_polecheck))
        else
          cpol = .not. (pole_diff(1)+pole_diff(2) < tolerance*10d0 .or. &
               mod(ret_code,10) == 7)
        end if
        if (tolerance < 0d0) cpol = .false.

        if (.not. cpol .and. firsttime) then
          write (*,*) '---- POLES CANCELLED ----'
          write (*,*) ' COEFFICIENT DOUBLE POLE:'
          write (*,*) '       MadFKS: ', madfks_double, &
               '          OLP: ', double_pole
          write (*,*) ' COEFFICIENT SINGLE POLE:'
          write (*,*) '       MadFKS: ', madfks_single, &
               '          OLP: ', single_pole
          if (iamp == 0) then
            write (*,*) ' FINITE:'
            write (*,*) '          OLP: ', virt_wgt
            write (*,*) '          BORN: ', born_wgt
            write (*,*) ' MOMENTA (Exyzm): '
            do i = 1, nexternal-1
              write (*,*) i, p(0,i), p(1,i), p(2,i), p(3,i), pmass(i)
            end do
          end if

          if (mc_hel /= 0) then
198         continue
            if (NHelForMCoverHels < 0) then
              mc_hel = 0
              go to 203
            end if
            open (unit=67, file='../MadLoop5_resources/HelFilter.dat', &
                 status='old', action='read', iostat=ioerr, err=201)
            hel(0) = 0
            j = 0
            do i = 1, max_bhel
              read (67, *, err=202) goodhel(i)
              if (goodhel(i) > -10000 .and. goodhel(i) /= 0) then
                j = j + 1
                goodhel(j) = goodhel(i)
                hel(0) = hel(0) + 1
                hel(j) = i
              end if
            end do
            go to 203

201         continue
            if (ioerr == 2 .and. ioerr_counter < 10) then
              ioerr_counter = ioerr_counter + 1
              write (*,*) 'File HelFilter.dat busy, retrying for' // &
                   ' the ', ioerr_counter, ' time.'
              call date_and_time(values=dt)
              call wait_retry_seconds(1+(dt(8)/200))
              go to 198
            end if
            write (*,*) 'Cannot do MC over hel:' // &
                 ' "HelFilter.dat" does not exist' // &
                 ' or does not have the correct format.' // &
                 ' Change NHelForMCoverHels in FKS_params.dat ' // &
                 'to explicitly summ over them instead.'
            stop

202         continue
            rewind(67)
            read (67, *, err=201) (include_hel(i),i=1,max_bhel)
            do i = 1, max_bhel
              if (include_hel(i) == 'T') then
                j = j + 1
                goodhel(j) = 1
                hel(0) = hel(0) + 1
                hel(j) = i
              end if
            end do

203         continue
            if (NHelForMCoverHels == -1) then
              write (*,*) 'Not doing MC over helicities: ' // &
                   'HelForMCoverHels=-1'
              mc_hel = 0
            else if (hel(0) < NHelForMCoverHels) then
              write (*,'(a,i3,a)') 'Only ', hel(0), &
                   ' independent helicities:' // &
                   ' switching to explicitly summing over them'
              mc_hel = 0
            end if
            close(67)
          end if
        else if (cpol .and. firsttime) then
          polecheck_passed = .false.
          write (*,*) 'POLES MISCANCELLATION, DIFFERENCE > ', &
               tolerance*10d0
          write (*,*) ' COEFFICIENT DOUBLE POLE:'
          write (*,*) '       MadFKS: ', madfks_double, &
               '          OLP: ', double_pole
          write (*,*) ' COEFFICIENT SINGLE POLE:'
          write (*,*) '       MadFKS: ', madfks_single, &
               '          OLP: ', single_pole
          if (iamp == 0) then
            write (*,*) ' FINITE:'
            write (*,*) '          OLP: ', virt_wgt
            write (*,*) '          BORN: ', born_wgt
            write (*,*) ' MOMENTA (Exyzm): '
            do i = 1, nexternal-1
              write (*,*) i, p(0,i), p(1,i), p(2,i), p(3,i), pmass(i)
            end do
          end if
          write (*,*)
          write (*,*) ' SCALE**2: ', qes2
          if (nbad < nbadmax) then
            nbad = nbad + 1
            write (*,*) ' Trying another PS point'
          else if (.not. force_polecheck) then
            write (*,*) ' TOO MANY FAILURES, QUITTING'
            stop
          end if
        end if
      end do
      firsttime = .false. .or. cpol
    end if

    ntot = ntot + 1
    if (ret_code/100 == 1) then
      nsun = nsun + 1
    else if (ret_code/100 == 2) then
      nsps = nsps + 1
    else if (ret_code/100 == 3) then
      nups = nups + 1
    else if (ret_code/100 == 4) then
      neps = neps + 1
    else
      n100 = n100 + 1
    end if

    if (mod(ret_code,100)/10 == 1 .or. &
        mod(ret_code,100)/10 == 3) then
      nddp = nddp + 1
      if (mod(ret_code,100)/10 == 3) nini = nini + 1
    else if (mod(ret_code,100)/10 == 2 .or. &
             mod(ret_code,100)/10 == 4) then
      nqdp = nqdp + 1
      if (mod(ret_code,100)/10 == 4) nini = nini + 1
    else
      n10 = n10 + 1
    end if
    n1(mod(ret_code,10)) = n1(mod(ret_code,10)) + 1

    do iamp = 1, amp_split_size
      if (.not. firsttime .and. (ret_code/100 == 4 .or. cpol .or. &
          prec_found(iamp) > 0.05d0 .or. &
          amp_split_finite_ml(iamp) /= amp_split_finite_ml(iamp))) then
        if (neps < 10) then
          if (neps == 1) then
            open (unit=78, file='UPS.log')
          else
            open (unit=78, file='UPS.log', status='unknown', &
                 position='append', action='write')
          end if
          write (78,*) '===== EPS #', neps, ' ====='
          write (78,*) 'mu_r    =', mu_r_value
          write (78,*) 'alpha_S =', alpha_s
          write (78,*) 'MadLoop return code, pole check and' // &
               ' accuracy reported', ret_code, cpol, prec_found
          if (mc_hel /= 0) then
            write (78,*) 'helicity (MadLoop only)', hel(i), mc_hel
          end if
          write (78,*) '1/eps**2 expected from MadFKS=', &
               amp_split_poles_fks(iamp,2)
          write (78,*) '1/eps**2 obtained in MadLoop =', &
               amp_split_poles_ml(iamp,2)
          write (78,*) '1/eps    expected from MadFKS=', &
               amp_split_poles_fks(iamp,1)
          write (78,*) '1/eps    obtained in MadLoop =', &
               amp_split_poles_ml(iamp,1)
          write (78,*) 'finite   obtained in MadLoop =', &
               amp_split_finite_ml(iamp)
          write (78,*) 'Accuracy estimated by MadLop =', prec_found(iamp)
          do i = 1, nexternal-1
            write (78,'(i2,1x,5e25.15)') i, p(0,i), p(1,i), &
                 p(2,i), p(3,i), pmass(i)
          end do
          close(78)
        end if
        if (prec_found(iamp) > 0.05d0 .or. &
            amp_split_finite_ml(iamp) /= amp_split_finite_ml(iamp)) then
          write (*,*) 'WARNING: unstable non-rescued phase-space' // &
               ' found for which the accuracy reported by' // &
               ' MadLoop is worse than 5%. Setting virtual to' // &
               ' zero for this PS point.'
          amp_split_finite_ml(iamp) = 0d0
        end if
      end if
    end do

    if (.not. firsttime .and. &
        (accuracies(0) > 0.05d0 .or. virt_wgt /= virt_wgt)) then
      virt_wgt = 0d0
    end if
    if ((mod(ret_code,100)/10 == 4 .or. &
         mod(ret_code,100)/10 == 3) .and. ret_code/100 == 1) then
      do iamp = 1, amp_split_size
        amp_split_finite_ml(iamp) = 0d0
      end do
      virt_wgt = 0d0
    end if
  end subroutine binoth_lha_eval


  subroutine binoth_lha_init_impl(filename)
    implicit none
    character(len=13), intent(in) :: filename

    ! The Rocket and BlackHat initialization examples remain intentionally
    ! inactive, as in the legacy MadLoop-selected backend.
  end subroutine binoth_lha_init_impl


  subroutine dr_to_cdr_impl(conversion, i_fks, j_fks, particle_type, &
       m_type, pmass)
    implicit none
    double precision, intent(out) :: conversion
    integer, intent(in) :: i_fks, j_fks
    integer, intent(in) :: particle_type(nexternal)
    integer, intent(in) :: m_type
    double precision, intent(in) :: pmass(nexternal)

    double precision, parameter :: ca = 3d0
    double precision, parameter :: cf = 4d0/3d0
    integer :: i, triplet, octet

    triplet = 0
    octet = 0
    conversion = 0d0
    do i = 1, nexternal
      if (i /= i_fks .and. i /= j_fks) then
        if (pmass(i) == 0d0) then
          if (abs(particle_type(i)) == 3) then
            conversion = conversion - cf/2d0
            triplet = triplet + 1
          else if (particle_type(i) == 8) then
            conversion = conversion - ca/6d0
            octet = octet + 1
          end if
        end if
      else if (i == min(i_fks,j_fks)) then
        if (pmass(j_fks) == 0d0 .and. pmass(i_fks) == 0d0) then
          if (m_type == 8) then
            conversion = conversion - ca/6d0
            octet = octet + 1
          else if (abs(m_type) == 3) then
            conversion = conversion - cf/2d0
            triplet = triplet + 1
          else
            write (*,*) 'Error in DRtoCDR, fks_mother must be' // &
                 'triplet or octet', i, m_type
            stop
          end if
        end if
      end if
    end do
    write (*,*) 'From DR to CDR conversion: ', octet, ' octets and ', &
         triplet, ' triplets in Born (both massless), sum =', conversion
  end subroutine dr_to_cdr_impl


  subroutine wait_retry_seconds(seconds)
    implicit none
    integer, intent(in) :: seconds
    integer :: start_count, current_count, count_rate

    if (seconds <= 0) return
    call system_clock(start_count, count_rate)
    if (count_rate <= 0) return
    do
      call system_clock(current_count)
      if (current_count < start_count) exit
      if (current_count-start_count >= seconds*count_rate) exit
    end do
  end subroutine wait_retry_seconds


  subroutine get_procnum_impl(filename, procnum)
    implicit none
    character(len=13), intent(in) :: filename
    integer, intent(out) :: procnum
    integer :: lookhere, procsize
    character(len=176) :: buff
    logical :: done

    open (unit=68, file=filename, status='old')
    done = .false.
    do while (.not. done)
      read (68, '(a)', end=889) buff
      if (index(buff, '->') /= 0) then
        if (lookhere /= 0 .and. lookhere < 170) then
          write (*,*) 'Read process number from contract file', procnum
          close(68)
          return
          done = .true.
        else
          write (*,*) 'syntax contract file not understandable', lookhere
          stop
        end if
      end if
    end do
    stop

    close(68)
    return

889 write (*,*) 'Error in contract file'
    stop
  end subroutine get_procnum_impl

end module binoth_lha_madloop_backend
