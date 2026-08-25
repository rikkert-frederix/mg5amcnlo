module cernlib_module
  implicit none
  private

  integer, parameter :: dp = selected_real_kind(15, 307)
  integer, parameter :: mtl_table_size = 132
  integer, parameter :: ker_table_size = 28

  character(len=6), allocatable :: mtl_codes(:)
  integer, allocatable :: mtl_message_limits(:)
  integer, allocatable :: mtl_return_limits(:)
  integer :: mtl_log_unit = 0

  character(len=6), allocatable :: ker_codes(:)
  integer, allocatable :: ker_message_limits(:)
  integer, allocatable :: ker_return_limits(:)
  integer :: ker_log_unit = 0

  public :: initialize_cernlib
  public :: finalize_cernlib
  public :: abend
  public :: lenocc
  public :: dlsqp2
  public :: dfint
  public :: mtlprt
  public :: mtlset
  public :: mtlmtr
  public :: kerset
  public :: kermtr
  public :: radmul
  public :: dadmul

contains

  subroutine initialize_cernlib()
    implicit none

    call initialize_mtl_monitor()
    call initialize_ker_monitor()
  end subroutine initialize_cernlib


  subroutine finalize_cernlib()
    implicit none

    if (allocated(mtl_codes)) deallocate(mtl_codes)
    if (allocated(mtl_message_limits)) deallocate(mtl_message_limits)
    if (allocated(mtl_return_limits)) deallocate(mtl_return_limits)
    if (allocated(ker_codes)) deallocate(ker_codes)
    if (allocated(ker_message_limits)) deallocate(ker_message_limits)
    if (allocated(ker_return_limits)) deallocate(ker_return_limits)
    mtl_log_unit = 0
    ker_log_unit = 0
  end subroutine finalize_cernlib


  subroutine initialize_mtl_monitor()
    implicit none

    if (allocated(mtl_codes)) return

    allocate(mtl_codes(mtl_table_size))
    allocate(mtl_message_limits(mtl_table_size))
    allocate(mtl_return_limits(mtl_table_size))

    mtl_codes = (/ &
      'B100.1', 'B300.1', 'B300.2', 'C200.0', 'C200.1', 'C200.2', &
      'C200.3', 'C201.0', 'C202.0', 'C202.1', 'C202.2', 'C205.1', &
      'C205.2', 'C207.0', 'C208.0', 'C209.0', 'C209.1', 'C209.2', &
      'C209.3', 'C210.1', 'C302.1', 'C303.1', 'C304.1', 'C305.1', &
      'C306.1', 'C307.1', 'C312.1', 'C313.1', 'C315.1', 'C316.1', &
      'C316.2', 'C320.1', 'C321.1', 'C323.1', 'C327.1', 'C328.1', &
      'C328.2', 'C328.3', 'C330.1', 'C330.2', 'C330.3', 'C331.1', &
      'C331.2', 'C334.1', 'C334.2', 'C334.3', 'C334.4', 'C334.5', &
      'C334.6', 'C336.1', 'C337.1', 'C338.1', 'C340.1', 'C343.1', &
      'C343.2', 'C343.3', 'C343.4', 'C344.1', 'C344.2', 'C344.3', &
      'C344.4', 'C345.1', 'C346.1', 'C346.2', 'C346.3', 'C347.1', &
      'C347.2', 'C347.3', 'C347.4', 'C347.5', 'C347.6', 'C348.1', &
      'C349.1', 'C349.2', 'C349.3', 'D101.1', 'D103.1', 'D104.1', &
      'D104.2', 'D105.1', 'D105.2', 'D107.1', 'D110.0', 'D110.1', &
      'D110.2', 'D110.3', 'D110.4', 'D110.5', 'D110.6', 'D113.1', &
      'D201.1', 'D202.1', 'D401.1', 'D601.1', 'E210.1', 'E210.2', &
      'E210.3', 'E210.4', 'E210.5', 'E210.6', 'E210.7', 'E211.0', &
      'E211.1', 'E211.2', 'E211.3', 'E211.4', 'E406.0', 'E406.1', &
      'E407.0', 'E408.0', 'E408.1', 'F500.0', 'F500.1', 'F500.2', &
      'F500.3', 'G100.1', 'G100.2', 'G101.1', 'G101.2', 'G105.1', &
      'G106.1', 'G106.2', 'G116.1', 'G116.2', 'H101.0', 'H101.1', &
      'H101.2', 'H301.1', 'U501.1', 'V202.1', 'V202.2', 'V202.3' /)
    mtl_message_limits = 255
    mtl_return_limits = 255
    mtl_log_unit = 0
  end subroutine initialize_mtl_monitor


  subroutine initialize_ker_monitor()
    implicit none

    if (allocated(ker_codes)) return

    allocate(ker_codes(ker_table_size))
    allocate(ker_message_limits(ker_table_size))
    allocate(ker_return_limits(ker_table_size))

    ker_codes = (/ &
      'C204.1', 'C204.2', 'C204.3', 'C205.1', 'C205.2', 'C205.3', &
      'C305.1', 'C308.1', 'C312.1', 'C313.1', 'C336.1', 'C337.1', &
      'C341.1', 'D103.1', 'D106.1', 'D209.1', 'D509.1', 'E100.1', &
      'E104.1', 'E105.1', 'E208.1', 'E208.2', 'F010.1', 'F011.1', &
      'F012.1', 'F406.1', 'G100.1', 'G100.2' /)
    ker_message_limits = 100
    ker_return_limits = 100
    ker_return_limits(23:26) = 0
    ker_log_unit = 0
  end subroutine initialize_ker_monitor


  subroutine abend()
    implicit none

    stop 7
  end subroutine abend


  integer function lenocc(chv) result(last_occupied)
    implicit none
    character(len=*), intent(in) :: chv
    integer :: jj

    do jj = len(chv), 1, -1
      if (chv(jj:jj) /= ' ') then
        last_occupied = jj
        return
      end if
    end do
    last_occupied = 0
  end function lenocc


  subroutine dlsqp2(n, x, y, a0, a1, a2, sd, ifail)
    implicit none
    integer, intent(in) :: n
    real(kind=dp), intent(in) :: x(*)
    real(kind=dp), intent(in) :: y(*)
    real(kind=dp), intent(out) :: a0
    real(kind=dp), intent(out) :: a1
    real(kind=dp), intent(out) :: a2
    real(kind=dp), intent(out) :: sd
    integer, intent(out) :: ifail

    real(kind=dp) :: determinant
    real(kind=dp) :: fn
    real(kind=dp) :: sxx
    real(kind=dp) :: sxxx
    real(kind=dp) :: sxxxx
    real(kind=dp) :: sy
    real(kind=dp) :: syy
    real(kind=dp) :: sxy
    real(kind=dp) :: sxxy
    real(kind=dp) :: xk
    real(kind=dp) :: xk2
    real(kind=dp) :: xm
    real(kind=dp) :: yk
    integer :: k

    a0 = 0.0_dp
    a1 = 0.0_dp
    a2 = 0.0_dp
    sd = 0.0_dp
    xm = 0.0_dp

    if (n <= 2) then
      ifail = 1
    else
      fn = real(n, kind=dp)
      do k = 1, n
        xm = xm + x(k)
      end do
      xm = xm / fn

      sxx = 0.0_dp
      sxxx = 0.0_dp
      sxxxx = 0.0_dp
      sy = 0.0_dp
      syy = 0.0_dp
      sxy = 0.0_dp
      sxxy = 0.0_dp
      do k = 1, n
        xk = x(k) - xm
        yk = y(k)
        xk2 = xk**2
        sxx = sxx + xk2
        sxxx = sxxx + xk2 * xk
        sxxxx = sxxxx + xk2**2
        sy = sy + yk
        syy = syy + yk**2
        sxy = sxy + xk * yk
        sxxy = sxxy + xk2 * yk
      end do

      determinant = (fn * sxxxx - sxx**2) * sxx - fn * sxxx**2
      if (determinant > 0.0_dp) then
        a2 = (sxx * (fn * sxxy - sxx * sy) - &
              fn * sxxx * sxy) / determinant
        a1 = (sxy - sxxx * a2) / sxx
        a0 = (sy - sxx * a2) / fn
        ifail = 0
      else
        ifail = -1
      end if
    end if

    if (ifail == 0 .and. n > 3) then
      sd = sqrt(max(0.0_dp, syy - a0 * sy - a1 * sxy - a2 * sxxy) / &
                real(n - 3, kind=dp))
    end if
    a0 = a0 + xm * (xm * a2 - a1)
    a1 = a1 - 2.0_dp * xm * a2
  end subroutine dlsqp2


  real(kind=dp) function dfint(narg, arg, nent, ent, table)
    implicit none
    integer, intent(in) :: narg
    real(kind=dp), intent(in) :: arg(*)
    integer, intent(in) :: nent(*)
    real(kind=dp), intent(in) :: ent(*)
    real(kind=dp), intent(in) :: table(*)

    integer :: table_index(32)
    real(kind=dp) :: weight(32)
    real(kind=dp) :: eta
    real(kind=dp) :: h
    real(kind=dp) :: x
    integer :: i
    integer :: ishift
    integer :: istep
    integer :: k
    integer :: knots
    integer :: loca
    integer :: locb
    integer :: locc
    integer :: lmax
    integer :: lmin
    integer :: ndim
    integer :: log_unit
    logical :: exact_knot
    logical :: message_flag
    logical :: return_flag

    dfint = 0.0_dp
    if (narg < 1 .or. narg > 5) then
      call kermtr('E104.1', log_unit, message_flag, return_flag)
      if (message_flag) then
        if (log_unit == 0) then
          write(*, 1000) narg
        else
          write(log_unit, 1000) narg
        end if
      end if
      if (.not. return_flag) call abend()
      return
    end if

    lmax = 0
    istep = 1
    knots = 1
    table_index(1) = 1
    weight(1) = 1.0_dp

    do i = 1, narg
      x = arg(i)
      ndim = nent(i)
      loca = lmax
      lmin = lmax + 1
      lmax = lmax + ndim
      exact_knot = .false.

      if (ndim <= 2) then
        if (ndim == 1) then
          exact_knot = .true.
          ishift = 0
        else
          h = x - ent(lmin)
          if (h == 0.0_dp) then
            exact_knot = .true.
            ishift = 0
          else if (x - ent(lmin + 1) == 0.0_dp) then
            exact_knot = .true.
            ishift = istep
          else
            ishift = 0
            eta = h / (ent(lmin + 1) - ent(lmin))
          end if
        end if
      else
        locb = lmax + 1
        do
          locc = (loca + locb) / 2
          if (x < ent(locc)) then
            locb = locc
          else if (x == ent(locc)) then
            exact_knot = .true.
            ishift = (locc - lmin) * istep
            exit
          else
            loca = locc
          end if
          if (locb - loca <= 1) exit
        end do
        if (.not. exact_knot) then
          loca = min(max(loca, lmin), lmax - 1)
          ishift = (loca - lmin) * istep
          eta = (x - ent(loca)) / (ent(loca + 1) - ent(loca))
        end if
      end if

      if (exact_knot) then
        do k = 1, knots
          table_index(k) = table_index(k) + ishift
        end do
      else
        do k = 1, knots
          table_index(k) = table_index(k) + ishift
          table_index(k + knots) = table_index(k) + istep
          weight(k + knots) = weight(k) * eta
          weight(k) = weight(k) - weight(k + knots)
        end do
        knots = 2 * knots
      end if
      istep = istep * ndim
    end do

    do k = 1, knots
      dfint = dfint + weight(k) * table(table_index(k))
    end do
1000 format(7x, 'FUNCTION DFINT ... NARG =', i6, ' NOT WITHIN RANGE')
  end function dfint


  subroutine mtlset(error_code, log_unit, message_limit, return_limit)
    implicit none
    character(len=*), intent(in) :: error_code
    integer, intent(in) :: log_unit
    integer, intent(in) :: message_limit
    integer, intent(in) :: return_limit

    integer :: code_length
    integer :: i
    logical :: matches_code

    call initialize_mtl_monitor()
    mtl_log_unit = log_unit
    code_length = min(len_trim(error_code), len(mtl_codes(1)))

    do i = 1, size(mtl_codes)
      if (code_length == 0) then
        matches_code = .true.
      else
        matches_code = &
          mtl_codes(i)(1:code_length) == error_code(1:code_length)
      end if
      if (matches_code) then
        if (message_limit >= 0) mtl_message_limits(i) = message_limit
        if (return_limit >= 0) mtl_return_limits(i) = return_limit
      end if
    end do
  end subroutine mtlset


  subroutine mtlmtr(error_code, log_unit, message_flag, return_flag)
    implicit none
    character(len=*), intent(in) :: error_code
    integer, intent(out) :: log_unit
    logical, intent(out) :: message_flag
    logical, intent(out) :: return_flag

    integer :: i

    call initialize_mtl_monitor()
    log_unit = mtl_log_unit

    do i = 1, size(mtl_codes)
      if (error_code == mtl_codes(i)) then
        message_flag = mtl_message_limits(i) >= 1
        return_flag = mtl_return_limits(i) >= 1
        if (message_flag .and. mtl_message_limits(i) < 255) then
          mtl_message_limits(i) = mtl_message_limits(i) - 1
        end if
        if (return_flag .and. mtl_return_limits(i) < 255) then
          mtl_return_limits(i) = mtl_return_limits(i) - 1
        end if
        if (.not. return_flag) then
          if (mtl_log_unit < 1) then
            write(*, 101) mtl_codes(i)
          else
            write(mtl_log_unit, 101) mtl_codes(i)
          end if
        end if
        return
      end if
    end do

    write(*, 100) error_code
    call abend()

100 format(7x, '***** CERN N002 MTLSET ... ERROR N002: ', &
           'ERROR CODE ', a6, &
           ' NOT RECOGNIZED BY ERROR MONITOR. RUN ABORTED.')
101 format(7x, '***** CERN N002 MTLSET ... ERROR NOO2.1: ', &
           'RUN TERMINATED BY LIBRARY ERROR CONDITION ', a6)
  end subroutine mtlmtr


  subroutine kerset(error_code, log_unit, message_limit, return_limit)
    implicit none
    character(len=*), intent(in) :: error_code
    integer, intent(in) :: log_unit
    integer, intent(in) :: message_limit
    integer, intent(in) :: return_limit

    integer :: code_length
    integer :: i
    logical :: matches_code

    call initialize_ker_monitor()
    ker_log_unit = log_unit
    code_length = min(len_trim(error_code), len(ker_codes(1)))

    do i = 1, size(ker_codes)
      if (code_length == 0) then
        matches_code = .true.
      else
        matches_code = &
          ker_codes(i)(1:code_length) == error_code(1:code_length)
      end if
      if (matches_code) then
        ker_message_limits(i) = message_limit
        ker_return_limits(i) = return_limit
      end if
    end do
  end subroutine kerset


  subroutine kermtr(error_code, log_unit, message_flag, return_flag)
    implicit none
    character(len=*), intent(in) :: error_code
    integer, intent(out) :: log_unit
    logical, intent(out) :: message_flag
    logical, intent(out) :: return_flag

    integer :: i

    call initialize_ker_monitor()
    log_unit = ker_log_unit

    do i = 1, size(ker_codes)
      if (error_code == ker_codes(i)) then
        return_flag = ker_return_limits(i) >= 1
        if (return_flag .and. ker_return_limits(i) < 100) then
          ker_return_limits(i) = ker_return_limits(i) - 1
        end if
        message_flag = ker_message_limits(i) >= 1
        if (message_flag .and. ker_message_limits(i) < 100) then
          ker_message_limits(i) = ker_message_limits(i) - 1
        end if

        if (.not. return_flag) then
          if (ker_log_unit < 1) then
            write(*, 1001) ker_codes(i)
          else
            write(ker_log_unit, 1001) ker_codes(i)
          end if
        end if
        if (message_flag .and. return_flag) then
          if (ker_log_unit < 1) then
            write(*, 1002) ker_codes(i)
          else
            write(ker_log_unit, 1002) ker_codes(i)
          end if
        end if
        return
      end if
    end do

    write(*, 1000) error_code
    call abend()

1000 format(' KERNLIB LIBRARY ERROR. ' / &
            ' ERROR CODE ', a6, ' NOT RECOGNIZED BY KERMTR', &
            ' ERROR MONITOR. RUN ABORTED.')
1001 format(/ ' ***** RUN TERMINATED BY CERN LIBRARY ERROR ', &
            'CONDITION ', a6)
1002 format(/ ' ***** CERN LIBRARY ERROR CONDITION ', a6)
  end subroutine kermtr


  subroutine mtlprt(name, error_code, text)
    implicit none
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: error_code
    character(len=*), intent(in) :: text

    character(len=6) :: normalized_code
    integer :: log_unit
    integer :: occupied_length
    logical :: message_flag
    logical :: return_flag

    normalized_code = '      '
    normalized_code(1:min(len(error_code), 6)) = &
      error_code(1:min(len(error_code), 6))

    if (normalized_code(5:6) /= '.0') then
      call mtlmtr(normalized_code, log_unit, message_flag, return_flag)
    else
      call initialize_mtl_monitor()
      log_unit = mtl_log_unit
      message_flag = .true.
      return_flag = .false.
    end if

    if (message_flag) then
      occupied_length = lenocc(text)
      if (log_unit < 1) then
        write(*, 100) normalized_code(1:4), name, normalized_code, &
                      text(1:occupied_length)
      else
        write(log_unit, 100) normalized_code(1:4), name, &
                             normalized_code, text(1:occupied_length)
      end if
    end if
    if (.not. return_flag) call abend()

100 format(7x, '***** CERN ', a, 1x, a, ' ERROR ', a, ': ', a)
  end subroutine mtlprt


  subroutine radmul(f, n, a, b, minpts, maxpts, eps, wk, iwk, result, &
                    relerr, nfnevl, ifail)
    implicit none
    real, external :: f
    integer, intent(in) :: n
    real, intent(in) :: a(*)
    real, intent(in) :: b(*)
    integer, intent(in) :: minpts
    integer, intent(in) :: maxpts
    real, intent(in) :: eps
    real, intent(inout) :: wk(*)
    integer, intent(in) :: iwk
    real, intent(out) :: result
    real, intent(out) :: relerr
    integer, intent(out) :: nfnevl
    integer, intent(out) :: ifail

    call mtlprt('RADMUL', 'D120.0', &
                'not available on this machine - see documentation')
  end subroutine radmul


  subroutine dadmul(f, n, a, b, minpts, maxpts, eps, wk, iwk, result, &
                    relerr, nfnevl, ifail)
    implicit none
    interface
      function f(n, z) result(value)
        integer, intent(in) :: n
        double precision, intent(in) :: z(*)
        double precision :: value
      end function f
    end interface
    integer, intent(in) :: n
    real(kind=dp), intent(in) :: a(*)
    real(kind=dp), intent(in) :: b(*)
    integer, intent(in) :: minpts
    integer, intent(in) :: maxpts
    real(kind=dp), intent(in) :: eps
    real(kind=dp), intent(inout) :: wk(*)
    integer, intent(in) :: iwk
    real(kind=dp), intent(out) :: result
    real(kind=dp), intent(out) :: relerr
    integer, intent(out) :: nfnevl
    integer, intent(out) :: ifail

    real(kind=dp), parameter :: half = 0.5_dp
    real(kind=dp), parameter :: xl2 = 0.358568582800318073_dp
    real(kind=dp), parameter :: xl4 = 0.948683298050513796_dp
    real(kind=dp), parameter :: xl5 = 0.688247201611685289_dp
    real(kind=dp), parameter :: w2 = 980.0_dp / 6561.0_dp
    real(kind=dp), parameter :: w4 = 200.0_dp / 19683.0_dp
    real(kind=dp), parameter :: wp2 = 245.0_dp / 486.0_dp
    real(kind=dp), parameter :: wp4 = 25.0_dp / 729.0_dp

    real(kind=dp) :: absolute_error
    real(kind=dp) :: center(15)
    real(kind=dp) :: difference
    real(kind=dp) :: difference_max
    real(kind=dp) :: f2
    real(kind=dp) :: f3
    real(kind=dp) :: region_comparison
    real(kind=dp) :: region_error
    real(kind=dp) :: region_value
    real(kind=dp) :: region_volume
    real(kind=dp) :: sum1
    real(kind=dp) :: sum2
    real(kind=dp) :: sum3
    real(kind=dp) :: sum4
    real(kind=dp) :: sum5
    real(kind=dp) :: two_to_n
    real(kind=dp) :: w(2:15, 5)
    real(kind=dp) :: width(15)
    real(kind=dp) :: width_offset(15)
    real(kind=dp) :: wp(2:15, 3)
    real(kind=dp) :: z(15)
    integer :: axis
    integer :: child_slot
    integer :: evaluation_count
    integer :: evaluations_per_region
    integer :: heap_end
    integer :: heap_slot
    integer :: i
    integer :: j
    integer :: j1
    integer :: k
    integer :: l
    integer :: m
    integer :: region_stride
    integer :: selected_axis
    logical :: divide_region

    w = 0.0_dp
    wp = 0.0_dp
    w(2:15, 1) = (/ &
      -0.193872885230909911_dp, -0.555606360818980835_dp, &
      -0.876695625666819078_dp, -1.15714067977442459_dp, &
      -1.39694152314179743_dp, -1.59609815576893754_dp, &
      -1.75461057765584494_dp, -1.87247878880251983_dp, &
      -1.94970278920896201_dp, -1.98628257887517146_dp, &
      -1.98221815780114818_dp, -1.93750952598689219_dp, &
      -1.85215668343240347_dp, -1.72615963013768225_dp /)
    w(2:15, 3) = (/ &
       0.0518213686937966768_dp,  0.0314992633236803330_dp, &
       0.0111771579535639891_dp, -0.00914494741655235473_dp, &
      -0.0294670527866686986_dp, -0.0497891581567850424_dp, &
      -0.0701112635269013768_dp, -0.0904333688970177241_dp, &
      -0.110755474267134071_dp,  -0.131077579637250419_dp, &
      -0.151399685007366752_dp,  -0.171721790377483099_dp, &
      -0.192043895747599447_dp,  -0.212366001117715794_dp /)
    w(2:15, 5) = (/ &
      0.0871183254585174982_dp, 0.0435591627292587508_dp, &
      0.0217795813646293754_dp, 0.0108897906823146873_dp, &
      0.00544489534115734364_dp, 0.00272244767057867193_dp, &
      0.00136122383528933596_dp, 0.000680611917644667955_dp, &
      0.000340305958822333977_dp, 0.000170152979411166995_dp, &
      0.0000850764897055834977_dp, 0.0000425382448527917472_dp, &
      0.0000212691224263958736_dp, 0.0000106345612131979372_dp /)
    wp(2:15, 1) = (/ &
      -1.33196159122085045_dp, -2.29218106995884763_dp, &
      -3.11522633744855959_dp, -3.80109739368998611_dp, &
      -4.34979423868312742_dp, -4.76131687242798352_dp, &
      -5.03566529492455417_dp, -5.17283950617283939_dp, &
      -5.17283950617283939_dp, -5.03566529492455417_dp, &
      -4.76131687242798352_dp, -4.34979423868312742_dp, &
      -3.80109739368998611_dp, -3.11522633744855959_dp /)
    wp(2:15, 3) = (/ &
       0.0445816186556927292_dp, -0.0240054869684499309_dp, &
      -0.0925925925925925875_dp, -0.161179698216735251_dp, &
      -0.229766803840877915_dp,  -0.298353909465020564_dp, &
      -0.366941015089163228_dp,  -0.435528120713305891_dp, &
      -0.504115226337448555_dp,  -0.572702331961591218_dp, &
      -0.641289437585733882_dp,  -0.709876543209876532_dp, &
      -0.778463648834019195_dp,  -0.847050754458161859_dp /)

    result = 0.0_dp
    absolute_error = 0.0_dp
    ifail = 3
    nfnevl = 0
    relerr = 0.0_dp
    if (n < 2 .or. n > 15) return
    if (minpts > maxpts) return

    evaluation_count = 0
    divide_region = .false.
    axis = 1
    two_to_n = real(2**n, kind=dp)
    region_stride = 2 * n + 3
    evaluations_per_region = 2**n + 2 * n * (n + 1) + 1
    heap_slot = region_stride
    heap_end = region_stride
    if (maxpts < evaluations_per_region) return
    do j = 1, n
      center(j) = (b(j) + a(j)) * half
      width(j) = (b(j) - a(j)) * half
    end do

20  continue
    region_volume = two_to_n
    do j = 1, n
      region_volume = region_volume * width(j)
      z(j) = center(j)
    end do
    sum1 = f(n, z)

    difference_max = 0.0_dp
    sum2 = 0.0_dp
    sum3 = 0.0_dp
    do j = 1, n
      z(j) = center(j) - xl2 * width(j)
      f2 = f(n, z)
      z(j) = center(j) + xl2 * width(j)
      f2 = f2 + f(n, z)
      width_offset(j) = xl4 * width(j)
      z(j) = center(j) - width_offset(j)
      f3 = f(n, z)
      z(j) = center(j) + width_offset(j)
      f3 = f3 + f(n, z)
      sum2 = sum2 + f2
      sum3 = sum3 + f3
      difference = abs(7.0_dp * f2 - f3 - 12.0_dp * sum1)
      if (difference >= difference_max) then
        difference_max = difference
        axis = j
      end if
      z(j) = center(j)
    end do

    sum4 = 0.0_dp
    do j = 2, n
      j1 = j - 1
      do k = j, n
        do l = 1, 2
          width_offset(j1) = -width_offset(j1)
          z(j1) = center(j1) + width_offset(j1)
          do m = 1, 2
            width_offset(k) = -width_offset(k)
            z(k) = center(k) + width_offset(k)
            sum4 = sum4 + f(n, z)
          end do
        end do
        z(k) = center(k)
      end do
      z(j1) = center(j1)
    end do

    sum5 = 0.0_dp
    do j = 1, n
      width_offset(j) = -xl5 * width(j)
      z(j) = center(j) + width_offset(j)
    end do
90  continue
    sum5 = sum5 + f(n, z)
    do j = 1, n
      width_offset(j) = -width_offset(j)
      z(j) = center(j) + width_offset(j)
      if (width_offset(j) > 0.0_dp) goto 90
    end do

    region_comparison = region_volume * &
      (wp(n, 1) * sum1 + wp2 * sum2 + wp(n, 3) * sum3 + wp4 * sum4)
    region_value = w(n, 1) * sum1 + w2 * sum2 + &
                   w(n, 3) * sum3 + w4 * sum4 + w(n, 5) * sum5
    region_value = region_volume * region_value
    region_error = abs(region_value - region_comparison)
    result = result + region_value
    absolute_error = absolute_error + region_error
    evaluation_count = evaluation_count + evaluations_per_region

    if (divide_region) then
110   continue
      child_slot = 2 * heap_slot
      if (child_slot > heap_end) goto 160
      if (child_slot < heap_end) then
        i = child_slot + region_stride
        if (wk(child_slot) < wk(i)) child_slot = i
      end if
      if (region_error >= wk(child_slot)) goto 160
      do k = 0, region_stride - 1
        wk(heap_slot - k) = wk(child_slot - k)
      end do
      heap_slot = child_slot
      goto 110
    end if

140 continue
    child_slot = (heap_slot / (2 * region_stride)) * region_stride
    if (child_slot >= region_stride) then
      if (region_error > wk(child_slot)) then
        do k = 0, region_stride - 1
          wk(heap_slot - k) = wk(child_slot - k)
        end do
        heap_slot = child_slot
        goto 140
      end if
    end if

160 continue
    wk(heap_slot) = region_error
    wk(heap_slot - 1) = region_value
    wk(heap_slot - 2) = real(axis, kind=dp)
    do j = 1, n
      child_slot = heap_slot - 2 * j - 2
      wk(child_slot + 1) = center(j)
      wk(child_slot) = width(j)
    end do
    if (divide_region) then
      divide_region = .false.
      center(selected_axis) = center(selected_axis) + &
                              2.0_dp * width(selected_axis)
      heap_end = heap_end + region_stride
      heap_slot = heap_end
      goto 20
    end if

    relerr = absolute_error / abs(result)
    if (heap_end + region_stride > iwk) ifail = 2
    if (evaluation_count + 2 * evaluations_per_region > maxpts) ifail = 1
    if (relerr < eps .and. evaluation_count >= minpts) ifail = 0
    if (ifail == 3) then
      divide_region = .true.
      heap_slot = region_stride
      absolute_error = absolute_error - wk(heap_slot)
      result = result - wk(heap_slot - 1)
      selected_axis = int(wk(heap_slot - 2))
      do j = 1, n
        child_slot = heap_slot - 2 * j - 2
        center(j) = wk(child_slot + 1)
        width(j) = wk(child_slot)
      end do
      width(selected_axis) = half * width(selected_axis)
      center(selected_axis) = center(selected_axis) - width(selected_axis)
      goto 20
    end if
    nfnevl = evaluation_count
  end subroutine dadmul

end module cernlib_module


subroutine initialize_cernlib()
  use cernlib_module, only: module_initialize_cernlib => initialize_cernlib
  implicit none

  call module_initialize_cernlib()
end subroutine initialize_cernlib


subroutine finalize_cernlib()
  use cernlib_module, only: module_finalize_cernlib => finalize_cernlib
  implicit none

  call module_finalize_cernlib()
end subroutine finalize_cernlib


subroutine abend()
  use cernlib_module, only: module_abend => abend
  implicit none

  call module_abend()
end subroutine abend


integer function lenocc(chv) result(last_occupied)
  use cernlib_module, only: module_lenocc => lenocc
  implicit none
  character(len=*), intent(in) :: chv

  last_occupied = module_lenocc(chv)
end function lenocc


subroutine dlsqp2(n, x, y, a0, a1, a2, sd, ifail)
  use cernlib_module, only: module_dlsqp2 => dlsqp2
  implicit none
  integer, intent(in) :: n
  double precision, intent(in) :: x(*)
  double precision, intent(in) :: y(*)
  double precision, intent(out) :: a0
  double precision, intent(out) :: a1
  double precision, intent(out) :: a2
  double precision, intent(out) :: sd
  integer, intent(out) :: ifail

  call module_dlsqp2(n, x, y, a0, a1, a2, sd, ifail)
end subroutine dlsqp2


double precision function dfint(narg, arg, nent, ent, table)
  use cernlib_module, only: module_dfint => dfint
  implicit none
  integer, intent(in) :: narg
  double precision, intent(in) :: arg(*)
  integer, intent(in) :: nent(*)
  double precision, intent(in) :: ent(*)
  double precision, intent(in) :: table(*)

  dfint = module_dfint(narg, arg, nent, ent, table)
end function dfint


subroutine mtlset(error_code, log_unit, message_limit, return_limit)
  use cernlib_module, only: module_mtlset => mtlset
  implicit none
  character(len=6), intent(in) :: error_code
  integer, intent(in) :: log_unit
  integer, intent(in) :: message_limit
  integer, intent(in) :: return_limit

  call module_mtlset(error_code, log_unit, message_limit, return_limit)
end subroutine mtlset


subroutine mtlmtr(error_code, log_unit, message_flag, return_flag)
  use cernlib_module, only: module_mtlmtr => mtlmtr
  implicit none
  character(len=6), intent(in) :: error_code
  integer, intent(out) :: log_unit
  logical, intent(out) :: message_flag
  logical, intent(out) :: return_flag

  call module_mtlmtr(error_code, log_unit, message_flag, return_flag)
end subroutine mtlmtr


subroutine kerset(error_code, log_unit, message_limit, return_limit)
  use cernlib_module, only: module_kerset => kerset
  implicit none
  character(len=6), intent(in) :: error_code
  integer, intent(in) :: log_unit
  integer, intent(in) :: message_limit
  integer, intent(in) :: return_limit

  call module_kerset(error_code, log_unit, message_limit, return_limit)
end subroutine kerset


subroutine kermtr(error_code, log_unit, message_flag, return_flag)
  use cernlib_module, only: module_kermtr => kermtr
  implicit none
  character(len=6), intent(in) :: error_code
  integer, intent(out) :: log_unit
  logical, intent(out) :: message_flag
  logical, intent(out) :: return_flag

  call module_kermtr(error_code, log_unit, message_flag, return_flag)
end subroutine kermtr


subroutine mtlprt(name, error_code, text)
  use cernlib_module, only: module_mtlprt => mtlprt
  implicit none
  character(len=*), intent(in) :: name
  character(len=*), intent(in) :: error_code
  character(len=*), intent(in) :: text

  call module_mtlprt(name, error_code, text)
end subroutine mtlprt


subroutine radmul(f, n, a, b, minpts, maxpts, eps, wk, iwk, result, &
                  relerr, nfnevl, ifail)
  use cernlib_module, only: module_radmul => radmul
  implicit none
  real, external :: f
  integer, intent(in) :: n
  real, intent(in) :: a(*)
  real, intent(in) :: b(*)
  integer, intent(in) :: minpts
  integer, intent(in) :: maxpts
  real, intent(in) :: eps
  real, intent(inout) :: wk(*)
  integer, intent(in) :: iwk
  real, intent(out) :: result
  real, intent(out) :: relerr
  integer, intent(out) :: nfnevl
  integer, intent(out) :: ifail

  call module_radmul(f, n, a, b, minpts, maxpts, eps, wk, iwk, &
                     result, relerr, nfnevl, ifail)
end subroutine radmul


subroutine dadmul(f, n, a, b, minpts, maxpts, eps, wk, iwk, result, &
                  relerr, nfnevl, ifail)
  use cernlib_module, only: module_dadmul => dadmul
  implicit none
  interface
    function f(n, z) result(value)
      integer, intent(in) :: n
      double precision, intent(in) :: z(*)
      double precision :: value
    end function f
  end interface
  integer, intent(in) :: n
  double precision, intent(in) :: a(*)
  double precision, intent(in) :: b(*)
  integer, intent(in) :: minpts
  integer, intent(in) :: maxpts
  double precision, intent(in) :: eps
  double precision, intent(inout) :: wk(*)
  integer, intent(in) :: iwk
  double precision, intent(out) :: result
  double precision, intent(out) :: relerr
  integer, intent(out) :: nfnevl
  integer, intent(out) :: ifail

  call module_dadmul(f, n, a, b, minpts, maxpts, eps, wk, iwk, &
                     result, relerr, nfnevl, ifail)
end subroutine dadmul
