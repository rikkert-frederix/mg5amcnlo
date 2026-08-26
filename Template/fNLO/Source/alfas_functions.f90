module alfas_functions_module
  use fnlo_runtime_common, only: cmass, bmass, asmz, nloop
  implicit none
  private

  integer, parameter :: dp = kind(1d0)
  integer, parameter :: gridsize_low = 255
  integer, parameter :: gridsize_mid = 255
  integer, parameter :: gridsize_high = 255
  real(dp), parameter :: cut_off_low = -0.69314718055994530942d0
  real(dp), parameter :: cut_off_high = 11.512925464970228420d0
  real(dp), parameter :: zmass = 91.188d0

  real(dp), parameter :: b0(3:5) = (/ &
       0.716197243913527d0, 0.66314559621623d0, &
       0.61009394851893d0 /)
  real(dp), parameter :: c1(3:5) = (/ &
       0.565884242104515d0, 0.49019722472304d0, &
       0.40134724779695d0 /)
  real(dp), parameter :: c2(3:5) = (/ &
       0.453013579178645d0, 0.30879037953664d0, &
       0.14942733137107d0 /)
  real(dp), parameter :: del(3:5) = (/ &
       1.22140465909230d0, 0.99743079911360d0, &
       0.66077962451190d0 /)

  logical :: qmass_initialized = .false.

  real(dp) :: cached_asmz = 0d0
  integer :: cached_nloop = 0
  real(dp) :: cached_amb, cached_amc

  logical :: first_grid_fill = .true.
  real(dp) :: log_cmass, log_bmass
  real(dp) :: grid_low(0:gridsize_low)
  real(dp) :: grid_mid(0:gridsize_mid)
  real(dp) :: grid_high(0:gridsize_high)
  real(dp) :: grid_low_q(0:gridsize_low)
  real(dp) :: grid_mid_q(0:gridsize_mid)
  real(dp) :: grid_high_q(0:gridsize_high)
  real(dp) :: grid_low_sf(0:gridsize_low)
  real(dp) :: grid_mid_sf(0:gridsize_mid)
  real(dp) :: grid_high_sf(0:gridsize_high)

  public :: alphas

contains

  double precision function alphas(q)
    implicit none
    real(dp), intent(in) :: q

    alphas = alphas_from_grids_impl(q, asmz, nloop)
  end function alphas

  subroutine initialize_alfas_backend()
    implicit none

    if (qmass_initialized) return
    cmass = 1.42d0
    bmass = 4.7d0
    qmass_initialized = .true.
  end subroutine initialize_alfas_backend

  double precision function alphas_from_grids_impl(scale, asmz, nloop)
    implicit none
    real(dp), intent(in) :: scale, asmz
    integer, intent(in) :: nloop
    real(dp) :: lscale, q
    integer :: i, ii

    if (first_grid_fill) then
      write (*,*) 'Filling grid for alpha_S computation ... '
      log_cmass = log(1.42d0)
      log_bmass = log(4.7d0)

      i = 0
      q = cut_off_low
      grid_low_q(i) = q
      grid_low(i) = alphas_not_timed_impl(exp(q), asmz, nloop)
      do i = 1, gridsize_low
        q = cut_off_low + (log_cmass - cut_off_low) * &
             dble(i) / dble(gridsize_low)
        grid_low_q(i) = q
        grid_low(i) = alphas_not_timed_impl(exp(q), asmz, nloop)
        grid_low_sf(i - 1) = (grid_low(i) - grid_low(i - 1)) / &
             (grid_low_q(i) - grid_low_q(i - 1))
      end do

      i = 0
      q = log_cmass
      grid_mid_q(i) = q
      grid_mid(i) = alphas_not_timed_impl(exp(q), asmz, nloop)
      do i = 1, gridsize_mid
        q = log_cmass + (log_bmass - log_cmass) * &
             dble(i) / dble(gridsize_mid)
        grid_mid_q(i) = q
        grid_mid(i) = alphas_not_timed_impl(exp(q), asmz, nloop)
        grid_mid_sf(i - 1) = (grid_mid(i) - grid_mid(i - 1)) / &
             (grid_mid_q(i) - grid_mid_q(i - 1))
      end do

      i = 0
      q = log_bmass
      grid_high_q(i) = q
      grid_high(i) = alphas_not_timed_impl(exp(q), asmz, nloop)
      do i = 1, gridsize_high
        q = log_bmass + (cut_off_high - log_bmass) * &
             dble(i) / dble(gridsize_high)
        grid_high_q(i) = q
        grid_high(i) = alphas_not_timed_impl(exp(q), asmz, nloop)
        grid_high_sf(i - 1) = (grid_high(i) - grid_high(i - 1)) / &
             (grid_high_q(i) - grid_high_q(i - 1))
      end do

      first_grid_fill = .false.
      write (*,*) 'done... '
    end if

    lscale = log(scale)
    if (lscale < cut_off_low) then
      alphas_from_grids_impl = &
           alphas_not_timed_impl(scale, asmz, nloop)
    else if (lscale < log_cmass) then
      ii = int((lscale - cut_off_low) / &
           (log_cmass - cut_off_low) * gridsize_low)
      alphas_from_grids_impl = grid_low(ii) + &
           (lscale - grid_low_q(ii)) * grid_low_sf(ii)
    else if (lscale < log_bmass) then
      ii = int((lscale - log_cmass) / &
           (log_bmass - log_cmass) * gridsize_mid)
      alphas_from_grids_impl = grid_mid(ii) + &
           (lscale - grid_mid_q(ii)) * grid_mid_sf(ii)
    else if (lscale < cut_off_high) then
      ii = int((lscale - log_bmass) / &
           (cut_off_high - log_bmass) * gridsize_high)
      alphas_from_grids_impl = grid_high(ii) + &
           (lscale - grid_high_q(ii)) * grid_high_sf(ii)
    else
      alphas_from_grids_impl = &
           alphas_not_timed_impl(scale, asmz, nloop)
    end if
  end function alphas_from_grids_impl


  double precision function alphas_not_timed_impl(q, asmz, nloop)
    implicit none
    real(dp), intent(in) :: q, asmz
    integer, intent(in) :: nloop
    real(dp) :: as_out, t
    integer, parameter :: nf3 = 3, nf4 = 4, nf5 = 5

    call initialize_alfas_backend()

    if (q <= 0d0) then
      write (6,*) 'q .le. 0 in alphas'
      write (6,*) 'q= ', q
      stop
    end if
    if (asmz <= 0d0) then
      write (6,*) 'asmz .le. 0 in alphas', asmz
      stop
    end if
    if (cmass <= 0.3d0) then
      write (6,*) 'cmass .le. 0.3GeV in alphas', cmass
      stop
    end if
    if (bmass <= 0d0) then
      write (6,*) 'bmass .le. 0 in alphas', bmass
      write (6,*) 'COMMON/QMASS/CMASS,BMASS'
      stop
    end if

    if ((asmz /= cached_asmz) .or. (nloop /= cached_nloop)) then
      cached_asmz = asmz
      cached_nloop = nloop
      t = 2d0 * log(bmass / zmass)
      call newton1_impl(t, asmz, cached_amb, nloop, nf5)
      t = 2d0 * log(cmass / bmass)
      call newton1_impl(t, cached_amb, cached_amc, nloop, nf4)
    end if

    if (q < bmass) then
      if (q < cmass) then
        t = 2d0 * log(q / cmass)
        call newton1_impl(t, cached_amc, as_out, nloop, nf3)
      else
        t = 2d0 * log(q / bmass)
        call newton1_impl(t, cached_amb, as_out, nloop, nf4)
      end if
    else
      t = 2d0 * log(q / zmass)
      call newton1_impl(t, asmz, as_out, nloop, nf5)
    end if
    alphas_not_timed_impl = as_out
  end function alphas_not_timed_impl


  subroutine newton1_impl(t, a_in, a_out, nloop, nf)
    implicit none
    real(dp), intent(in) :: t, a_in
    real(dp), intent(out) :: a_out
    integer, intent(in) :: nloop, nf
    real(dp), parameter :: tol = 5d-4
    real(dp) :: as, delta, f, fp

    a_out = a_in / (1d0 + a_in * b0(nf) * t)
    if (nloop == 1) return

    a_out = a_in / (1d0 + b0(nf) * a_in * t + &
         c1(nf) * a_in * log(1d0 + a_in * b0(nf) * t))
    if (a_out < 0d0) as = 0.3d0

    do
      as = a_out
      if (nloop == 2) then
        f = b0(nf) * t + f2_value(a_in, nf) - f2_value(as, nf)
        fp = 1d0 / (as**2 * (1d0 + c1(nf) * as))
      end if
      if (nloop == 3) then
        f = b0(nf) * t + f3_value(a_in, nf) - f3_value(as, nf)
        fp = 1d0 / &
             (as**2 * (1d0 + c1(nf) * as + c2(nf) * as**2))
      end if
      a_out = as - f / fp
      delta = abs(f / fp / as)
      if (delta > tol) cycle
      exit
    end do
  end subroutine newton1_impl
  pure double precision function f2_value(as, nf)
    implicit none
    real(dp), intent(in) :: as
    integer, intent(in) :: nf

    f2_value = 1d0 / as + &
         c1(nf) * log((c1(nf) * as) / (1d0 + c1(nf) * as))
  end function f2_value


  pure double precision function f3_value(as, nf)
    implicit none
    real(dp), intent(in) :: as
    integer, intent(in) :: nf

    f3_value = 1d0 / as + 0.5d0 * c1(nf) * &
         log((c2(nf) * as**2) / &
         (1d0 + c1(nf) * as + c2(nf) * as**2)) - &
         (c1(nf)**2 - 2d0 * c2(nf)) / del(nf) * &
         atan((2d0 * c2(nf) * as + c1(nf)) / del(nf))
  end function f3_value

end module alfas_functions_module
