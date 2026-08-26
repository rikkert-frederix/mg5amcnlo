module polynomial_fit
  implicit none
  private

  integer, parameter :: maxdeg = 20

  integer :: ndim = 0
  integer :: nchan = 0
  integer :: maxpoint = 0
  integer :: n_ord_virt = 0
  integer, allocatable :: n(:)
  integer, allocatable :: ndeg(:, :, :)
  double precision, allocatable :: r(:)
  double precision, allocatable :: yp(:)
  double precision, allocatable :: a(:)
  double precision, allocatable :: y(:, :, :)
  double precision, allocatable :: w(:, :, :)
  double precision, allocatable :: x2d(:, :, :)
  double precision, allocatable :: a2d(:, :, :, :)
  logical, allocatable :: valid_ord_virt(:)
  logical :: fit_done = .false.
  logical :: initialized = .false.

  public :: init_polyfit
  public :: add_point_polyfit
  public :: do_polyfit
  public :: get_polyfit
  public :: save_polyfit
  public :: restore_polyfit

contains

  subroutine init_polyfit(ndims, nchans, k_ord_virt, npoints)
    implicit none
    integer, intent(in) :: ndims
    integer, intent(in) :: nchans
    integer, intent(in) :: k_ord_virt
    integer, intent(in) :: npoints

    if (ndims < 0 .or. nchans <= 0 .or. k_ord_virt < 0 .or. &
        npoints < 0) then
      write (*, *) 'Invalid dimensions passed to init_polyfit'
      stop 1
    end if

    call finalize_polyfit()

    ndim = ndims
    nchan = nchans
    maxpoint = npoints
    n_ord_virt = k_ord_virt

    allocate (ndeg(ndim, 0:n_ord_virt, nchan))
    allocate (n(nchan))
    allocate (a(maxpoint*3 + maxdeg*3 + 3))
    allocate (y(maxpoint, 0:n_ord_virt, nchan))
    allocate (w(maxpoint, 0:n_ord_virt, nchan))
    allocate (r(maxpoint))
    allocate (yp(maxdeg))
    allocate (x2d(maxpoint, ndim, nchan))
    allocate (a2d(maxdeg*3 + 3, ndim, 0:n_ord_virt, nchan))
    allocate (valid_ord_virt(0:n_ord_virt))

    n = 0
    fit_done = .false.
    valid_ord_virt = .false.
    initialized = .true.
  end subroutine init_polyfit

  subroutine finalize_polyfit()
    implicit none

    if (allocated(ndeg)) deallocate (ndeg)
    if (allocated(n)) deallocate (n)
    if (allocated(a)) deallocate (a)
    if (allocated(y)) deallocate (y)
    if (allocated(w)) deallocate (w)
    if (allocated(r)) deallocate (r)
    if (allocated(yp)) deallocate (yp)
    if (allocated(x2d)) deallocate (x2d)
    if (allocated(a2d)) deallocate (a2d)
    if (allocated(valid_ord_virt)) deallocate (valid_ord_virt)

    ndim = 0
    nchan = 0
    maxpoint = 0
    n_ord_virt = 0
    fit_done = .false.
    initialized = .false.
  end subroutine finalize_polyfit

  subroutine add_point_polyfit(ichan, k_ord_virt, x, virt, absxs)
    implicit none
    integer, intent(in) :: ichan
    integer, intent(in) :: k_ord_virt
    double precision, intent(in) :: x(:)
    double precision, intent(in) :: virt
    double precision, intent(in) :: absxs
    integer :: i

    call require_initialized('add_point_polyfit')
    if (ichan < 1 .or. ichan > nchan) then
      write (*, *) 'Invalid channel passed to add_point_polyfit:', ichan
      stop 1
    end if
    if (k_ord_virt < 0 .or. k_ord_virt > n_ord_virt) then
      write (*, *) 'Invalid order passed to add_point_polyfit:', k_ord_virt
      stop 1
    end if
    if (size(x) < ndim) then
      write (*, *) 'Too few coordinates passed to add_point_polyfit'
      stop 1
    end if

    valid_ord_virt(k_ord_virt) = .true.
    if (k_ord_virt == 0) then
      if (n(ichan) == maxpoint) call grow_point_storage()
      n(ichan) = n(ichan) + 1
      do i = 1, ndim
        x2d(n(ichan), i, ichan) = x(i)
      end do
      y(n(ichan), 0:n_ord_virt, ichan) = 0.0d0
      w(n(ichan), 0:n_ord_virt, ichan) = 0.0d0
    end if

    y(n(ichan), k_ord_virt, ichan) = virt
    w(n(ichan), k_ord_virt, ichan) = sqrt(abs(absxs))
  end subroutine add_point_polyfit

  subroutine grow_point_storage()
    implicit none
    integer :: i
    integer :: new_maxpoint
    double precision, allocatable :: temp_y(:, :, :)
    double precision, allocatable :: temp_w(:, :, :)
    double precision, allocatable :: temp_x2d(:, :, :)

    new_maxpoint = maxpoint + 1000
    allocate (temp_y(new_maxpoint, 0:n_ord_virt, nchan))
    allocate (temp_w(new_maxpoint, 0:n_ord_virt, nchan))
    allocate (temp_x2d(new_maxpoint, ndim, nchan))
    temp_y = 0.0d0
    temp_w = 0.0d0
    temp_x2d = 0.0d0

    do i = 1, nchan
      if (n(i) > 0) then
        temp_y(1:n(i), 0:n_ord_virt, i) = &
          y(1:n(i), 0:n_ord_virt, i)
        temp_w(1:n(i), 0:n_ord_virt, i) = &
          w(1:n(i), 0:n_ord_virt, i)
        temp_x2d(1:n(i), 1:ndim, i) = x2d(1:n(i), 1:ndim, i)
      end if
    end do

    deallocate (y)
    allocate (y(new_maxpoint, 0:n_ord_virt, nchan))
    y = temp_y
    deallocate (temp_y)

    deallocate (w)
    allocate (w(new_maxpoint, 0:n_ord_virt, nchan))
    w = temp_w
    deallocate (temp_w)

    deallocate (x2d)
    allocate (x2d(new_maxpoint, ndim, nchan))
    x2d = temp_x2d
    deallocate (temp_x2d)

    deallocate (r)
    allocate (r(new_maxpoint))
    deallocate (a)
    allocate (a(new_maxpoint*3 + maxdeg*3 + 3))
    maxpoint = new_maxpoint
  end subroutine grow_point_storage

  subroutine do_polyfit()
    implicit none
    integer :: ic
    integer :: ko
    integer :: i
    integer :: ierr
    double precision :: eps

    call require_initialized('do_polyfit')
    do ic = 1, nchan
      do ko = 0, n_ord_virt
        if (.not. valid_ord_virt(ko)) cycle
        do i = 1, ndim
          eps = -1.0d0
          call polfit(n(ic), x2d(1, i, ic), y(1, ko, ic), &
                      w(1, ko, ic), maxdeg, ndeg(i, ko, ic), eps, r, ierr, a)
          a2d(1:maxdeg*3 + 3, i, ko, ic) = a(1:maxdeg*3 + 3)
        end do
      end do
    end do
    fit_done = .true.
  end subroutine do_polyfit

  subroutine get_polyfit(ichan, k_ord_virt, x, fun_at_x)
    implicit none
    integer, intent(in) :: ichan
    integer, intent(in) :: k_ord_virt
    double precision, intent(in) :: x(:)
    double precision, intent(out) :: fun_at_x
    integer :: i
    double precision :: yfit
    double precision :: ave_fun

    call require_initialized('get_polyfit')
    if (ichan < 1 .or. ichan > nchan) then
      write (*, *) 'Invalid channel passed to get_polyfit:', ichan
      stop 1
    end if
    if (k_ord_virt < 0 .or. k_ord_virt > n_ord_virt) then
      write (*, *) 'Invalid order passed to get_polyfit:', k_ord_virt
      stop 1
    end if
    if (size(x) < ndim) then
      write (*, *) 'Too few coordinates passed to get_polyfit'
      stop 1
    end if

    if ((.not. fit_done) .or. &
        (.not. valid_ord_virt(k_ord_virt))) then
      fun_at_x = 0.0d0
      return
    end if

    fun_at_x = 0.0d0
    do i = 1, ndim
      a(1:maxdeg*3 + 3) = a2d(1:maxdeg*3 + 3, i, k_ord_virt, ichan)
      call pvalue(ndeg(i, k_ord_virt, ichan), 0, x(i), yfit, yp, a)
      call pvalue(0, 0, x(i), ave_fun, yp, a)
      fun_at_x = fun_at_x + (yfit - ave_fun) + ave_fun/dble(ndim)
    end do
  end subroutine get_polyfit

  subroutine save_polyfit(iunit)
    implicit none
    integer, intent(in) :: iunit
    integer :: ic
    integer :: i

    call require_initialized('save_polyfit')
    do ic = 1, nchan
      write (iunit, *) 'POL', n(ic)
    end do
    do ic = 1, nchan
      do i = 1, n(ic)
        write (iunit, *) 'POL', x2d(i, 1:ndim, ic), &
          y(i, 0:n_ord_virt, ic), w(i, 0:n_ord_virt, ic)
      end do
    end do
  end subroutine save_polyfit

  subroutine restore_polyfit(iunit)
    implicit none
    integer, intent(in) :: iunit
    integer :: ic
    integer :: i
    integer :: ko
    character(len=3) :: dummy

    call require_initialized('restore_polyfit')
    do ic = 1, nchan
      read (iunit, *) dummy, n(ic)
      if (n(ic) < 0 .or. n(ic) > maxpoint) then
        write (*, *) 'Invalid point count in polyfit checkpoint:', n(ic)
        stop 1
      end if
    end do
    do ic = 1, nchan
      do i = 1, n(ic)
        read (iunit, *) dummy, x2d(i, 1:ndim, ic), &
          y(i, 0:n_ord_virt, ic), w(i, 0:n_ord_virt, ic)
      end do
    end do
    do ic = 1, nchan
      do ko = 0, n_ord_virt
        do i = 1, n(ic)
          if (y(i, ko, ic) /= 0.0d0) valid_ord_virt(ko) = .true.
        end do
      end do
    end do
  end subroutine restore_polyfit

  subroutine require_initialized(caller)
    implicit none
    character(len=*), intent(in) :: caller

    if (.not. initialized) then
      write (*, *) trim(caller), ': polynomial fitter is not initialized'
      stop 1
    end if
  end subroutine require_initialized

  ! Double-precision adaptation of the SLATEC POLFIT routine.
  subroutine polfit(nin, x, ydata, weights, max_degree, fit_degree, &
                    eps, results, ierr, work)
    implicit none
    integer, intent(in) :: nin
    integer, intent(in) :: max_degree
    integer, intent(out) :: fit_degree
    integer, intent(out) :: ierr
    double precision, intent(in) :: x(*)
    double precision, intent(in) :: ydata(*)
    double precision, intent(inout) :: weights(*)
    double precision, intent(inout) :: eps
    double precision, intent(out) :: results(*)
    double precision, intent(inout) :: work(*)
    integer :: m
    integer :: mop1
    integer :: i
    integer :: idegf
    integer :: ksig
    integer :: k1
    integer :: k2
    integer :: k3
    integer :: k4
    integer :: k5
    integer :: k4pi
    integer :: k5pi
    integer :: j
    integer :: jp1
    integer :: k1pj
    integer :: k2pj
    integer :: k3pi
    integer :: nfail
    integer :: jpas
    integer :: nder
    double precision :: temd1
    double precision :: temd2
    double precision :: co(4, 3)
    double precision :: dummy_yp(0)
    double precision :: xm
    double precision :: etst
    double precision :: w11
    double precision :: sigj
    double precision :: sigjm1
    double precision :: w1
    double precision :: temp
    double precision :: degf
    double precision :: den
    double precision :: fcrit
    double precision :: fstat
    double precision :: sigpas
    double precision :: sig
    save co
    data co/-13.086850, -2.4648165, -3.3846535, -1.2973162, &
      -3.3381146, -1.7812271, -3.2578406, -1.6589279, &
      -1.6282703, -1.3152745, -3.2640179, -1.9829776/

    m = abs(nin)
    if (m == 0) go to 30
    if (max_degree < 0) go to 30
    work(1) = max_degree
    mop1 = max_degree + 1
    if (m < mop1) go to 30
    if (eps < 0.0d0 .and. m == mop1) go to 30
    xm = m
    etst = eps*eps*xm
    if (weights(1) < 0.0d0) go to 2
    do i = 1, m
      if (weights(i) <= 0.0d0) go to 30
    end do
    go to 4

2   do i = 1, m
      weights(i) = 1.0d0
    end do

4   if (eps >= 0.0d0) go to 8
    if (eps > -0.55d0) go to 5
    idegf = m - max_degree - 1
    ksig = 1
    if (idegf < 10) ksig = 2
    if (idegf < 5) ksig = 3
    go to 8

5   ksig = 1
    if (eps < -0.03d0) ksig = 2
    if (eps < -0.07d0) ksig = 3

8   k1 = max_degree + 1
    k2 = k1 + max_degree
    k3 = k2 + max_degree + 2
    k4 = k3 + m
    k5 = k4 + m
    do i = 2, k4
      work(i) = 0.0d0
    end do
    w11 = 0.0d0
    if (nin < 0) go to 11

    do i = 1, m
      k4pi = k4 + i
      work(k4pi) = 1.0d0
      w11 = w11 + weights(i)
    end do
    go to 13

11  do i = 1, m
      k4pi = k4 + i
      w11 = w11 + weights(i)*work(k4pi)**2
    end do

13  temd1 = 0.0d0
    do i = 1, m
      k4pi = k4 + i
      temd1 = temd1 + weights(i)*ydata(i)*work(k4pi)
    end do
    temd1 = temd1/w11
    work(k2 + 1) = temd1
    sigj = 0.0d0
    do i = 1, m
      k4pi = k4 + i
      k5pi = k5 + i
      temd2 = temd1*work(k4pi)
      results(i) = temd2
      work(k5pi) = temd2 - results(i)
      sigj = sigj + weights(i)*((ydata(i) - results(i)) - work(k5pi))**2
    end do
    j = 0

    if (eps < 0.0d0) then
      go to 24
    else if (eps == 0.0d0) then
      go to 26
    else
      go to 27
    end if

16  j = j + 1
    jp1 = j + 1
    k1pj = k1 + j
    k2pj = k2 + j
    sigjm1 = sigj
    if (j > 1) work(k1pj) = w11/w1

    temd1 = 0.0d0
    do i = 1, m
      k4pi = k4 + i
      temd2 = work(k4pi)
      temd1 = temd1 + x(i)*weights(i)*temd2*temd2
    end do
    work(jp1) = temd1/w11

    w1 = w11
    w11 = 0.0d0
    do i = 1, m
      k3pi = k3 + i
      k4pi = k4 + i
      temp = work(k3pi)
      work(k3pi) = work(k4pi)
      work(k4pi) = (x(i) - work(jp1))*work(k3pi) - work(k1pj)*temp
      w11 = w11 + weights(i)*work(k4pi)**2
    end do

    temd1 = 0.0d0
    do i = 1, m
      k4pi = k4 + i
      k5pi = k5 + i
      temd2 = weights(i)*((ydata(i) - results(i)) - work(k5pi))* &
              work(k4pi)
      temd1 = temd1 + temd2
    end do
    temd1 = temd1/w11
    work(k2pj + 1) = temd1

    sigj = 0.0d0
    do i = 1, m
      k4pi = k4 + i
      k5pi = k5 + i
      temd2 = results(i) + work(k5pi) + temd1*work(k4pi)
      results(i) = temd2
      work(k5pi) = temd2 - results(i)
      sigj = sigj + weights(i)*((ydata(i) - results(i)) - work(k5pi))**2
    end do

    if (eps < 0.0d0) then
      go to 23
    else if (eps == 0.0d0) then
      go to 26
    else
      go to 27
    end if

23  if (sigj == 0.0d0) go to 29
    degf = m - j - 1
    den = (co(4, ksig)*degf + 1.0d0)*degf
    fcrit = (((co(3, ksig)*degf) + co(2, ksig))*degf + co(1, ksig))/den
    fcrit = fcrit*fcrit
    fstat = (sigjm1 - sigj)*degf/sigj
    if (fstat < fcrit) go to 25

24  sigpas = sigj
    jpas = j
    nfail = 0
    if (max_degree == j) go to 32
    go to 16

25  nfail = nfail + 1
    if (nfail >= 3) go to 29
    if (max_degree == j) go to 32
    go to 16

26  if (max_degree == j) go to 28
    go to 16

27  if (sigj <= etst) go to 28
    if (max_degree == j) go to 31
    go to 16

28  ierr = 1
    fit_degree = j
    sig = sigj
    go to 33

29  ierr = 1
    fit_degree = jpas
    sig = sigpas
    go to 33

30  ierr = 2
    write (*, *) 'SLATEC POLFIT INVALID INPUT PARAMETER.', 2, 1
    stop 1

31  ierr = 3
    fit_degree = max_degree
    sig = sigj
    go to 33

32  ierr = 4
    fit_degree = jpas
    sig = sigpas

33  work(k3) = fit_degree
    if (eps >= 0.0d0 .or. fit_degree == max_degree) go to 36
    nder = 0
    do i = 1, m
      call pvalue(fit_degree, nder, x(i), results(i), dummy_yp, work)
    end do

36  eps = sqrt(sig/xm)
  end subroutine polfit

  ! Double-precision adaptation of the SLATEC PVALUE routine.
  subroutine pvalue(l, nder, x, yfit, derivatives, work)
    implicit none
    integer, intent(in) :: l
    integer, intent(in) :: nder
    double precision, intent(in) :: x
    double precision, intent(out) :: yfit
    double precision, intent(inout) :: derivatives(*)
    double precision, intent(inout) :: work(*)
    integer :: ndo
    integer :: maxord
    integer :: k1
    integer :: k2
    integer :: k3
    integer :: nord
    integer :: k4
    integer :: i
    integer :: ndp1
    integer :: k3p1
    integer :: k4p1
    integer :: lp1
    integer :: lm1
    integer :: ilo
    integer :: iup
    integer :: kc
    integer :: indx
    integer :: inp1
    integer :: k1i
    integer :: ic
    integer :: deriv
    integer :: k3pn
    integer :: k4pn
    double precision :: val
    double precision :: cc
    double precision :: dif
    character(len=8) :: xern1
    character(len=8) :: xern2

    if (l < 0) go to 12
    ndo = max(nder, 0)
    ndo = min(ndo, l)
    maxord = nint(work(1))
    k1 = maxord + 1
    k2 = k1 + maxord
    k3 = k2 + maxord + 2
    nord = nint(work(k3))
    if (l > nord) go to 11
    k4 = k3 + l + 1
    if (nder < 1) go to 2
    do i = 1, nder
      derivatives(i) = 0.0d0
    end do

2   if (l >= 2) go to 4
    if (l == 1) go to 3
    val = work(k2 + 1)
    go to 10

3   cc = work(k2 + 2)
    val = work(k2 + 1) + (x - work(2))*cc
    if (nder >= 1) derivatives(1) = cc
    go to 10

4   ndp1 = ndo + 1
    k3p1 = k3 + 1
    k4p1 = k4 + 1
    lp1 = l + 1
    lm1 = l - 1
    ilo = k3 + 3
    iup = k4 + ndp1
    do i = ilo, iup
      work(i) = 0.0d0
    end do
    dif = x - work(lp1)
    kc = k2 + lp1
    work(k4p1) = work(kc)
    work(k3p1) = work(kc - 1) + dif*work(k4p1)
    work(k3 + 2) = work(k4p1)

    do i = 1, lm1
      indx = l - i
      inp1 = indx + 1
      k1i = k1 + inp1
      ic = k2 + indx
      dif = x - work(inp1)
      val = work(ic) + dif*work(k3p1) - work(k1i)*work(k4p1)
      if (ndo <= 0) go to 8
      do deriv = 1, ndo
        k3pn = k3p1 + deriv
        k4pn = k4p1 + deriv
        derivatives(deriv) = dif*work(k3pn) + &
                             deriv*work(k3pn - 1) - work(k1i)*work(k4pn)
      end do
      do deriv = 1, ndo
        k3pn = k3p1 + deriv
        k4pn = k4p1 + deriv
        work(k4pn) = work(k3pn)
        work(k3pn) = derivatives(deriv)
      end do
8     work(k4p1) = work(k3p1)
      work(k3p1) = val
    end do

10  yfit = val
    return

11  write (xern1, '(i8)') l
    write (xern2, '(i8)') nord
    write (*, *) 'SLATEC PVALUE '// &
      'THE ORDER OF POLYNOMIAL EVALUATION, L = '//xern1// &
      ' REQUESTED EXCEEDS THE HIGHEST ORDER FIT, NORD = '//xern2// &
      ', COMPUTED BY POLFIT -- EXECUTION TERMINATED.', 8, 2
    stop 1

12  write (*, *) 'SLATEC PVALUE '// &
      'INVALID INPUT PARAMETER. ORDER OF POLYNOMIAL EVALUATION '// &
      'REQUESTED IS NEGATIVE -- EXECUTION TERMINATED.', 2, 2
    stop 1
  end subroutine pvalue

end module polynomial_fit
