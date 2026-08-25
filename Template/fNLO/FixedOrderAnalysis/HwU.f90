! Histograms with uncertainties (HwU).
!
! The histogram implementation is kept in modules.  External entry points and
! process-generated PineAPPL declarations are isolated in HwU_bridge.f.

module HwU_module
  implicit none
  private

  integer, parameter :: wgts_info_len = 80
  integer :: max_plots = 0
  integer :: max_points = 0
  integer :: max_bins = 0
  integer :: nwgts = 0
  integer :: np = 0
  integer :: error_estimation = 3

  logical, allocatable :: booked(:)
  integer, allocatable :: nbin(:)
  integer, allocatable :: histi(:, :)
  integer, allocatable :: p_bin(:)
  integer, allocatable :: p_label(:)
  character(len=50), allocatable :: title(:)
  character(len=wgts_info_len), allocatable :: wgts_info(:)
  double precision, allocatable :: histy(:, :, :)
  double precision, allocatable :: histy_acc(:, :, :)
  double precision, allocatable :: histy2(:, :)
  double precision, allocatable :: histy_err(:, :)
  double precision, allocatable :: histxl(:, :)
  double precision, allocatable :: histxm(:, :)
  double precision, allocatable :: step(:)
  double precision, allocatable :: p_wgts(:, :)

  double precision :: accumulated_iterations = 0d0

  interface
    subroutine hwu_pineappl_inithist()
      implicit none
    end subroutine hwu_pineappl_inithist

    subroutine hwu_pineappl_book(label, title_l, nbin_l, xmin, xmax)
      implicit none
      integer, intent(in) :: label, nbin_l
      character(len=*), intent(in) :: title_l
      double precision, intent(in) :: xmin, xmax
    end subroutine hwu_pineappl_book

    subroutine hwu_pineappl_fill(label, x)
      implicit none
      integer, intent(in) :: label
      double precision, intent(in) :: x
    end subroutine hwu_pineappl_fill
  end interface

  public :: HwU_inithist
  public :: HwU_book
  public :: HwU_fill
  public :: HwU_add_points
  public :: HwU_accum_iter
  public :: HwU_output
  public :: accum
  public :: addfil

contains

  subroutine HwU_inithist(nweights, wgt_info)
    integer, intent(in) :: nweights
    character(len=*), intent(in) :: wgt_info(*)
    integer :: i

    call hwu_pineappl_inithist()
    call HwU_deallocate_all()

    max_plots = 0
    max_points = 0
    max_bins = 0
    np = 0
    nwgts = nweights
    accumulated_iterations = 0d0

    allocate(wgts_info(nwgts))
    do i = 1, nwgts
      wgts_info(i) = wgt_info(i)
    end do
  end subroutine HwU_inithist


  subroutine set_error_estimation(input)
    integer, intent(in) :: input

    if (input >= 0 .and. input <= 3) then
      error_estimation = input
    else
      write (*, *) 'unknown error estimation', input
      stop 1
    end if
  end subroutine set_error_estimation


  subroutine HwU_book(label, title_l, nbin_l, xmin, xmax)
    integer, intent(in) :: label
    character(len=*), intent(in) :: title_l
    integer, intent(in) :: nbin_l
    double precision, intent(in) :: xmin
    double precision, intent(in) :: xmax
    integer :: i, j

    call hwu_pineappl_book(label, title_l, nbin_l, xmin, xmax)
    call HwU_allocate_histo(label, nbin_l)

    booked(label) = .true.
    title(label) = title_l
    nbin(label) = nbin_l
    step(label) = (xmax - xmin) / dble(nbin(label))

    do i = 1, nbin(label)
      histxl(label, i) = xmin + step(label) * dble(i - 1)
      histxm(label, i) = xmin + step(label) * dble(i)
      do j = 1, nwgts
        histy(j, label, i) = 0d0
        histy_acc(j, label, i) = 0d0
      end do
      histi(label, i) = 0
      histy2(label, i) = 0d0
      histy_err(label, i) = 0d0
    end do
  end subroutine HwU_book


  subroutine HwU_fill(label, x, wgts)
    integer, intent(in) :: label
    double precision, intent(in) :: x
    double precision, intent(in) :: wgts(*)
    integer :: i, j, bin

    call hwu_pineappl_fill(label, x)

    if (wgts(1) == 0d0) return
    if (x < histxl(label, 1) .or. x > histxm(label, nbin(label))) return

    bin = int((x - histxl(label, 1)) / step(label)) + 1
    if (bin < 1 .or. bin > nbin(label)) return

    do i = 1, np
      if (p_label(i) == label .and. p_bin(i) == bin) then
        do j = 1, nwgts
          p_wgts(j, i) = p_wgts(j, i) + wgts(j)
        end do
        return
      end if
    end do

    np = np + 1
    call HwU_allocate_p()
    p_label(np) = label
    p_bin(np) = bin
    do j = 1, nwgts
      p_wgts(j, np) = wgts(j)
    end do
  end subroutine HwU_fill


  subroutine HwU_add_points()
    integer :: i, j

    do i = 1, np
      do j = 1, nwgts
        histy(j, p_label(i), p_bin(i)) = &
          histy(j, p_label(i), p_bin(i)) + p_wgts(j, i)
      end do
      histi(p_label(i), p_bin(i)) = histi(p_label(i), p_bin(i)) + 1
      histy2(p_label(i), p_bin(i)) = &
        histy2(p_label(i), p_bin(i)) + p_wgts(1, i)**2
    end do
    np = 0
  end subroutine HwU_add_points


  subroutine HwU_accum_iter(inclde, nPSpoints, values)
    logical, intent(in) :: inclde
    integer, intent(in) :: nPSpoints
    double precision, intent(in) :: values(2)
    integer :: label, i, j
    double precision :: nPSinv

    nPSinv = 1d0 / dble(nPSpoints)
    if (inclde) accumulated_iterations = accumulated_iterations + 1d0

    do label = 1, max_plots
      if (.not. booked(label)) cycle
      if (inclde) call accumulate_results(label, nPSinv, &
           accumulated_iterations, values)

      do i = 1, nbin(label)
        do j = 1, nwgts
          histy(j, label, i) = 0d0
        end do
        histy2(label, i) = 0d0
        histi(label, i) = 0
      end do
    end do
  end subroutine HwU_accum_iter


  subroutine finalize_histograms(nPSpoints)
    integer, intent(in) :: nPSpoints
    integer :: label, i, j
    double precision :: nPSinv, niter, dummy(2)

    nPSinv = 1d0 / dble(nPSpoints)
    niter = 1d0
    do label = 1, max_plots
      if (.not. booked(label)) cycle
      do i = 1, nbin(label)
        do j = 1, nwgts
          histy_acc(j, label, i) = 0d0
        end do
        histy_err(label, i) = 0d0
      end do
      call accumulate_results(label, nPSinv, niter, dummy)
    end do
  end subroutine finalize_histograms


  subroutine accumulate_results(label, nPSinv, niter, values)
    integer, intent(in) :: label
    double precision, intent(in) :: nPSinv
    double precision, intent(in) :: niter
    double precision, intent(in) :: values(2)
    integer :: i, j
    double precision :: etot, y_squared, a1, a2
    double precision, allocatable :: vtot(:)

    allocate(vtot(nwgts))

    if (error_estimation == 2) then
      do i = 1, nbin(label)
        if (histi(label, i) == 0) cycle
        do j = 1, nwgts
          vtot(j) = histy(j, label, i) * nPSinv
        end do

        etot = sqrt(abs(histy2(label, i) * nPSinv - vtot(1)**2) * nPSinv)
        if (histi(label, i) > 1) then
          etot = etot * sqrt(dble(histi(label, i)) / &
            (dble(histi(label, i)) - 1.5d0))
        else
          etot = abs(vtot(1)) * 10d0
        end if

        if (histy_err(label, i) == 0d0) then
          do j = 1, nwgts
            histy_acc(j, label, i) = vtot(j)
          end do
          histy_err(label, i) = etot
        else
          do j = 1, nwgts
            histy_acc(j, label, i) = &
              (histy_acc(j, label, i) / histy_err(label, i) + &
               vtot(j) / etot) / &
              (1d0 / histy_err(label, i) + 1d0 / etot)
          end do
          histy_err(label, i) = 1d0 / sqrt(1d0 / histy_err(label, i)**2 + &
            1d0 / etot**2)
        end if
      end do

    else if (error_estimation == 3) then
      do i = 1, nbin(label)
        if (histi(label, i) == 0) cycle
        do j = 1, nwgts
          vtot(j) = histy(j, label, i) * nPSinv
        end do

        etot = sqrt(abs(histy2(label, i) * nPSinv - vtot(1)**2) * nPSinv)
        if (histi(label, i) > 1) then
          etot = etot * sqrt(dble(histi(label, i)) / &
            (dble(histi(label, i)) - 1.5d0))
        else
          etot = abs(vtot(1)) * 10d0
        end if

        if (histy_err(label, i) == 0d0) then
          do j = 1, nwgts
            histy_acc(j, label, i) = vtot(j)
          end do
          histy_err(label, i) = etot
        else
          do j = 1, nwgts
            histy_acc(j, label, i) = &
              (histy_acc(j, label, i) / values(2) + vtot(j) / values(1)) / &
              (1d0 / values(2) + 1d0 / values(1))
          end do
          a1 = ((1d0 / values(1)) / &
            ((1d0 / values(1)) + 1d0 / values(2)))**2
          a2 = ((1d0 / values(2)) / &
            ((1d0 / values(1)) + 1d0 / values(2)))**2
          histy_err(label, i) = sqrt(a2 * histy_err(label, i)**2 + &
            a1 * etot**2)
        end if
      end do

    else if (error_estimation == 1) then
      do i = 1, nbin(label)
        if (histi(label, i) == 0 .and. histy_acc(1, label, i) == 0d0) cycle
        if (niter /= 1d0) then
          y_squared = ((niter - 1d0) * histy_err(label, i))**2 + &
            (niter - 1d0) * histy_acc(1, label, i)**2
        end if
        do j = 1, nwgts
          vtot(j) = histy(j, label, i) * nPSinv
          histy_acc(j, label, i) = &
            (histy_acc(j, label, i) * (niter - 1d0) + vtot(j)) / niter
        end do
        if (niter == 1d0) then
          histy_err(label, i) = 0d0
        else
          histy_err(label, i) = sqrt(((y_squared + vtot(1)**2) / niter - &
            histy_acc(1, label, i)**2) / niter)
        end if
      end do

    else if (error_estimation == 0) then
      do i = 1, nbin(label)
        if (histi(label, i) == 0 .and. histy_acc(1, label, i) == 0d0) cycle
        do j = 1, nwgts
          vtot(j) = histy(j, label, i) * nPSinv
          histy_acc(j, label, i) = &
            (histy_acc(j, label, i) * (niter - 1d0) + vtot(j)) / niter
        end do
        if (histi(label, i) /= 0) then
          etot = sqrt(histy2(label, i)) * nPSinv
          histy_err(label, i) = sqrt(((niter - 1d0) * histy_err(label, i))**2 + &
            etot**2) / niter
        else
          histy_err(label, i) = (niter - 1d0) / niter * histy_err(label, i)
        end if
      end do
    end if

    deallocate(vtot)
  end subroutine accumulate_results


  subroutine HwU_output(unit, xnorm)
    integer, intent(in) :: unit
    double precision, intent(in) :: xnorm
    integer :: i, j, label
    character(len=4) :: str_nbin

    write (unit, '(a)', advance='no') '##& xmin'
    write (unit, '(a)', advance='no') ' & xmax'
    write (unit, '(a)', advance='no') ' & ' // trim(adjustl(wgts_info(1)))
    write (unit, '(a)', advance='no') ' & dy'
    do j = 2, nwgts
      write (unit, '(a)', advance='no') ' & ' // trim(adjustl(wgts_info(j)))
    end do
    write (unit, '(a)') ''
    write (unit, '(a)') ''

    do label = 1, max_plots
      if (.not. booked(label)) cycle
      write (str_nbin, '(i4.4)') nbin(label)
      write (unit, '(12a,4a,3a,a,2a)') &
        '<histogram> ', str_nbin, ' " ', trim(title(label)), ' "'
      do i = 1, nbin(label)
        write (unit, '(2x,e14.7)', advance='no') histxl(label, i)
        write (unit, '(2x,e14.7)', advance='no') histxm(label, i)
        write (unit, '(2x,e14.7)', advance='no') &
          histy_acc(1, label, i) * xnorm
        write (unit, '(2x,e14.7)', advance='no') histy_err(label, i) * xnorm
        do j = 2, nwgts
          write (unit, '(2x,e14.7)', advance='no') &
            histy_acc(j, label, i) * xnorm
        end do
        write (unit, *) ''
      end do
      write (unit, '(12a)') '<\histogram>'
      write (unit, '(a)') ''
      write (unit, '(a)') ''
    end do
  end subroutine HwU_output


  subroutine HwU_deallocate_all()
    if (allocated(wgts_info)) deallocate(wgts_info)
    if (allocated(booked)) deallocate(booked)
    if (allocated(title)) deallocate(title)
    if (allocated(nbin)) deallocate(nbin)
    if (allocated(step)) deallocate(step)
    if (allocated(histxl)) deallocate(histxl)
    if (allocated(histxm)) deallocate(histxm)
    if (allocated(histy)) deallocate(histy)
    if (allocated(histy_acc)) deallocate(histy_acc)
    if (allocated(histi)) deallocate(histi)
    if (allocated(histy2)) deallocate(histy2)
    if (allocated(histy_err)) deallocate(histy_err)
    if (allocated(p_bin)) deallocate(p_bin)
    if (allocated(p_label)) deallocate(p_label)
    if (allocated(p_wgts)) deallocate(p_wgts)
  end subroutine HwU_deallocate_all


  subroutine HwU_allocate_p()
    integer :: new_size
    integer, allocatable :: itemp(:)
    double precision, allocatable :: rtemp(:, :)

    if (.not. allocated(p_bin)) then
      max_points = max(max_plots, 1)
      allocate(p_bin(max_points), p_label(max_points))
      allocate(p_wgts(nwgts, max_points))
      p_bin = 0
      p_label = 0
      p_wgts = 0d0
    else if (np > max_points) then
      new_size = np + max(max_plots, 1)

      allocate(itemp(new_size))
      itemp = 0
      itemp(1:max_points) = p_bin
      deallocate(p_bin)
      allocate(p_bin(new_size))
      p_bin = itemp

      itemp = 0
      itemp(1:max_points) = p_label
      deallocate(p_label)
      allocate(p_label(new_size))
      p_label = itemp
      deallocate(itemp)

      allocate(rtemp(nwgts, new_size))
      rtemp = 0d0
      rtemp(:, 1:max_points) = p_wgts
      deallocate(p_wgts)
      allocate(p_wgts(nwgts, new_size))
      p_wgts = rtemp
      deallocate(rtemp)

      max_points = new_size
    end if
  end subroutine HwU_allocate_p


  subroutine HwU_allocate_histo(label, nbin_l)
    integer, intent(in) :: label
    integer, intent(in) :: nbin_l
    integer :: label_max, nbin_max
    logical, allocatable :: ltemp(:)
    integer, allocatable :: itemp1(:), itemp2(:, :)
    character(len=50), allocatable :: ctemp(:)
    double precision, allocatable :: rtemp1(:), rtemp2(:, :), rtemp3(:, :, :)

    if (.not. allocated(booked)) then
      max_plots = max(label, 1)
      max_bins = nbin_l
      allocate(booked(max_plots), title(max_plots), nbin(max_plots))
      allocate(step(max_plots))
      allocate(histxl(max_plots, max_bins), histxm(max_plots, max_bins))
      allocate(histy(nwgts, max_plots, max_bins))
      allocate(histy_acc(nwgts, max_plots, max_bins))
      allocate(histi(max_plots, max_bins))
      allocate(histy2(max_plots, max_bins), histy_err(max_plots, max_bins))

      booked = .false.
      title = ''
      nbin = 0
      step = 0d0
      histxl = 0d0
      histxm = 0d0
      histy = 0d0
      histy_acc = 0d0
      histi = 0
      histy2 = 0d0
      histy_err = 0d0
      return
    end if

    if (label > max_plots .or. nbin_l > max_bins) then
      label_max = max(label, max_plots)
      nbin_max = max(nbin_l, max_bins)

      allocate(ltemp(label_max))
      ltemp = .false.
      ltemp(1:max_plots) = booked
      deallocate(booked)
      allocate(booked(label_max))
      booked = ltemp
      deallocate(ltemp)

      allocate(ctemp(label_max))
      ctemp = ''
      ctemp(1:max_plots) = title
      deallocate(title)
      allocate(title(label_max))
      title = ctemp
      deallocate(ctemp)

      allocate(itemp1(label_max))
      itemp1 = 0
      itemp1(1:max_plots) = nbin
      deallocate(nbin)
      allocate(nbin(label_max))
      nbin = itemp1
      deallocate(itemp1)

      allocate(rtemp1(label_max))
      rtemp1 = 0d0
      rtemp1(1:max_plots) = step
      deallocate(step)
      allocate(step(label_max))
      step = rtemp1
      deallocate(rtemp1)

      allocate(rtemp2(label_max, nbin_max))
      rtemp2 = 0d0
      rtemp2(1:max_plots, 1:max_bins) = histxl
      deallocate(histxl)
      allocate(histxl(label_max, nbin_max))
      histxl = rtemp2

      rtemp2 = 0d0
      rtemp2(1:max_plots, 1:max_bins) = histxm
      deallocate(histxm)
      allocate(histxm(label_max, nbin_max))
      histxm = rtemp2

      rtemp2 = 0d0
      rtemp2(1:max_plots, 1:max_bins) = histy2
      deallocate(histy2)
      allocate(histy2(label_max, nbin_max))
      histy2 = rtemp2

      rtemp2 = 0d0
      rtemp2(1:max_plots, 1:max_bins) = histy_err
      deallocate(histy_err)
      allocate(histy_err(label_max, nbin_max))
      histy_err = rtemp2
      deallocate(rtemp2)

      allocate(rtemp3(nwgts, label_max, nbin_max))
      rtemp3 = 0d0
      rtemp3(:, 1:max_plots, 1:max_bins) = histy
      deallocate(histy)
      allocate(histy(nwgts, label_max, nbin_max))
      histy = rtemp3

      rtemp3 = 0d0
      rtemp3(:, 1:max_plots, 1:max_bins) = histy_acc
      deallocate(histy_acc)
      allocate(histy_acc(nwgts, label_max, nbin_max))
      histy_acc = rtemp3
      deallocate(rtemp3)

      allocate(itemp2(label_max, nbin_max))
      itemp2 = 0
      itemp2(1:max_plots, 1:max_bins) = histi
      deallocate(histi)
      allocate(histi(label_max, nbin_max))
      histi = itemp2
      deallocate(itemp2)

      max_plots = label_max
      max_bins = nbin_max
    else if (booked(label)) then
      write (*, *) 'ERROR in HwU.f90: histogram already booked', label
      stop
    end if
  end subroutine HwU_allocate_histo


  subroutine accum(idummy)
    integer, intent(in) :: idummy
    if (idummy /= 0) return
  end subroutine accum


  subroutine addfil(string)
    character(len=*), intent(in) :: string
    if (len(string) /= 0) return
  end subroutine addfil

end module HwU_module
