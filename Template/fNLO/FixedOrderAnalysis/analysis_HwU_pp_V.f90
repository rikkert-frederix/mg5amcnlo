module analysis_hwu_pp_v_module
  use process_dimensions, only: event_capacity
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  use kin_functions_module, only: dot => dot_impl
  implicit none
  private
  public :: analysis_begin, analysis_end, analysis_fill

contains

!
! Example analysis for "p p > w+ [QCD]" process.
! Example analysis for "p p > w- [QCD]" process.
! Example analysis for "p p > z [QCD]" process.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_begin(nwgt, weights_info)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer nwgt
    character(len=*) weights_info(*)
    integer j, l
    character(len=6) cc(2)
    data cc/'|T@NLO', '|T@LO '/
    double precision xmi, xms
    call HwU_inithist(nwgt, weights_info)
    xmi = 40.d0
    xms = 120.d0
    do j = 1, 2
      l = (j - 1)*5
      call HwU_book(l + 1, 'V pt     '//cc(j), 100, 0.d0, 200.d0)
      call HwU_book(l + 2, 'V log pt '//cc(j), 100, 0.d0, 5.d0)
      call HwU_book(l + 3, 'V y      '//cc(j), 78, -9.d0, 9.d0)
      call HwU_book(l + 4, 'V eta    '//cc(j), 78, -9.d0, 9.d0)
      call HwU_book(l + 5, 'mV       '//cc(j), 80, xmi, xms)
    end do
    return
  end subroutine analysis_begin

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_end()
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    use open_output_files_module, only: HwU_write_file
    implicit none
    call HwU_write_file
    return
  end subroutine analysis_end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_fill(p, ipdg, wgts, ibody)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none

    integer iPDG(event_capacity)
    double precision p(0:4, event_capacity)
    double precision wgts(*)
    integer ibody
    integer i, l
    double precision pv(0:3), xmv, ptv, yv, etav
    if (event_capacity .ne. 4) then
      write (*, *) 'error #1 in analysis_fill: '// &
      & 'only for process "p p > V [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(1)) .le. 5 .or. ipdg(1) .eq. 21)) then
      write (*, *) 'error #2 in analysis_fill: '// &
      & 'only for process "p p > V [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(2)) .le. 5 .or. ipdg(2) .eq. 21)) then
      write (*, *) 'error #3 in analysis_fill: '// &
      & 'only for process "p p > V [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(4)) .le. 5 .or. ipdg(4) .eq. 21)) then
      write (*, *) 'error #4 in analysis_fill: '// &
      & 'only for process "p p > V [QCD]"'
      stop 1
    end if
    if (abs(ipdg(3)) .ne. 24 .and. ipdg(3) .ne. 23) then
      write (*, *) 'error #5 in analysis_fill: '// &
      & 'only for process "p p > V [QCD]"'
      stop 1
    end if
!
    do i = 0, 3
      pv(i) = p(i, 3)
    end do
    xmv = sqrt(max(dot(pv, pv), 0d0))
    ptv = sqrt(max(pv(1)**2 + pv(2)**2, 0d0))
    yv = getrapidity(pv(0), pv(3))
    etav = getpseudorap(pv(1), pv(2), pv(3))
!
    do i = 1, 2
      l = (i - 1)*5
      if (ibody .ne. 3 .and. i .eq. 2) cycle
      call HwU_fill(l + 1, ptv, WGTS)
      if (ptv .gt. 0) call HwU_fill(l + 2, log10(ptv), WGTS)
      call HwU_fill(l + 3, yv, WGTS)
      call HwU_fill(l + 4, etav, WGTS)
      call HwU_fill(l + 5, xmv, WGTS)
    end do
!
    return
  end subroutine analysis_fill

  function getrapidity(en, pl)
    implicit none
    double precision getrapidity, en, pl, tiny, xplus, xminus, y
    parameter(tiny=1.d-8)
    xplus = en + pl
    xminus = en - pl
    if (xplus .gt. tiny .and. xminus .gt. tiny) then
      if ((xplus/xminus) .gt. tiny .and. (xminus/xplus) .gt. tiny) then
        y = 0.5d0*log(xplus/xminus)
      else
        y = sign(1.d0, pl)*1.d8
      end if
    else
      y = sign(1.d0, pl)*1.d8
    end if
    getrapidity = y
    return
  end function getrapidity

  function getpseudorap(ptx, pty, pl)
    implicit none
    double precision getpseudorap, ptx, pty, pl, tiny, pt, eta, th
    parameter(tiny=1.d-5)
!
    pt = sqrt(ptx**2 + pty**2)
    if (pt .lt. tiny .and. abs(pl) .lt. tiny) then
      eta = sign(1.d0, pl)*1.d8
    else
      th = atan2(pt, pl)
      eta = -log(tan(th/2.d0))
    end if
    getpseudorap = eta
    return
  end function getpseudorap

end module analysis_hwu_pp_v_module
