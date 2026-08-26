module analysis_hwu_pp_lplm_module
  use process_dimensions, only: nexternal
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  implicit none
  private
  public :: analysis_begin, analysis_end, analysis_fill

contains

!
! Example analysis for "p p > l+ l- [QCD]" process.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_begin(nwgt, weights_info)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer nwgt
    character(len=*) weights_info(*)
    integer j, l
    character(len=5) cc(2)
    data cc/'     ', 'cuts '/
    double precision xmi, xms, pi
    parameter(pi=3.141592653589793d0)
    call HwU_inithist(nwgt, weights_info)
    xmi = 50.d0
    xms = 130.d0
    do j = 1, 2
      l = (j - 1)*21
      call HwU_book(l + 1, 'V pt      '//cc(j), 100, 0.d0, 200.d0)
      call HwU_book(l + 2, 'V pt      '//cc(j), 100, 0.d0, 1000.d0)
      call HwU_book(l + 3, 'V log[pt] '//cc(j), 98, 0.1d0, 5.d0)
      call HwU_book(l + 4, 'V y       '//cc(j), 90, -9.d0, 9.d0)
      call HwU_book(l + 5, 'V eta     '//cc(j), 90, -9.d0, 9.d0)
      call HwU_book(l + 6, 'mV        '//cc(j), 100, xmi, xms)
!
      call HwU_book(l + 7, 'lm pt      '//cc(j), 100, 0.d0, 200.d0)
      call HwU_book(l + 8, 'lm pt      '//cc(j), 100, 0.d0, 1000.d0)
      call HwU_book(l + 9, 'lm log[pt] '//cc(j), 98, 0.1d0, 5.d0)
      call HwU_book(l + 10, 'lm eta     '//cc(j), 90, -9.d0, 9.d0)
      call HwU_book(l + 11, 'lp pt      '//cc(j), 100, 0.d0, 200.d0)
      call HwU_book(l + 12, 'lp pt      '//cc(j), 100, 0.d0, 1000.d0)
      call HwU_book(l + 13, 'lp log[pt] '//cc(j), 98, 0.1d0, 5.d0)
      call HwU_book(l + 14, 'lp eta     '//cc(j), 90, -9.d0, 9.d0)
!
      call HwU_book(l + 15, 'lmlp delta eta     '//cc(j), 90, -9.d0, 9.d0)
      call HwU_book(l + 16, 'lmlp azimt         '//cc(j), 20, 0.d0, pi)
      call HwU_book(l + 17, 'lmlp log[pi-azimt] '//cc(j), 82, -4.d0, 0.1d0)
      call HwU_book(l + 18, 'lmlp inv m         '//cc(j), 100, xmi, xms)
      call HwU_book(l + 19, 'lmlp pt            '//cc(j), 100, 0.d0, 200.d0)
      call HwU_book(l + 20, 'lmlp log[pt]       '//cc(j), 98, 0.1d0, 5.d0)
!
      call HwU_book(l + 21, 'total'//cc(j), 2, -1.d0, 1.d0)
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
  subroutine analysis_fill(p, ipdg, wgts)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none

    integer iPDG(nexternal)
    double precision p(0:4, nexternal)
    double precision wgts(*)
    integer i, l, nwgt_analysis
    common/c_analysis/nwgt_analysis
    double precision ppl(0:3), pplb(0:3), ppv(0:3), ycut, xmv, ptv, yv &
    & , etav, ptl, yl, etal, ptlb, ylb, etalb, ptpair, azi, azinorm, xmll &
    & , detallb
    double precision pi
    parameter(pi=3.141592653589793d0)
    if (nexternal .ne. 5) then
      write (*, *) 'error #1 in analysis_fill: '// &
      & 'only for process "p p > l+ l- [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(1)) .le. 5 .or. ipdg(1) .eq. 21)) then
      write (*, *) 'error #2 in analysis_fill: '// &
      & 'only for process "p p > l+ l- [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(2)) .le. 5 .or. ipdg(2) .eq. 21)) then
      write (*, *) 'error #3 in analysis_fill: '// &
      & 'only for process "p p > l+ l- [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(5)) .le. 5 .or. ipdg(5) .eq. 21)) then
      write (*, *) 'error #4 in analysis_fill: '// &
      & 'only for process "p p > l+ l- [QCD]"'
      stop 1
    end if
    if (ipdg(3) .ne. -11 .and. ipdg(3) .ne. -13 .and. ipdg(3) .ne. -15) then
      write (*, *) 'error #5 in analysis_fill: '// &
      & 'only for process "p p > l+ l- [QCD]"'
      stop 1
    end if
    if (ipdg(4) .ne. 11 .and. ipdg(4) .ne. 13 .and. ipdg(4) .ne. 15) then
      write (*, *) 'error #6 in analysis_fill: '// &
      & 'only for process "p p > l+ l- [QCD]"'
      stop 1
    end if

    do i = 0, 3
      ppl(i) = p(i, 4)
      pplb(i) = p(i, 3)
      ppv(i) = ppl(i) + pplb(i)
    end do

! FILL THE HISTOS
    YCUT = 2.5d0
! Variables of the vector boson
    xmv = getinvm(ppv(0), ppv(1), ppv(2), ppv(3))
    ptv = sqrt(ppv(1)**2 + ppv(2)**2)
    yv = getrapidity(ppv(0), ppv(3))
    etav = getpseudorap(ppv(1), ppv(2), ppv(3))
! Variables of the leptons
    ptl = sqrt(ppl(1)**2 + ppl(2)**2)
    yl = getrapidity(ppl(0), ppl(3))
    etal = getpseudorap(ppl(1), ppl(2), ppl(3))
!
    ptlb = sqrt(pplb(1)**2 + pplb(2)**2)
    ylb = getrapidity(pplb(0), pplb(3))
    etalb = getpseudorap(pplb(1), pplb(2), pplb(3))
!
    ptpair = ptv
    azi = getdelphi(ppl(1), ppl(2), pplb(1), pplb(2))
    azinorm = (pi - azi)/pi
    xmll = xmv
    detallb = etal - etalb
!
    l = 0
    call HwU_fill(l + 1, (ptv), WGTS)
    call HwU_fill(l + 2, (ptv), WGTS)
    if (ptv .gt. 0.d0) call HwU_fill(l + 3, (log10(ptv)), WGTS)
    call HwU_fill(l + 4, (yv), WGTS)
    call HwU_fill(l + 5, (etav), WGTS)
    call HwU_fill(l + 6, (xmv), WGTS)
!
    call HwU_fill(l + 7, (ptl), WGTS)
    call HwU_fill(l + 8, (ptl), WGTS)
    if (ptl .gt. 0.d0) call HwU_fill(l + 9, (log10(ptl)), WGTS)
    call HwU_fill(l + 10, (etal), WGTS)
    call HwU_fill(l + 11, (ptlb), WGTS)
    call HwU_fill(l + 12, (ptlb), WGTS)
    if (ptlb .gt. 0.d0) call HwU_fill(l + 13, (log10(ptlb)), WGTS)
    call HwU_fill(l + 14, (etalb), WGTS)
!
    call HwU_fill(l + 15, (detallb), WGTS)
    call HwU_fill(l + 16, (azi), WGTS)
    if (azinorm .gt. 0.d0) call HwU_fill(l + 17, (log10(azinorm)), WGTS)
    call HwU_fill(l + 18, (xmll), WGTS)
    call HwU_fill(l + 19, (ptpair), WGTS)
    if (ptpair .gt. 0) call HwU_fill(l + 20, (log10(ptpair)), WGTS)
    call HwU_fill(l + 21, (0d0), WGTS)
!
    l = l + 21

    if (abs(etav) .lt. ycut) then
      call HwU_fill(l + 1, (ptv), WGTS)
      call HwU_fill(l + 2, (ptv), WGTS)
      if (ptv .gt. 0.d0) call HwU_fill(l + 3, (log10(ptv)), WGTS)
    end if
    if (ptv .gt. 20.d0) then
      call HwU_fill(l + 4, (yv), WGTS)
      call HwU_fill(l + 5, (etav), WGTS)
    end if
    if (abs(etav) .lt. ycut .and. ptv .gt. 20.d0) then
      call HwU_fill(l + 6, (xmv), WGTS)
      call HwU_fill(l + 21, (0d0), WGTS)
    end if
!
    if (abs(etal) .lt. ycut) then
      call HwU_fill(l + 7, (ptl), WGTS)
      call HwU_fill(l + 8, (ptl), WGTS)
      if (ptl .gt. 0.d0) call HwU_fill(l + 9, (log10(ptl)), WGTS)
    end if
    if (ptl .gt. 20.d0) call HwU_fill(l + 10, (etal), WGTS)
    if (abs(etalb) .lt. ycut) then
      call HwU_fill(l + 11, (ptlb), WGTS)
      call HwU_fill(l + 12, (ptlb), WGTS)
      if (ptlb .gt. 0.d0) call HwU_fill(l + 13, (log10(ptlb)), WGTS)
    end if
    if (ptlb .gt. 20.d0) call HwU_fill(l + 14, (etalb), WGTS)
!
    if (abs(etal) .lt. ycut .and. abs(etalb) .lt. ycut .and. &
    & ptl .gt. 20.d0 .and. ptlb .gt. 20.d0) then
      call HwU_fill(l + 15, (detallb), WGTS)
      call HwU_fill(l + 16, (azi), WGTS)
      if (azinorm .gt. 0.d0) call HwU_fill(l + 17, (log10(azinorm)), WGTS)
      call HwU_fill(l + 18, (xmll), WGTS)
      call HwU_fill(l + 19, (ptpair), WGTS)
      if (ptpair .gt. 0) call HwU_fill(l + 20, (log10(ptpair)), WGTS)
    end if

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

  function getinvm(en, ptx, pty, pl)
    implicit none
    double precision getinvm, en, ptx, pty, pl, tiny, tmp
    parameter(tiny=1.d-5)
!
    tmp = en**2 - ptx**2 - pty**2 - pl**2
    if (tmp .gt. 0.d0) then
      tmp = sqrt(tmp)
    elseif (tmp .gt. -tiny) then
      tmp = 0.d0
    else
      write (*, *) 'Attempt to compute a negative mass'
      stop
    end if
    getinvm = tmp
    return
  end function getinvm

  function getdelphi(ptx1, pty1, ptx2, pty2)
    implicit none
    double precision getdelphi, ptx1, pty1, ptx2, pty2, tiny, pt1, pt2, tmp
    parameter(tiny=1.d-5)
!
    pt1 = sqrt(ptx1**2 + pty1**2)
    pt2 = sqrt(ptx2**2 + pty2**2)
    if (pt1 .gt. 0.d0 .and. pt2 .gt. 0.d0) then
      tmp = ptx1*ptx2 + pty1*pty2
      tmp = tmp/(pt1*pt2)
      if (abs(tmp) .gt. 1.d0 + tiny) then
        write (*, *) 'Cosine larger than 1'
        stop
      elseif (abs(tmp) .ge. 1.d0) then
        tmp = sign(1.d0, tmp)
      end if
      tmp = acos(tmp)
    else
      tmp = 1.d8
    end if
    getdelphi = tmp
    return
  end function getdelphi

end module analysis_hwu_pp_lplm_module
