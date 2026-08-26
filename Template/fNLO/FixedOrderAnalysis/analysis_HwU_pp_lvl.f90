module analysis_hwu_pp_lvl_module
  use process_dimensions, only: nexternal
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  implicit none
  private
  public :: analysis_begin, analysis_end, analysis_fill

contains

!
! Example analysis for "p p > e+ ve [QCD]" process.
! Example analysis for "p p > e- ve~ [QCD]" process.
! Example analysis for "p p > mu+ vm [QCD]" process.
! Example analysis for "p p > mu- vm~ [QCD]" process.
! Example analysis for "p p > ta+ vt [QCD]" process.
! Example analysis for "p p > ta- vt~ [QCD]" process.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_begin(nwgt, weights_info)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer nwgt
    character(len=*) weights_info(*)
    integer i, l
    character(len=6) cc(2)
    data cc/'|T@NLO', '|T@LO '/
    call HwU_inithist(nwgt, weights_info)
    do i = 1, 2
      l = (i - 1)*8
      call HwU_book(l + 1, 'total rate '//cc(i), 5, 0.5d0, 5.5d0)
      call HwU_book(l + 2, 'lep rapidity '//cc(i), 20, -5d0, 5d0)
      call HwU_book(l + 3, 'lep pt '//cc(i), 20, 0d0, 200d0)
      call HwU_book(l + 4, 'et miss '//cc(i), 20, 0d0, 200d0)
      call HwU_book(l + 5, 'trans. mass '//cc(i), 40, 0d0, 200d0)
      call HwU_book(l + 6, 'w rapidity '//cc(i), 20, -5d0, 5d0)
      call HwU_book(l + 7, 'w pt '//cc(i), 20, 0d0, 200d0)
      call HwU_book(l + 8, 'cphi[l,vl] '//cc(i), 40, -1d0, 1d0)
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

    integer iPDG(nexternal)
    double precision p(0:4, nexternal)
    double precision wgts(*)
    integer ibody
    double precision var
    integer i, l
    double precision pw(0:3), pe(0:3), pn(0:3), ye, yw, pte, etmiss, mtr, ptw, cphi
    if (nexternal .ne. 5) then
      write (*, *) 'error #1 in analysis_fill: '// &
      & 'only for process "p p > l vl [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(1)) .le. 4 .or. ipdg(1) .eq. 21)) then
      write (*, *) 'error #2 in analysis_fill: '// &
      & 'only for process "p p > l vl [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(2)) .le. 4 .or. ipdg(2) .eq. 21)) then
      write (*, *) 'error #3 in analysis_fill: '// &
      & 'only for process "p p > l vl [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(5)) .le. 4 .or. ipdg(5) .eq. 21)) then
      write (*, *) 'error #4 in analysis_fill: '// &
      & 'only for process "p p > l vl [QCD]"'
      stop 1
    end if
    if ((abs(abs(ipdg(3)) - abs(ipdg(4))) .ne. 1) .or. &
    & (sign(1, ipdg(3)) .eq. sign(1, ipdg(4))) .or. &
    & (abs(ipdg(3)) .le. 10 .or. abs(ipdg(3)) .ge. 16) .or. &
    & (abs(ipdg(4)) .le. 10 .or. abs(ipdg(4)) .ge. 16)) then
      write (*, *) 'analysis not suited for this process', ipdg(3), ipdg(4)
    end if
    do i = 0, 3
      if (abs(ipdg(3)) .eq. 11 .or. abs(ipdg(3)) .eq. 13 .or. abs(ipdg(3)) .eq. 15) then
        pe(i) = p(i, 3)
        pn(i) = p(i, 4)
      else
        pe(i) = p(i, 4)
        pn(i) = p(i, 3)
      end if
      pw(i) = pe(i) + pn(i)
    end do
    ye = getrapidity(pe(0), pe(3))
    yw = getrapidity(pw(0), pw(3))
    pte = dsqrt(pe(1)**2 + pe(2)**2)
    ptw = dsqrt(pw(1)**2 + pw(2)**2)
    etmiss = dsqrt(pn(1)**2 + pn(2)**2)
    mtr = dsqrt(2d0*pte*etmiss - 2d0*pe(1)*pn(1) - 2d0*pe(2)*pn(2))
    cphi = (pe(1)*pn(1) + pe(2)*pn(2))/pte/etmiss
    var = 1.d0
    do i = 1, 2
      l = (i - 1)*8
      if (ibody .ne. 3 .and. i .eq. 2) cycle
      call HwU_fill(l + 1, var, wgts)
      call HwU_fill(l + 2, ye, wgts)
      call HwU_fill(l + 3, pte, wgts)
      call HwU_fill(l + 4, etmiss, wgts)
      call HwU_fill(l + 5, mtr, wgts)
      call HwU_fill(l + 6, yw, wgts)
      call HwU_fill(l + 7, ptw, wgts)
      call HwU_fill(l + 8, cphi, wgts)
    end do
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

end module analysis_hwu_pp_lvl_module
