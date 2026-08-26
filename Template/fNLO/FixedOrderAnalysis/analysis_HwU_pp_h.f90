module analysis_hwu_pp_h_module
  use process_dimensions, only: nexternal
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  implicit none
  private
  public :: analysis_begin, analysis_end, analysis_fill

contains

!
! Example analysis for "p p > h [QCD]" process.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_begin(nwgt, weights_info)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer nwgt
    character(len=*) weights_info(*)
    integer l
    l = 0
    call HwU_inithist(nwgt, weights_info)
    call HwU_book(l + 1, 'Higgs pT ', 100, 0.d0, 200.d0)
    call HwU_book(l + 2, 'Higgs pT ', 100, 0.d0, 500.d0)
    call HwU_book(l + 3, 'Higgs log[pT] ', 98, 0.1d0, 5.d0)
    call HwU_book(l + 4, 'Higgs pT,abs(y_H)<2 ', 100, 0.d0, 200.d0)
    call HwU_book(l + 5, 'Higgs pT,abs(y_H)<2 ', 100, 0.d0, 500.d0)
    call HwU_book(l + 6, 'Higgs log[pT],abs(y_H)<2 ', 98, 0.1d0, 5.d0)

    call HwU_book(l + 7, 'j1 pT ', 100, 0.d0, 200.d0)
    call HwU_book(l + 8, 'j1 pT ', 100, 0.d0, 500.d0)
    call HwU_book(l + 9, 'j1 log[pT] ', 98, 0.1d0, 5.d0)
    call HwU_book(l + 10, 'j1 pT,abs(y_j1)<2 ', 100, 0.d0, 200.d0)
    call HwU_book(l + 11, 'j1 pT,abs(y_j1)<2 ', 100, 0.d0, 500.d0)
    call HwU_book(l + 12, 'j1 log[pT],abs(y_j1)<2 ', 98, 0.1d0, 5.d0)

    call HwU_book(l + 13, 'Inc j pT ', 100, 0.d0, 200.d0)
    call HwU_book(l + 14, 'Inc j pT ', 100, 0.d0, 500.d0)
    call HwU_book(l + 15, 'Inc j log[pT] ', 98, 0.1d0, 5.d0)
    call HwU_book(l + 16, 'Inc j pT,abs(y_Ij)<2 ', 100, 0.d0, 2.d2)
    call HwU_book(l + 17, 'Inc j pT,abs(y_Ij)<2 ', 100, 0.d0, 5.d2)
    call HwU_book(l + 18, 'Inc j log[pT],abs(y_Ij)<2', 98, 0.1d0, 5.d0)

    call HwU_book(l + 19, 'Higgs y ', 60, -6.d0, 6.d0)
    call HwU_book(l + 20, 'Higgs y,pT_H>10GeV ', 100, -6.d0, 6.d0)
    call HwU_book(l + 21, 'Higgs y,pT_H>30GeV ', 100, -6.d0, 6.d0)
    call HwU_book(l + 22, 'Higgs y,pT_H>50GeV ', 100, -6.d0, 6.d0)
    call HwU_book(l + 23, 'Higgs y,pT_H>70GeV ', 100, -6.d0, 6.d0)
    call HwU_book(l + 24, 'Higgs y,pt_H>90GeV ', 100, -6.d0, 6.d0)

    call HwU_book(l + 25, 'j1 y ', 60, -6.d0, 6.d0)
    call HwU_book(l + 26, 'j1 y,pT_j1>10GeV ', 100, -6.d0, 6.d0)
    call HwU_book(l + 27, 'j1 y,pT_j1>30GeV ', 100, -6.d0, 6.d0)
    call HwU_book(l + 28, 'j1 y,pT_j1>50GeV ', 100, -6.d0, 6.d0)
    call HwU_book(l + 29, 'j1 y,pT_j1>70GeV ', 100, -6.d0, 6.d0)
    call HwU_book(l + 30, 'j1 y,pT_j1>90GeV ', 100, -6.d0, 6.d0)

    call HwU_book(l + 31, 'H-j1 y ', 60, -6.d0, 6.d0)
    call HwU_book(l + 32, 'H-j1 y,pT_j1>10GeV ', 100, -6.d0, 6.d0)
    call HwU_book(l + 33, 'H-j1 y,pT_j1>30GeV ', 100, -6.d0, 6.d0)
    call HwU_book(l + 34, 'H-j1 y,pT_j1>50GeV ', 100, -6.d0, 6.d0)
    call HwU_book(l + 35, 'H-j1 y,pT_j1>70GeV ', 100, -6.d0, 6.d0)
    call HwU_book(l + 36, 'H-j1 y,pT_j1>90GeV ', 100, -6.d0, 6.d0)

    call HwU_book(l + 37, 'njets ', 11, -0.5d0, 10.5d0)
    call HwU_book(l + 38, 'njets,abs(y_j)<2.5 ', 11, -0.5d0, 10.5d0)
    call HwU_book(l + 39, 'xsec ', 3, -0.5d0, 2.5d0)
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
    integer i, l
    double precision pph(0:3), ppj1(0:3), pth, yh, ptj1, yj1, y_central &
    & , ptcut, njdble, njcdble
    integer njet, njet_central, nj
    if (nexternal .ne. 4) then
      write (*, *) 'error #1 in analysis_fill: '// &
      & 'only for process "p p > h [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(1)) .le. 5 .or. ipdg(1) .eq. 21)) then
      write (*, *) 'error #2 in analysis_fill: '// &
      & 'only for process "p p > h [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(2)) .le. 5 .or. ipdg(2) .eq. 21)) then
      write (*, *) 'error #3 in analysis_fill: '// &
      & 'only for process "p p > h [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(4)) .le. 5 .or. ipdg(4) .eq. 21)) then
      write (*, *) 'error #4 in analysis_fill: '// &
      & 'only for process "p p > h [QCD]"'
      stop 1
    end if
    if (ipdg(3) .ne. 25) then
      write (*, *) 'error #5 in analysis_fill: '// &
      & 'only for process "p p > h [QCD]"'
      stop 1
    end if
!
    do i = 0, 3
      pph(i) = p(i, 3)
      ppj1(i) = p(i, 4)
    end do
! Higgs variables
    pth = sqrt(pph(1)**2 + pph(2)**2)
    yh = getrapidity(pph(0), pph(3))
! hardest jet variables
    ptj1 = sqrt(ppj1(1)**2 + ppj1(2)**2)
    yj1 = getrapidity(ppj1(0), ppj1(3))
    njet = 0
    njet_central = 0
    y_central = 2.5d0
    ptcut = 1d1
    if (ptj1 .gt. ptcut) njet = 1
    if (ptj1 .gt. ptcut .and. abs(yj1) .le. y_central) njet_central = 1
    njdble = dble(njet)
    njcdble = dble(njet_central)
!
    l = 0
    call HwU_fill(l + 1, pth, WGTS)
    call HwU_fill(l + 2, pth, WGTS)
    if (pth .gt. 0.d0) call HwU_fill(l + 3, log10(pth), WGTS)
    if (abs(yh) .le. 2.d0) then
      call HwU_fill(l + 4, pth, WGTS)
      call HwU_fill(l + 5, pth, WGTS)
      if (pth .gt. 0.d0) call HwU_fill(l + 6, log10(pth), WGTS)
    end if
!
    if (njet .ge. 1) then
      call HwU_fill(l + 7, ptj1, WGTS)
      call HwU_fill(l + 8, ptj1, WGTS)
      if (ptj1 .gt. 0.d0) call HwU_fill(l + 9, log10(ptj1), WGTS)
      if (abs(yj1) .le. 2.d0) then
        call HwU_fill(l + 10, ptj1, WGTS)
        call HwU_fill(l + 11, ptj1, WGTS)
        if (ptj1 .gt. 0.d0) call HwU_fill(l + 12, log10(ptj1), WGTS)
      end if
!
      do nj = 1, njet
        call HwU_fill(l + 13, ptj1, WGTS)
        call HwU_fill(l + 14, ptj1, WGTS)
        if (ptj1 .gt. 0.d0) call HwU_fill(l + 15, log10(ptj1), WGTS)
        if (abs(yj1) .le. 2.d0) then
          call HwU_fill(l + 16, ptj1, WGTS)
          call HwU_fill(l + 17, ptj1, WGTS)
          if (ptj1 .gt. 0d0) call HwU_fill(l + 18, log10(ptj1), WGTS)
        end if
      end do
    end if
!
    call HwU_fill(l + 19, yh, WGTS)
    if (pth .ge. 10.d0) call HwU_fill(l + 20, yh, WGTS)
    if (pth .ge. 30.d0) call HwU_fill(l + 21, yh, WGTS)
    if (pth .ge. 50.d0) call HwU_fill(l + 22, yh, WGTS)
    if (pth .ge. 70.d0) call HwU_fill(l + 23, yh, WGTS)
    if (pth .ge. 90.d0) call HwU_fill(l + 24, yh, WGTS)
!
    if (njet .ge. 1) then
      call HwU_fill(l + 25, yj1, WGTS)
      if (ptj1 .ge. 10.d0) call HwU_fill(l + 26, yj1, WGTS)
      if (ptj1 .ge. 30.d0) call HwU_fill(l + 27, yj1, WGTS)
      if (ptj1 .ge. 50.d0) call HwU_fill(l + 28, yj1, WGTS)
      if (ptj1 .ge. 70.d0) call HwU_fill(l + 29, yj1, WGTS)
      if (ptj1 .ge. 90.d0) call HwU_fill(l + 30, yj1, WGTS)
!
      call HwU_fill(l + 31, yh - yj1, WGTS)
      if (ptj1 .ge. 10.d0) call HwU_fill(l + 32, yh - yj1, WGTS)
      if (ptj1 .ge. 30.d0) call HwU_fill(l + 33, yh - yj1, WGTS)
      if (ptj1 .ge. 50.d0) call HwU_fill(l + 34, yh - yj1, WGTS)
      if (ptj1 .ge. 70.d0) call HwU_fill(l + 35, yh - yj1, WGTS)
      if (ptj1 .ge. 90.d0) call HwU_fill(l + 36, yh - yj1, WGTS)
    end if
!
    call HwU_fill(l + 37, njdble, WGTS)
    call HwU_fill(l + 38, njcdble, WGTS)
    call HwU_fill(l + 39, 1d0, WGTS)
!
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
end module analysis_hwu_pp_h_module
