module analysis_hwu_pp_hjj_module
  use process_dimensions, only: event_capacity, nincoming
  use HwU_module, only: HwU_inithist, HwU_book, HwU_fill
  implicit none
  private

  public :: analysis_begin, analysis_end, analysis_fill

  double precision, parameter :: vetomin = 0d0
  double precision, parameter :: vetomax = 100d0
  integer, parameter :: nbinveto = 50
  double precision, allocatable :: veto_weights(:)

  interface
    subroutine amcatnlo_fastjetppgenkt(pqcd, nn, rfj, sycut, palg, &
                                       pjet, njet, jet)
      implicit none
      integer, intent(in) :: nn
      double precision, intent(in) :: pqcd(0:3, nn)
      double precision, intent(in) :: rfj, sycut, palg
      double precision, intent(out) :: pjet(0:3, nn)
      integer, intent(out) :: njet, jet(nn)
    end subroutine amcatnlo_fastjetppgenkt
  end interface

contains

!
! Example analysis for "p p > h j j [QCD]" process.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine analysis_begin(nwgt, weights_info)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    implicit none
    integer, intent(in) :: nwgt
    character(len=*), intent(in) :: weights_info(*)
    integer i, l
    character(len=8) cc(2)
    data cc/'        ', 'vbfcuts '/
    double precision pi
    parameter(PI=3.141592653589793d0)
!
    if (allocated(veto_weights)) deallocate (veto_weights)
    allocate (veto_weights(nwgt))
    call HwU_inithist(nwgt, weights_info)
    do i = 1, 2
      l = (i - 1)*54
      call HwU_book(l + 1, 'total '//cc(i), 5, 0.5d0, 5.5d0)

      call HwU_book(l + 2, 'Higgs pT '//cc(i), 50, 0.d0, 400.d0)
      call HwU_book(l + 3, 'Higgs pT '//cc(i), 50, 0.d0, 800.d0)
      call HwU_book(l + 4, 'Higgs logpT '//cc(i), 50, 0.d0, 4.d0)
      call HwU_book(l + 5, 'Higgs eta '//cc(i), 50, -6.d0, 6.d0)
      call HwU_book(l + 6, 'Higgs y '//cc(i), 50, -6.d0, 6.d0)

      call HwU_book(l + 7, 'j1 pT '//cc(i), 50, 0.d0, 400.d0)
      call HwU_book(l + 8, 'j1 pT '//cc(i), 50, 0.d0, 800.d0)
      call HwU_book(l + 9, 'j1 logpT '//cc(i), 50, 0.d0, 4.d0)
      call HwU_book(l + 10, 'j1 eta '//cc(i), 50, -6.d0, 6.d0)
      call HwU_book(l + 11, 'j1 y '//cc(i), 50, -6.d0, 6.d0)

      call HwU_book(l + 12, 'j2 pT '//cc(i), 50, 0.d0, 400.d0)
      call HwU_book(l + 13, 'j2 pT '//cc(i), 50, 0.d0, 800.d0)
      call HwU_book(l + 14, 'j2 logpT '//cc(i), 50, 0.d0, 4.d0)
      call HwU_book(l + 15, 'j2 eta '//cc(i), 50, -6.d0, 6.d0)
      call HwU_book(l + 16, 'j2 y '//cc(i), 50, -6.d0, 6.d0)

      call HwU_book(l + 17, 'j3 pT '//cc(i), 50, 0.d0, 400.d0)
      call HwU_book(l + 18, 'j3 pT '//cc(i), 50, 0.d0, 800.d0)
      call HwU_book(l + 19, 'j3 logpT '//cc(i), 50, 0.d0, 4.d0)
      call HwU_book(l + 20, 'j3 eta '//cc(i), 50, -6.d0, 6.d0)
      call HwU_book(l + 21, 'j3 y '//cc(i), 50, -6.d0, 6.d0)

      call HwU_book(l + 22, 'H+j1 pT '//cc(i), 50, 0.d0, 400.d0)
      call HwU_book(l + 23, 'H+j1 pT '//cc(i), 50, 0.d0, 800.d0)
      call HwU_book(l + 24, 'H+j1 logpT '//cc(i), 50, 0.d0, 4.d0)
      call HwU_book(l + 25, 'H+j1 eta '//cc(i), 50, -6.d0, 6.d0)
      call HwU_book(l + 26, 'H+j1 y '//cc(i), 50, -6.d0, 6.d0)

      call HwU_book(l + 27, 'j1+j2 pT '//cc(i), 50, 0.d0, 400.d0)
      call HwU_book(l + 28, 'j1+j2 pT '//cc(i), 50, 0.d0, 800.d0)
      call HwU_book(l + 29, 'j1+j2 logpT '//cc(i), 50, 0.d0, 4.d0)
      call HwU_book(l + 30, 'j1+j2 eta '//cc(i), 50, -6.d0, 6.d0)
      call HwU_book(l + 31, 'j1+j2 y '//cc(i), 50, -6.d0, 6.d0)

      call HwU_book(l + 32, 'syst pT '//cc(i), 50, 0.d0, 400.d0)
      call HwU_book(l + 33, 'syst pT '//cc(i), 50, 0.d0, 800.d0)
      call HwU_book(l + 34, 'syst logpT '//cc(i), 50, 0.d0, 4.d0)
      call HwU_book(l + 35, 'syst eta '//cc(i), 50, -10.d0, 10.d0)
      call HwU_book(l + 36, 'syst y '//cc(i), 50, -6.d0, 6.d0)

      call HwU_book(l + 37, 'Dphi H-j1 '//cc(i), 50, 0d0, pi)
      call HwU_book(l + 38, 'Dphi H-j2 '//cc(i), 50, 0d0, pi)
      call HwU_book(l + 39, 'Dphi j1-j2 '//cc(i), 50, 0d0, pi)

      call HwU_book(l + 40, 'DR H-j1 '//cc(i), 50, 0d0, 10.d0)
      call HwU_book(l + 41, 'DR H-j2 '//cc(i), 50, 0d0, 10.d0)
      call HwU_book(l + 42, 'DR j1-j2 '//cc(i), 50, 0d0, 10.d0)

      call HwU_book(l + 43, 'mj1j2 '//cc(i), 50, 0d0, 3000.d0)

! Nason-Oleari plots (hep-ph/0911.5299)
      call HwU_book(l + 44, '|yj1-yj2| '//cc(i), 25, 0.d0, 10.d0)
      call HwU_book(l + 45, 'yj3_rel '//cc(i), 50, -6.d0, 6.d0)
      call HwU_book(l + 46, 'njets '//cc(i), 100, -0.5d0, 9.5d0)
      call HwU_book(l + 47, 'ptrel_j1 '//cc(i), 50, 0.d0, 200.d0)
      call HwU_book(l + 48, 'ptrel_j2 '//cc(i), 50, 0.d0, 200.d0)
      call HwU_book(l + 49, 'P-veto '//cc(i), nbinveto, vetomin &
      & , vetomax)
      call HwU_book(l + 50, 'jveto pT '//cc(i), nbinveto, vetomin &
      & , vetomax)
      call HwU_book(l + 51, 'jveto pT '//cc(i), nbinveto, &
      & vetomin, 2d0*vetomax)
      call HwU_book(l + 52, 'jveto logpT '//cc(i), 50, 0.d0, 4.d0)
      call HwU_book(l + 53, 'jveto eta '//cc(i), 50, -6.d0, 6.d0)
      call HwU_book(l + 54, 'jveto y '//cc(i), 50, -6.d0, 6.d0)

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
    use boostwdir2_module, only: boostwdir2
    implicit none
    integer, intent(in) :: iPDG(event_capacity)
    double precision, intent(in) :: p(0:4, event_capacity)
    double precision, intent(in) :: wgts(*)
    integer i, j, k, l
    double precision pQCD(0:3, event_capacity), palg, rfj, sycut, yjmax &
    & , pjet(0:3, event_capacity), ptjet(event_capacity), yjet(event_capacity) &
    & , etajet(event_capacity), ptj_tag, deltay12, mj1j2min, ph(0:3), pj1(0:3) &
    & , pj2(0:3), pj3(0:3), pjj(0:3), pHj(0:3), psyst(0:3), pjveto(0:3) &
    & , ptH, etaH, yH, njdble, ptj1, etaj1, yj1, ptHj, etaHj, yHj, DphiHj1 &
    & , DRHj1, ptj2, etaj2, yj2, ptjj, etajj, yjj, ptsyst, etasyst, ysyst &
    & , DphiHj2, Dphij1j2, DRHj2, DRj1j2, mj1j2, Dyj1j2, ptj3, etaj3, yj3 &
    & , yj3rel, chy1, shy1, chy1mo, chy2, shy2, chy2mo, ptrel_j1, ptrel_j2 &
    & , ppboost(0:3, event_capacity), prel_j1(0:3), prel_j2(0:3) &
    & , pj1boost(0:3), pj2boost(0:3), pt_veto, previous_pt
    logical pass_tag_cuts, flag
    integer nQCD, jet(event_capacity), ij1y, ij2y, ij3y, njet, njety, ijveto &
    & , ijvetoy, ij1, ij2, ij3
    double precision xd(1:3)
    data(xd(i), i=1, 3)/0d0, 0d0, 1d0/
    if (event_capacity .ne. 6) then
      write (*, *) 'error #1 in analysis_fill: '// &
      & 'only for process "p p > h j j [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(1)) .le. 5 .or. ipdg(1) .eq. 21)) then
      write (*, *) 'error #2 in analysis_fill: '// &
      & 'only for process "p p > h j j [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(2)) .le. 5 .or. ipdg(2) .eq. 21)) then
      write (*, *) 'error #3 in analysis_fill: '// &
      & 'only for process "p p > h j j [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(4)) .le. 5 .or. ipdg(4) .eq. 21)) then
      write (*, *) 'error #4 in analysis_fill: '// &
      & 'only for process "p p > h j j [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(5)) .le. 5 .or. ipdg(5) .eq. 21)) then
      write (*, *) 'error #5 in analysis_fill: '// &
      & 'only for process "p p > h j j [QCD]"'
      stop 1
    end if
    if (.not. (abs(ipdg(6)) .le. 5 .or. ipdg(6) .eq. 21)) then
      write (*, *) 'error #6 in analysis_fill: '// &
      & 'only for process "p p > h j j [QCD]"'
      stop 1
    end if
    if (ipdg(3) .ne. 25) then
      write (*, *) 'error #7 in analysis_fill: '// &
      & 'only for process "p p > h j j [QCD]"'
      stop 1
    end if
!

! Put all (light) QCD partons in momentum array for jet clustering.
    nQCD = 0
    do j = nincoming + 1, event_capacity
      if (abs(ipdg(j)) .le. 5 .or. ipdg(j) .eq. 21) then
        nQCD = nQCD + 1
        do i = 0, 3
          pQCD(i, nQCD) = p(i, j)
        end do
      end if
    end do

!---CLUSTER THE EVENT
    palg = -1.d0
    rfj = 0.4d0
    sycut = 20d0
    yjmax = 4.5d0
    do i = 1, event_capacity
      do j = 0, 3
        pjet(j, i) = 0d0
      end do
      ptjet(i) = 0d0
      yjet(i) = 0d0
      etajet(i) = 0d0
      jet(i) = 0
    end do
    ij1y = 0
    ij2y = 0
    ij3y = 0
    njet = 0
    njety = 0
    ijveto = 0
    ijvetoy = 0

!******************************************************************************
!     call FASTJET to get all the jets
!
!     INPUT:
!     input momenta:               pQCD(0:3,event_capacity), energy is 0th component
!     number of input momenta:     nQCD
!     radius parameter:            rfj
!     minumum jet pt:              sycut
!     jet algorithm:               palg, 1.0=kt, 0.0=C/A, -1.0 = anti-kt
!
!     OUTPUT:
!     jet momenta:                             pjet(0:3,event_capacity), E is 0th cmpnt
!     the number of jets (with pt > SYCUT):    njet
!     the jet for a given particle 'i':        jet(i),   note that this is
!     the particle in pQCD, which doesn't necessarily correspond to the particle
!     label in the process
!
    call amcatnlo_fastjetppgenkt(pQCD, nQCD, rfj, sycut, palg, pjet, njet &
    & , jet)
!
!******************************************************************************
    previous_pt = huge(previous_pt)
    do i = 1, njet
      ptjet(i) = getptv4(pjet(0, i))
      if (i .gt. 1) then
        if (ptjet(i) .gt. previous_pt) then
          write (*, *) "Error 1: jets should be ordered in pt"
          stop
        end if
      end if
      previous_pt = ptjet(i)
      yjet(i) = getrapidityv4(pjet(0, i))
      etajet(i) = getpseudorapv4(pjet(0, i))
! look for veto jet without y cuts
      if (i .gt. 2 .and. yjet(i) .gt. min(yjet(1), yjet(2)) .and. &
      & yjet(i) .lt. max(yjet(1), yjet(2)) .and. ijveto .eq. 0) ijveto = i

! now look for jets within the rapidity cuts
      if (dabs(yjet(i)) .lt. yjmax) then
        njety = njety + 1
        if (ij1y .eq. 0) then
          ij1y = i
        else if (ij2y .eq. 0) then
          ij2y = i
        else if (ij3y .eq. 0) then
          ij3y = i
        end if
! look for veto jet with y cuts
        if (ij3y .gt. 0 .and. &
        & yjet(i) .gt. min(yjet(ij1y), yjet(ij2y)) .and. &
        & yjet(i) .lt. max(yjet(ij1y), yjet(ij2y)) .and. ijvetoy .eq. 0) &
        & ijvetoy = i
      end if
    end do

! Nason-Oleari cuts (hep-ph/0911.5299)
    ptj_tag = 20d0
    deltay12 = 4.d0
    mj1j2min = 600d0

! this is the loop for w-o / w vbf cuts
    do i = 1, 2
      if (i .eq. 1) then
        ij1 = 1
        ij2 = 2
        ij3 = 3
      end if
      if (i .eq. 2) then
        njet = njety
        ijveto = ijvetoy
        ij1 = ij1y
        ij2 = ij2y
        ij3 = ij3y
      end if

! Load momenta
      do k = 0, 3
        pH(k) = p(k, 3)
        pj1(k) = pjet(k, ij1)
        pj2(k) = pjet(k, ij2)
        pj3(k) = pjet(k, ij3)
        pjj(k) = pjet(k, ij1) + pjet(k, ij2)
        pHj(k) = pjet(k, ij1) + pH(k)
        psyst(k) = pjet(k, ij1) + pjet(k, ij2) + pH(k)
        pjveto(k) = pjet(k, ijveto)
      end do

! Define observables
! Higgs
      ptH = getptv4(pH)
      etaH = getpseudorapv4(pH)
      yH = getrapidityv4(pH)
      njdble = dble(njet)
! At least one jet
      if (njet .ge. 1) then
        ptj1 = getptv4(pj1)
        etaj1 = getpseudorapv4(pj1)
        yj1 = getrapidityv4(pj1)
        ptHj = getptv4(pHj)
        etaHj = getpseudorapv4(pHj)
        yHj = getrapidityv4(pHj)
        DphiHj1 = getdelphiv4(pH, pj1)
        DRHj1 = getdrv4(pH, pj1)
      end if
! At least two jets
      if (njet .ge. 2) then
        ptj2 = getptv4(pj2)
        etaj2 = getpseudorapv4(pj2)
        yj2 = getrapidityv4(pj2)
        ptjj = getptv4(pjj)
        etajj = getpseudorapv4(pjj)
        yjj = getrapidityv4(pjj)
        ptsyst = getptv4(psyst)
        etasyst = getpseudorapv4(psyst)
        ysyst = getrapidityv4(psyst)
        DphiHj2 = getdelphiv4(pH, pj2)
        Dphij1j2 = getdelphiv4(pj1, pj2)
        DRHj2 = getdrv4(pH, pj2)
        DRj1j2 = getdrv4(pj1, pj2)
        mj1j2 = getinvmv4(pjj)
        Dyj1j2 = abs(yj1 - yj2)
      end if
! At least three jets
      if (njet .ge. 3) then
        ptj3 = getptv4(pj3)
        etaj3 = getpseudorapv4(pj3)
        yj3 = getrapidityv4(pj3)
        yj3rel = yj3 - (yj1 + yj2)/2d0
      end if
!
      chy1 = cosh(yj1)
      shy1 = sinh(yj1)
      chy1mo = chy1 - 1.d0
      chy2 = cosh(yj2)
      shy2 = sinh(yj2)
      chy2mo = chy2 - 1.d0

      call boostwdir2(chy1, shy1, chy1mo, xd, pj1, pj1boost)
      call boostwdir2(chy2, shy2, chy2mo, xd, pj2, pj2boost)
      ptrel_j1 = 0d0
      ptrel_j2 = 0d0

      pass_tag_cuts = njety .ge. 2 .and. &
      & ptj1 .ge. ptj_tag .and. &
      & ptj2 .ge. ptj_tag .and. &
      & abs(yj1 - yj2) .ge. deltay12 .and. &
      & yj1*yj2 .le. 0d0 .and. &
      & mj1j2 .ge. mj1j2min

      if (i .eq. 1) then
        flag = .true.
      end if

      if (i .eq. 2) then
        flag = pass_tag_cuts
      end if

      do j = 1, nQCD
        if (njet .ge. 1 .and. jet(j) .eq. 1) then
          call boostwdir2(chy1, shy1, chy1mo, xd, pQCD(0, j), ppboost(0, j))
          call getwedge(ppboost(0, j), pj1boost, prel_j1)
          ptrel_j1 = ptrel_j1 + getmod(prel_j1)/getmod(pj1boost)
        elseif (njet .ge. 2 .and. jet(j) .eq. 2) then
          call boostwdir2(chy2, shy2, chy2mo, xd, pQCD(0, j), ppboost(0, j))
          call getwedge(ppboost(0, j), pj2boost, prel_j2)
          ptrel_j2 = ptrel_j2 + getmod(prel_j2)/getmod(pj2boost)
        end if
      end do

      l = (i - 1)*54
      if (flag) then
        call HwU_fill(l + 1, 1d0, wgts)
        call HwU_fill(l + 2, ptH, wgts)
        call HwU_fill(l + 3, ptH, wgts)
        if (ptH .gt. 0d0) call HwU_fill(l + 4, log10(ptH), wgts)
        call HwU_fill(l + 5, etaH, wgts)
        call HwU_fill(l + 6, yH, wgts)
        call HwU_fill(l + 46, njdble, wgts)

        if (njet .ge. 1) then
          call HwU_fill(l + 7, ptj1, wgts)
          call HwU_fill(l + 8, ptj1, wgts)
          if (ptj1 .gt. 0d0) call HwU_fill(l + 9, log10(ptj1), wgts)
          call HwU_fill(l + 10, etaj1, wgts)
          call HwU_fill(l + 11, yj1, wgts)
          call HwU_fill(l + 22, ptHj, wgts)
          call HwU_fill(l + 23, ptHj, wgts)
          if (ptHj .gt. 0d0) call HwU_fill(l + 24, log10(ptHj), wgts)
          call HwU_fill(l + 25, etaHj, wgts)
          call HwU_fill(l + 26, yHj, wgts)
          call HwU_fill(l + 37, DphiHj1, wgts)
          call HwU_fill(l + 40, DRHj1, wgts)
          call HwU_fill(l + 47, ptrel_j1, wgts)
        end if

        if (njet .ge. 2) then
          call HwU_fill(l + 12, ptj2, wgts)
          call HwU_fill(l + 13, ptj2, wgts)
          if (ptj2 .gt. 0d0) call HwU_fill(l + 14, log10(ptj2), wgts)
          call HwU_fill(l + 15, etaj2, wgts)
          call HwU_fill(l + 16, yj2, wgts)
          call HwU_fill(l + 27, ptjj, wgts)
          call HwU_fill(l + 28, ptjj, wgts)
          if (ptjj .gt. 0d0) call HwU_fill(l + 29, log10(ptjj), wgts)
          call HwU_fill(l + 30, etajj, wgts)
          call HwU_fill(l + 31, yjj, wgts)
          call HwU_fill(l + 32, ptsyst, wgts)
          call HwU_fill(l + 33, ptsyst, wgts)
          if (ptsyst .gt. 0d0) call HwU_fill(l + 34, log10(ptsyst), wgts)
          call HwU_fill(l + 35, etasyst, wgts)
          call HwU_fill(l + 36, ysyst, wgts)
          call HwU_fill(l + 38, DphiHj2, wgts)
          call HwU_fill(l + 39, Dphij1j2, wgts)
          call HwU_fill(l + 41, DRHj2, wgts)
          call HwU_fill(l + 42, DRj1j2, wgts)
          call HwU_fill(l + 43, mj1j2, wgts)
          call HwU_fill(l + 44, Dyj1j2, wgts)
          call HwU_fill(l + 48, ptrel_j2, wgts)
        end if

        if (njet .ge. 3) then
          call HwU_fill(l + 17, ptj3, wgts)
          call HwU_fill(l + 18, ptj3, wgts)
          if (ptj3 .gt. 0d0) call HwU_fill(l + 19, log10(ptj3), wgts)
          call HwU_fill(l + 20, etaj3, wgts)
          call HwU_fill(l + 21, yj3, wgts)
          call HwU_fill(l + 45, yj3rel, wgts)
        end if
        if (ijveto .gt. 0) then
          pt_veto = getptv4(pjveto)
          if (.not. allocated(veto_weights)) then
            write (*, *) 'analysis_fill called before analysis_begin'
            stop 1
          end if
          veto_weights = wgts(1:size(veto_weights))
          do k = 1, nbinveto
            if (pt_veto .gt. (vetomin + (vetomax - vetomin)*dble(k - 1) &
            & /dble(nbinveto))) then
              call HwU_fill(l + 49, (vetomax - vetomin)*dble(k) &
              & /dble(nbinveto)*0.99d0, veto_weights)
            end if
          end do
          call HwU_fill(l + 50, pt_veto, wgts)
          call HwU_fill(l + 51, pt_veto, wgts)
          if (pt_veto .gt. 0d0) call HwU_fill(l + 52, dlog10(pt_veto), wgts)
          call HwU_fill(l + 53, getpseudorapv4(pjveto), wgts)
          call HwU_fill(l + 54, getrapidityv4(pjveto), wgts)
        end if
      end if
    end do

    return
  end subroutine analysis_fill

  function getrapidity(en, pl)
    implicit none
    double precision getrapidity, en, pl, tiny, xplus, xminus, y
    parameter(tiny=1.d-8)
!
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

  function getdr(ptx1, pty1, pl1, ptx2, pty2, pl2)
    implicit none
    double precision getdr, ptx1, pty1, pl1, ptx2, pty2, pl2, deta, dphi
!
    deta = getpseudorap(ptx1, pty1, pl1) - getpseudorap(ptx2, pty2, pl2)
    dphi = getdelphi(ptx1, pty1, ptx2, pty2)
    getdr = sqrt(dphi**2 + deta**2)
    return
  end function getdr

  function getptv4(p)
    implicit none
    double precision getptv4, p(0:3)
!
    getptv4 = sqrt(p(1)**2 + p(2)**2)
    return
  end function getptv4

  function getpseudorapv4(p)
    implicit none
    double precision getpseudorapv4, p(0:3)
!
    getpseudorapv4 = getpseudorap(p(1), p(2), p(3))
    return
  end function getpseudorapv4

  function getrapidityv4(p)
    implicit none
    double precision getrapidityv4, p(0:3)
!
    getrapidityv4 = getrapidity(p(0), p(3))
    return
  end function getrapidityv4

  function getdrv4(p1, p2)
    implicit none
    double precision getdrv4, p1(0:3), p2(0:3)
!
    getdrv4 = getdr(p1(1), p1(2), p1(3), p2(1), p2(2), p2(3))
    return
  end function getdrv4

  function getinvmv4(p)
    implicit none
    double precision getinvmv4, p(0:3)
!
    getinvmv4 = getinvm(p(0), p(1), p(2), p(3))
    return
  end function getinvmv4

  function getdelphiv4(p1, p2)
    implicit none
    double precision getdelphiv4, p1(0:3), p2(0:3)
!
    getdelphiv4 = getdelphi(p1(1), p1(2), &
    & p2(1), p2(2))
    return
  end function getdelphiv4

  function getmod(p)
    implicit none
    double precision p(0:3), getmod

    getmod = sqrt(p(1)**2 + p(2)**2 + p(3)**2)

    return
  end function getmod

  subroutine getwedge(p1, p2, pout)
    implicit none
    double precision p1(0:3), p2(0:3), pout(0:3)

    pout(1) = p1(2)*p2(3) - p1(3)*p2(2)
    pout(2) = p1(3)*p2(1) - p1(1)*p2(3)
    pout(3) = p1(1)*p2(2) - p1(2)*p2(1)
    pout(0) = 0d0

    return
  end subroutine getwedge

end module analysis_hwu_pp_hjj_module
