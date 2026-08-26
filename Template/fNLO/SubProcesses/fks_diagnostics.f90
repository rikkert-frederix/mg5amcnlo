module fks_diagnostics
  use process_dimensions, only: nexternal, nincoming
  use fnlo_process_common, only: i_momcmp_count, xratmax
  implicit none
  private

  public :: xmom_compare, xprintout, checkres2

contains

  subroutine xmom_compare(i_fks, j_fks, jac_cnt, p1_cnt, pmass, pass)
    implicit none
    double precision, intent(in) :: p1_cnt(0:3, nexternal, 0:2)
    double precision, intent(in) :: jac_cnt(0:2), pmass(nexternal)
    integer i_fks, j_fks
    integer izero, ione, itwo, iunit, isum
    logical verbose, pass, pass0
    parameter(izero=0)
    parameter(ione=1)
    parameter(itwo=2)
    parameter(iunit=6)
    parameter(verbose=.false.)
!
    isum = 0
    if (jac_cnt(0) .gt. 0.d0) isum = isum + 1
    if (jac_cnt(1) .gt. 0.d0) isum = isum + 2
    if (jac_cnt(2) .gt. 0.d0) isum = isum + 4
    pass = .true.
!
    if (isum .eq. 0 .or. isum .eq. 1 .or. isum .eq. 2 .or. isum .eq. 4) then
! Nothing to be done: 0 or 1 configurations computed
      if (verbose) write (iunit, *) 'none'
    elseif (isum .eq. 3 .or. isum .eq. 5 .or. isum .eq. 7) then
! Soft is taken as reference
      if (isum .eq. 7) then
        if (verbose) then
          write (iunit, *) 'all'
          write (iunit, *) '    '
          write (iunit, *) 'C/S'
        end if
        call xmcompare(verbose, pass0, ione, izero, i_fks, j_fks, p1_cnt, pmass)
        pass = pass .and. pass0
        if (verbose) then
          write (iunit, *) '    '
          write (iunit, *) 'SC/S'
        end if
        call xmcompare(verbose, pass0, itwo, izero, i_fks, j_fks, p1_cnt, pmass)
        pass = pass .and. pass0
      elseif (isum .eq. 3) then
        if (verbose) then
          write (iunit, *) 'C+S'
          write (iunit, *) '    '
          write (iunit, *) 'C/S'
        end if
        call xmcompare(verbose, pass0, ione, izero, i_fks, j_fks, p1_cnt, pmass)
        pass = pass .and. pass0
      elseif (isum .eq. 5) then
        if (verbose) then
          write (iunit, *) 'SC+S'
          write (iunit, *) '    '
          write (iunit, *) 'SC/S'
        end if
        call xmcompare(verbose, pass0, itwo, izero, i_fks, j_fks, p1_cnt, pmass)
        pass = pass .and. pass0
      end if
    elseif (isum .eq. 6) then
! Collinear is taken as reference
      if (verbose) then
        write (iunit, *) 'SC+C'
        write (iunit, *) '    '
        write (iunit, *) 'SC/C'
      end if
      call xmcompare(verbose, pass0, itwo, ione, i_fks, j_fks, p1_cnt, pmass)
      pass = pass .and. pass0
    else
      write (6, *) 'Fatal error in xmom_compare', isum
      stop
    end if
    if (.not. pass) i_momcmp_count = i_momcmp_count + 1
!
    return
  end subroutine xmom_compare

  subroutine xmcompare(verbose, pass0, inum, iden, i_fks, j_fks, p1_cnt, pmass)
    implicit none
    double precision, intent(in) :: p1_cnt(0:3, nexternal, 0:2)
    double precision, intent(in) :: pmass(nexternal)
    logical verbose, pass0
    integer inum, iden, i_fks, j_fks, iunit, ipart, i, j, k
    double precision tiny, vtiny, xnum, xden, xrat
    parameter(iunit=6)
    parameter(tiny=1.d-4)
    parameter(vtiny=1.d-10)
!
    pass0 = .true.
    do ipart = 1, nexternal
      do i = 0, 3
        xnum = p1_cnt(i, ipart, inum)
        xden = p1_cnt(i, ipart, iden)
        if (verbose) then
          if (i .eq. 0) then
            write (iunit, *) ' '
            write (iunit, *) 'part=', ipart
          end if
          call xprintout(iunit, xnum, xden)
        else
          if (ipart .ne. i_fks .and. ipart .ne. j_fks) then
            if (xden .ne. 0.d0) then
              xrat = abs(1 - xnum/xden)
            else
              xrat = abs(xnum)
            end if
            if (abs(xnum) .eq. 0d0 .and. abs(xden) .le. vtiny) xrat = 0d0
! The following line solves some problem as well, but before putting
! it as the standard, one should think a bit about it
            if (abs(xnum) .le. vtiny .and. abs(xden) .le. vtiny) xrat = 0d0
            if (xrat .gt. tiny .and. (pmass(ipart) .eq. 0d0 .or. xnum/pmass(ipart) .gt. vtiny)) then
              write (*, *) 'Kinematics of counterevents'
              write (*, *) inum, iden
              write (*, *) 'is different. Particle:', ipart
              write (*, *) xrat, xnum, xden
              do j = 1, nexternal
                write (*, *) j, (p1_cnt(k, j, inum), k=0, 3)
              end do
              do j = 1, nexternal
                write (*, *) j, (p1_cnt(k, j, iden), k=0, 3)
              end do
              xratmax = max(xratmax, xrat)
              pass0 = .false.
            end if
          end if
        end if
      end do
    end do
    do i = 0, 3
      if (j_fks .gt. nincoming) then
        xnum = p1_cnt(i, i_fks, inum) + p1_cnt(i, j_fks, inum)
        xden = p1_cnt(i, i_fks, iden) + p1_cnt(i, j_fks, iden)
      else
        xnum = -p1_cnt(i, i_fks, inum) + p1_cnt(i, j_fks, inum)
        xden = -p1_cnt(i, i_fks, iden) + p1_cnt(i, j_fks, iden)
      end if
      if (verbose) then
        if (i .eq. 0) then
          write (iunit, *) ' '
          write (iunit, *) 'part=i+j'
        end if
        call xprintout(iunit, xnum, xden)
      else
        if (xden .ne. 0.d0) then
          xrat = abs(1 - xnum/xden)
        else
          xrat = abs(xnum)
        end if
        if (xrat .gt. tiny) then
          write (*, *) 'Kinematics of counterevents'
          write (*, *) inum, iden
          write (*, *) 'is different. Particle i+j'
          xratmax = max(xratmax, xrat)
          pass0 = .false.
        end if
      end if
    end do
    return
  end subroutine xmcompare

  subroutine xprintout(iunit, xv, xlim)
    implicit none
    integer iunit
    double precision xv, xlim
!
    if (abs(xlim) .gt. 1.d-30) then
      write (iunit, *) xv/xlim, xv, xlim
    else
      write (iunit, *) xv, xlim
    end if
    return
  end subroutine xprintout

! The following has been derived with minor modifications from the
! analogous routine written for VBF
  subroutine checkres2(xsecvc, xseclvc, wgt, wgtl, xp, lxp, iflag, imax, iev, i_fks, j_fks, iret)
!     same as checkres, but also limits are arrays.
    implicit none
    double precision xsecvc(15), xseclvc(15), wgt(15), wgtl(15), lxp(0:3, nexternal + 1), xp(15, 0:3, nexternal + 1)
    double precision ckc(15), rckc(15), rat
    integer iflag, imax, iev, i_fks, j_fks, iret, ithrs, istop, iwrite, i, k, l, imin, icount
    parameter(ithrs=3)
    parameter(istop=0)
    parameter(iwrite=1)
!
    if (imax .gt. 15) then
      write (6, *) 'Error in checkres: imax is too large', imax
      stop
    end if
    do i = 1, imax
      if (xseclvc(i) .eq. 0.d0) then
        ckc(i) = abs(xsecvc(i))
      else
        ckc(i) = abs(xsecvc(i)/xseclvc(i) - 1.d0)
      end if
    end do
    if (iflag .eq. 0) then
      rat = 4.d0
    elseif (iflag .eq. 1) then
      rat = 2.d0
    else
      write (6, *) 'Error in checkres: iflag=', iflag
      write (6, *) ' Must be 0 for soft, 1 for collinear'
      stop
    end if
!
    i = 1
    do while (ckc(i) .gt. 0.1d0 .and. xseclvc(i) .ne. 0d0)
      i = i + 1
    end do
    imin = i
    do i = imin, imax - 1
      if (ckc(i + 1) .ne. 0.d0) then
        rckc(i) = ckc(i)/ckc(i + 1)
      else
        rckc(i) = 1.d8
      end if
    end do
    icount = 0
    i = imin
    do while (icount .lt. ithrs .and. i .lt. imax)
      if (rckc(i) .gt. rat) then
        icount = icount + 1
      else
        icount = 0
      end if
      i = i + 1
    end do
!
    iret = 0
    if (icount .ne. ithrs) then
      iret = 1
      if (istop .eq. 1) then
        write (6, *) 'Test failed', iflag
        write (6, *) 'Event #', iev
        stop
      end if
      if (iwrite .eq. 1) then
        write (77, *) '    '
        if (iflag .eq. 0) then
          write (77, *) 'Soft #', iev
        elseif (iflag .eq. 1) then
          write (77, *) 'Collinear #', iev
        end if
        write (77, *) 'ME*wgt:'
        do i = 1, imax
          call xprintout(77, xsecvc(i), xseclvc(i))
        end do
        write (77, *) 'wgt:'
        do i = 1, imax
          call xprintout(77, wgt(i), wgtl(i))
        end do
!
        write (78, *) '    '
        if (iflag .eq. 0) then
          write (78, *) 'Soft #', iev
        elseif (iflag .eq. 1) then
          write (78, *) 'Collinear #', iev
        end if
        do k = 1, nexternal
          write (78, *) ''
          write (78, *) 'part:', k
          do l = 0, 3
            write (78, *) 'comp:', l
            do i = 1, imax
              call xprintout(78, xp(i, l, k), lxp(l, k))
            end do
          end do
        end do
        if (iflag .eq. 0) then
          write (78, *) ''
          write (78, *) 'part: i_fks reduced'
          do l = 0, 3
            write (78, *) 'comp:', l
            do i = 1, imax
              call xprintout(78, xp(i, l, nexternal + 1), lxp(l, nexternal + 1))
            end do
          end do
          write (78, *) ''
          write (78, *) 'part: i_fks full/reduced'
          do l = 0, 3
            write (78, *) 'comp:', l
            do i = 1, imax
              call xprintout(78, xp(i, l, i_fks), xp(i, l, nexternal + 1))
            end do
          end do
        elseif (iflag .eq. 1) then
          write (78, *) ''
          write (78, *) 'part: i_fks+j_fks'
          do l = 0, 3
            write (78, *) 'comp:', l
            do i = 1, imax
              call xprintout(78, xp(i, l, i_fks) + xp(i, l, j_fks), lxp(l, i_fks) + lxp(l, j_fks))
            end do
          end do
        end if
      end if
    end if
    return
  end subroutine checkres2

end module fks_diagnostics
