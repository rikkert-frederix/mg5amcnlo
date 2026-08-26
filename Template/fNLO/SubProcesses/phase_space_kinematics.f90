module phase_space_kinematics
  use process_dimensions, only: nexternal, nincoming, max_particles, &
                                max_branch
  use kin_functions_module, only: dot => dot_impl, rho => rho_impl
  implicit none
  private

  public :: rotate_invar, getaziangles
  public :: phspncheck_born, phspncheck_nocms
  public :: phase_space_lambda

contains

  subroutine rotate_invar(pin, pout, cth, sth, cphi, sphi)
! Given the four momentum pin, returns the four momentum pout (in the same
! Lorentz frame) by performing a three-rotation of an angle theta
! (cos(theta)=cth) around the y axis, followed by a three-rotation of an
! angle phi (cos(phi)=cphi) along the z axis. The components of pin
! and pout are given along these axes
    implicit none
    double precision cth, sth, cphi, sphi, pin(0:3), pout(0:3)
    double precision q1, q2, q3
!
    q1 = pin(1)
    q2 = pin(2)
    q3 = pin(3)
    pout(1) = q1*cphi*cth - q2*sphi + q3*cphi*sth
    pout(2) = q1*sphi*cth + q2*cphi + q3*sphi*sth
    pout(3) = -q1*sth + q3*cth
    pout(0) = pin(0)
    return
  end subroutine rotate_invar

  subroutine getaziangles(p, cphi, sphi)
    implicit none
    double precision p(0:3), cphi, sphi
    double precision xlength, cth, sth
!
    xlength = rho(p)
    if (xlength .ne. 0.d0) then
      cth = p(3)/xlength
      sth = sqrt(1 - cth**2)
      if (sth .ne. 0.d0) then
        cphi = p(1)/(xlength*sth)
        sphi = p(2)/(xlength*sth)
      else
        cphi = 1.d0
        sphi = 0.d0
      end if
    else
      cphi = 1.d0
      sphi = 0.d0
    end if
    return
  end subroutine getaziangles

  subroutine phspncheck_born(ecm, xmass, xmom, pass)
! Checks four-momentum conservation.
! WARNING: works only in the partonic c.m. frame
    implicit none
    double precision ecm, xmass(nexternal - 1), xmom(0:3, nexternal - 1)
    double precision tiny, xm, xsum(0:3), xsuma(0:3), xrat(0:3), ptmp(0:3)
    parameter(tiny=5.d-3)
    integer jflag, npart, i, j, jj
    logical pass
!
    pass = .true.
    jflag = 0
    npart = nexternal - 1
    do i = 0, 3
      xsum(i) = 0.d0
      xsuma(i) = 0.d0
      do j = nincoming + 1, npart
        xsum(i) = xsum(i) + xmom(i, j)
        xsuma(i) = xsuma(i) + abs(xmom(i, j))
      end do
      if (i .eq. 0) xsum(i) = xsum(i) - ecm
      if (xsuma(i) .lt. 1.d0) then
        xrat(i) = abs(xsum(i))
      else
        xrat(i) = abs(xsum(i))/xsuma(i)
      end if
      if (xrat(i) .gt. tiny .and. jflag .eq. 0) then
        write (*, *) 'Momentum is not conserved'
        write (*, *) 'i=', i
        do j = 1, npart
          write (*, '(4(d14.8,1x))') (xmom(jj, j), jj=0, 3)
        end do
        jflag = 1
      end if
    end do
    if (jflag .eq. 1) then
      write (*, '(4(d14.8,1x))') (xsum(jj), jj=0, 3)
      write (*, '(4(d14.8,1x))') (xrat(jj), jj=0, 3)
      pass = .false.
      return
    end if
!
    do j = 1, npart
      do i = 0, 3
        ptmp(i) = xmom(i, j)
      end do
      xm = xlen4(ptmp)
      if (abs(xm - xmass(j))/ptmp(0) .gt. tiny .and. abs(xm - xmass(j)) .gt. tiny) then
        write (*, *) 'Mass shell violation'
        write (*, *) 'j=', j
        write (*, *) 'mass=', xmass(j)
        write (*, *) 'mass computed=', xm
        write (*, '(4(d14.8,1x))') (xmom(jj, j), jj=0, 3)
        pass = .false.
        return
      end if
    end do
    return
  end subroutine phspncheck_born

  subroutine phspncheck_nocms(npart, ecm, xmass, xmom, pass)
! Checks four-momentum conservation. Derived from phspncheck;
! works in any frame
    implicit none
    integer npart
    double precision ecm, xmass(-max_branch:max_particles), xmom(0:3, nexternal)
    double precision tiny, vtiny, xm, den, ecmtmp, xsum(0:3), xsuma(0:3), xrat(0:3), ptmp(0:3)
    parameter(tiny=5.d-3)
    parameter(vtiny=1.d-6)
    integer jflag, i, j, jj
    logical pass
!
    pass = .true.
    jflag = 0
    do i = 0, 3
      if (nincoming .eq. 2) then
        xsum(i) = -xmom(i, 1) - xmom(i, 2)
        xsuma(i) = abs(xmom(i, 1)) + abs(xmom(i, 2))
      elseif (nincoming .eq. 1) then
        xsum(i) = -xmom(i, 1)
        xsuma(i) = abs(xmom(i, 1))
      end if
      do j = nincoming + 1, npart
        xsum(i) = xsum(i) + xmom(i, j)
        xsuma(i) = xsuma(i) + abs(xmom(i, j))
      end do
      if (xsuma(i) .lt. 1.d0) then
        xrat(i) = abs(xsum(i))
      else
        xrat(i) = abs(xsum(i))/xsuma(i)
      end if
      if (xrat(i) .gt. tiny .and. jflag .eq. 0) then
        write (*, *) 'Momentum is not conserved (/nocms/)'
        write (*, *) 'i=', i
        do j = 1, npart
          write (*, '(i2,1x,4(d14.8,1x))') j, (xmom(jj, j), jj=0, 3)
        end do
        jflag = 1
      end if
    end do
    if (jflag .eq. 1) then
      write (*, '(a3,1x,4(d14.8,1x))') 'sum', (xsum(jj), jj=0, 3)
      write (*, '(a3,1x,4(d14.8,1x))') 'rat', (xrat(jj), jj=0, 3)
      pass = .false.
      return
    end if
!
    do j = 1, npart
      do i = 0, 3
        ptmp(i) = xmom(i, j)
      end do
      xm = xlen4(ptmp)
      if (ptmp(0) .ge. 1.d0) then
        den = ptmp(0)
      else
        den = 1.d0
      end if
      if (abs(xm - xmass(j))/den .gt. tiny .and. abs(xm - xmass(j)) .gt. tiny) then
        write (*, *) 'Mass shell violation (/nocms/)'
        write (*, *) 'j=', j
        write (*, *) 'mass=', xmass(j)
        write (*, *) 'mass computed=', xm
        write (*, '(4(d14.8,1x))') (xmom(jj, j), jj=0, 3)
        pass = .false.
        return
      end if
    end do
!
    if (nincoming .eq. 2) then
      ecmtmp = sqrt(2d0*dot(xmom(0, 1), xmom(0, 2)))
    elseif (nincoming .eq. 1) then
      ecmtmp = xmom(0, 1)
    end if
    if (abs(ecm - ecmtmp) .gt. vtiny) then
      write (*, *) 'Inconsistent shat (/nocms/)'
      write (*, *) 'ecm given=   ', ecm
      write (*, *) 'ecm computed=', ecmtmp
      write (*, '(4(d14.8,1x))') (xmom(jj, 1), jj=0, 3)
      write (*, '(4(d14.8,1x))') (xmom(jj, 2), jj=0, 3)
      pass = .false.
      return
    end if

    return
  end subroutine phspncheck_nocms

  function xlen4(v)
    implicit none
    double precision xlen4, tmp, v(0:3)
!
    tmp = v(0)**2 - v(1)**2 - v(2)**2 - v(3)**2
    xlen4 = sign(1.d0, tmp)*sqrt(abs(tmp))
    return
  end function xlen4

  double precision function phase_space_lambda(S, MA2, MB2)
    implicit none
!****************************************************************************
!     THIS IS THE phase_space_lambda FUNCTION FROM VERNONS BOOK COLLIDER PHYSICS P 662
!     MA2 AND MB2 ARE THE MASS SQUARED OF THE FINAL STATE PARTICLES
!     2-D PHASE SPACE = .5*PI*SQRT(1.,MA2/S^2,MB2/S^2)*(D(OMEGA)/4PI)
!****************************************************************************
    double precision MA2, MB2, S, tiny, tmp, rat
    parameter(tiny=1.d-8)
!
    tmp = S**2 + MA2**2 + MB2**2 - 2d0*S*MA2 - 2d0*MA2*MB2 - 2d0*S*MB2
    if (tmp .le. 0.d0) then
      if (ma2 .lt. 0.d0 .or. mb2 .lt. 0.d0) then
        write (6, *) 'Error #1 in function phase_space_lambda:', s, ma2, mb2
        stop
      end if
      rat = 1 - (sqrt(ma2) + sqrt(mb2))/s
      if (rat .gt. -tiny) then
        tmp = 0.d0
      else
        write (6, *) 'Error #2 in function phase_space_lambda:', s, ma2, mb2, rat
      end if
    end if
    phase_space_lambda = tmp
    return
  end function phase_space_lambda

end module phase_space_kinematics
