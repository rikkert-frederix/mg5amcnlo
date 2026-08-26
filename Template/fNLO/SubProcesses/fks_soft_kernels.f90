module fks_soft_kernels
  use process_dimensions, only: nexternal
  use kin_functions_module, only: dot => dot_impl
  implicit none
  private

  public :: eikonal_reduced, eikonal_Ireg

contains

  subroutine eikonal_reduced(pp, m, n, i_fks, j_fks, xi_i_fks, &
                               y_ij_fks, p_i_fks_cnt, pmass, &
                               sqrtshat, eik)
!     Returns the eikonal factor
    implicit none
    double precision eik, pp(0:3, nexternal), xi_i_fks, y_ij_fks
    double precision, intent(in) :: p_i_fks_cnt(0:3, 0:2)
    double precision, intent(in) :: pmass(nexternal), sqrtshat
    double precision dotnm, dotni, dotmi, fact
    integer n, m, i_fks, j_fks, i
    integer softcol


    double precision phat_i_fks(0:3)

    double precision zero, tiny
    parameter(zero=0d0)
    parameter(tiny=1d-6)
! Define the reduced momentum for i_fks
    softcol = 0
    if (1d0 - y_ij_fks .lt. tiny) softcol = 2
    if (p_i_fks_cnt(0, softcol) .lt. 0d0) then
      if (xi_i_fks .eq. 0.d0) then
        write (*, *) 'Error #1 in eikonal_reduced', softcol, xi_i_fks, y_ij_fks
        stop
      end if
      if (pp(0, i_fks) .ne. 0.d0) then
        write (*, *) 'WARNING in eikonal_reduced: no cnt momenta', softcol, xi_i_fks, y_ij_fks
        do i = 0, 3
          phat_i_fks(i) = pp(i, i_fks)/xi_i_fks
        end do
      else
        write (*, *) 'Error #2 in eikonal_reduced', softcol, xi_i_fks, y_ij_fks
        stop
      end if
    else
      do i = 0, 3
        phat_i_fks(i) = p_i_fks_cnt(i, softcol)
      end do
    end if
! Calculate the eikonal factor
    dotnm = dot(pp(0, n), pp(0, m))
    if ((m .ne. j_fks .and. n .ne. j_fks) .or. pmass(j_fks) .ne. ZERO) then
      dotmi = dot(pp(0, m), phat_i_fks)
      dotni = dot(pp(0, n), phat_i_fks)
      fact = 1d0 - y_ij_fks
    elseif (m .eq. j_fks .and. n .ne. j_fks .and. pmass(j_fks) .eq. ZERO) then
      dotni = dot(pp(0, n), phat_i_fks)
      dotmi = sqrtshat/2d0*pp(0, j_fks)
      fact = 1d0
    elseif (m .ne. j_fks .and. n .eq. j_fks .and. pmass(j_fks) .eq. ZERO) then
      dotni = sqrtshat/2d0*pp(0, j_fks)
      dotmi = dot(pp(0, m), phat_i_fks)
      fact = 1d0
    else
      write (*, *) 'Error #3 in eikonal_reduced'
      stop
    end if

    eik = dotnm/(dotni*dotmi)*fact
    return
  end subroutine eikonal_reduced

  subroutine eikonal_Ireg(p, m, n, xicut_used, pmass, shat, &
                            qes2, abrv, eikIreg)
    implicit none
    double precision pi, pi2
    parameter(pi=3.1415926535897932385d0)
    parameter(pi2=pi**2)
    double precision p(0:3, nexternal), xicut_used, eikIreg
    double precision, intent(in) :: pmass(nexternal), shat, qes2
    character(len=4), intent(in) :: abrv
    integer m, n



    double precision Ei, Ej, kikj, rij, tmp, xmj, betaj, betai
    double precision xmi2, xmj2, vij, xi0, alij, tHVvl, tHVv
    double precision arg1, arg2, arg3, arg4, xi1a

    tmp = 0.d0
    if (pmass(m) .eq. 0.d0 .and. pmass(n) .eq. 0.d0) then
      if (m .eq. n) then
        write (*, *) 'Error #2 in eikonal_Ireg', m, n
        stop
      end if
      Ei = p(0, n)
      Ej = p(0, m)
      kikj = dot(p(0, n), p(0, m))
      rij = kikj/(2*Ei*Ej)
      if (abs(rij - 1.d0) .gt. 1.d-6) then
        if (abrv .ne. 'virt') then
! 1+2+3+4
          tmp = 1d0/2d0*dlog(xicut_used**2*shat/QES2)**2 &
                + dlog(xicut_used**2*shat/QES2)*dlog(rij) &
                - ddilog(rij) + 1d0/2d0*dlog(rij)**2 &
                - dlog(1 - rij)*dlog(rij)
        else
          write (*, *) 'Error #11 in eikonal_Ireg', abrv
          stop
        end if
      else
        if (abrv .ne. 'virt') then
! 1+2+3+4
          tmp = 1d0/2d0*dlog(xicut_used**2*shat/QES2)**2 - pi2/6.d0
        else
          write (*, *) 'Error #12 in eikonal_Ireg', abrv
          stop
        end if
      end if
    elseif ((pmass(m) .ne. 0.d0 .and. pmass(n) .eq. 0.d0) .or. (pmass(m) .eq. 0.d0 .and. pmass(n) .ne. 0.d0)) then
      if (m .eq. n) then
        write (*, *) 'Error #3 in eikonal_Ireg', m, n
        stop
      end if
      if (pmass(m) .ne. 0.d0 .and. pmass(n) .eq. 0.d0) then
        Ei = p(0, n)
        Ej = p(0, m)
        xmj = pmass(m)
        betaj = sqrt(1 - xmj**2/Ej**2)
      else
        Ei = p(0, m)
        Ej = p(0, n)
        xmj = pmass(n)
        betaj = sqrt(1 - xmj**2/Ej**2)
      end if
      kikj = dot(p(0, n), p(0, m))
      rij = kikj/(2*Ei*Ej)

      if (abrv .ne. 'virt') then
! 1+2+3+4
        tmp = dlog(xicut_used) &
              *(dlog(xicut_used*shat/QES2) + 2*dlog(kikj/(xmj*Ei))) &
              - ddilog(1 - (1 + betaj)/(2*rij)) &
              + ddilog(1 - 2*rij/(1 - betaj)) &
              + 1/2.d0*log(2*rij/(1 - betaj))**2 &
              + dlog(shat/QES2)*dlog(kikj/(xmj*Ei)) - pi2/12.d0 &
              + 1/4.d0*dlog(shat/QES2)**2 &
              - 1/4.d0*dlog((1 + betaj)/(1 - betaj))**2
      else
        write (*, *) 'Error #13 in eikonal_Ireg', abrv
        stop
      end if
    elseif (pmass(m) .ne. 0.d0 .and. pmass(n) .ne. 0.d0) then
      if (n .eq. m) then
        Ei = p(0, n)
        betai = sqrt(1 - pmass(n)**2/Ei**2)
        if (abrv .ne. 'virt') then
! 1+2+3+4
          if (betai .gt. 1d-6) then
            tmp = dlog(xicut_used**2*shat/QES2) - 1/betai*dlog((1 + betai)/(1 - betai))
          else
            tmp = dlog(xicut_used**2*shat/QES2) - 2d0*(1d0 + betai**2/3d0 + betai**4/5d0)
          end if
        else
          write (*, *) 'Error #14 in eikonal_Ireg', abrv
          stop
        end if
      else
        Ei = p(0, n)
        Ej = p(0, m)
        betai = sqrt(1 - pmass(n)**2/Ei**2)
        betaj = sqrt(1 - pmass(m)**2/Ej**2)
        xmi2 = pmass(n)**2
        xmj2 = pmass(m)**2
        kikj = dot(p(0, n), p(0, m))
        vij = sqrt(1 - xmi2*xmj2/kikj**2)
        alij = kikj*(1 + vij)/xmi2
        tHVvl = (alij**2*xmi2 - xmj2)/2.d0
        tHVv = tHVvl/(alij*Ei - Ej)
        arg1 = alij*Ei
        arg2 = arg1*betai
        arg3 = Ej
        arg4 = arg3*betaj
        if (vij .lt. 1d0) then
          xi0 = 1/vij*log((1 + vij)/(1 - vij))
        else
          xi0 = dlog(4d0*kikj**2/(xmi2*xmj2))
        end if
!          xi0=1/vij*log((1+vij)/(1-vij))
        xi1a = kikj**2*(1 + vij)/xmi2*(xj1a(arg1, arg2, tHVv, tHVvl) - xj1a(arg3, arg4, tHVv, tHVvl))

        if (abrv .ne. 'virt') then
! 1+2+3+4
          tmp = 1/2.d0*xi0*dlog(xicut_used**2*shat/QES2) + 1/2.d0*xi1a
        else
          write (*, *) 'Error #15 in eikonal_Ireg', abrv
          stop
        end if
      end if
    else
      write (*, *) 'Error #4 in eikonal_Ireg', m, n, pmass(m), pmass(n)
      stop
    end if
    eikIreg = tmp
    return
  end subroutine eikonal_Ireg

  function xj1a(x, y, tHVv, tHVvl)
    implicit none
    double precision xj1a, x, y, tHVv, tHVvl
!
    xj1a = 1/(2*tHVvl)*(dlog((x - y)/(x + y))**2 + 4*ddilog(1 - (x + y)/tHVv) + 4*ddilog(1 - (x - y)/tHVv))
    return
  end function xj1a

  function DDILOG(X)
!
! $Id: imp64.inc,v 1.1.1.1 1996/04/01 15:02:59 mclareni Exp $
!
! $Log: imp64.inc,v $
! Revision 1.1.1.1  1996/04/01 15:02:59  mclareni
! Mathlib gen
!
!
! imp64.inc
!
    implicit none
    integer i
    double precision ddilog, x, y, s, a, t, h, alfa, b0, b1, b2, c(0:19)
    double precision z1, hf, pi, pi3, pi6, pi12
    parameter(Z1=1, HF=Z1/2)
    parameter(PI=3.14159265358979324d0)
    parameter(PI3=PI**2/3, PI6=PI**2/6, PI12=PI**2/12)
    data C(0)/0.42996693560813697d0/
    data C(1)/0.40975987533077105d0/
    data C(2)/-0.01858843665014592d0/
    data C(3)/0.00145751084062268d0/
    data C(4)/-0.00014304184442340d0/
    data C(5)/0.00001588415541880d0/
    data C(6)/-0.00000190784959387d0/
    data C(7)/0.00000024195180854d0/
    data C(8)/-0.00000003193341274d0/
    data C(9)/0.00000000434545063d0/
    data C(10)/-0.00000000060578480d0/
    data C(11)/0.00000000008612098d0/
    data C(12)/-0.00000000001244332d0/
    data C(13)/0.00000000000182256d0/
    data C(14)/-0.00000000000027007d0/
    data C(15)/0.00000000000004042d0/
    data C(16)/-0.00000000000000610d0/
    data C(17)/0.00000000000000093d0/
    data C(18)/-0.00000000000000014d0/
    data C(19)/+0.00000000000000002d0/
    if (X .eq. 1) then
      H = PI6
    elseif (X .eq. -1) then
      H = -PI12
    else
      T = -X
      if (T .le. -2) then
        Y = -1/(1 + T)
        S = 1
        A = -PI3 + HF*(log(-T)**2 - log(1 + 1/T)**2)
      elseif (T .lt. -1) then
        Y = -1 - T
        S = -1
        A = log(-T)
        A = -PI6 + A*(A + log(1 + 1/T))
      else if (T .le. -HF) then
        Y = -(1 + T)/T
        S = 1
        A = log(-T)
        A = -PI6 + A*(-HF*A + log(1 + T))
      else if (T .lt. 0) then
        Y = -T/(1 + T)
        S = -1
        A = HF*log(1 + T)**2
      else if (T .le. 1) then
        Y = T
        S = 1
        A = 0
      else
        Y = 1/T
        S = -1
        A = PI6 + HF*log(T)**2
      end if
      H = Y + Y - 1
      ALFA = H + H
      B1 = 0
      B2 = 0
      do I = 19, 0, -1
        B0 = C(I) + ALFA*B1 - B2
        B2 = B1
        B1 = B0
      end do
      H = -(S*(B0 - H*B2) + A)
    end if
    DDILOG = H
    return
  end function DDILOG

end module fks_soft_kernels
