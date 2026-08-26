module fks_qcd_splitting
  implicit none
  private

  public :: AP_reduced, AP_reduced_prime
  public :: Qterms_reduced_timelike, Qterms_reduced_spacelike

contains

  subroutine AP_reduced(col1, col2, t, z, gs, ap)
! Returns Altarelli-Parisi splitting function summed/averaged over helicities
! times prefactors such that |M_n+1|^2 = ap * |M_n|^2. This means
!    AP_reduced = (1-z) P_{S(part1,part2)->part1+part2}(z) * g^2/t
! Therefore, the labeling conventions for particle IDs are not as in FKS:
! part1 and part2 are the two particles emerging from the branching.
! part1 and part2 can be either gluon (8) or (anti-)quark (+-3). z is the
! fraction of the energy of part1 and t is the invariant mass of the mother.
    implicit none

    integer col1, col2
    double precision z, ap(2), t, gs

    double precision CA, TR, CF
    parameter(CA=3d0, TR=1d0/2d0, CF=4d0/3d0)
    ap(2) = 0d0

    if (col1 .eq. 8 .and. col2 .eq. 8) then
! g->gg splitting
      ap(1) = 2d0*CA*((1d0 - z)**2/z + z + z*(1d0 - z)**2)
      ap(2) = 0d0

    elseif (abs(col1) .eq. 3 .and. abs(col2) .eq. 3) then
! g->qqbar splitting
      ap(1) = TR*(z**2 + (1d0 - z)**2)*(1d0 - z)

    elseif (abs(col1) .eq. 3 .and. col2 .eq. 8) then
! q->qg splitting
      ap(1) = CF*(1d0 + z**2)

    elseif (col1 .eq. 8 .and. abs(col2) .eq. 3) then
! q->gq splitting
      ap(1) = CF*(1d0 + (1d0 - z)**2)*(1d0 - z)/z

    else
      write (*, *) 'Fatal Error in AP_reduced', col1, col2
      stop
    end if

    ap(1) = ap(1)*gs**2/t
    return
  end subroutine AP_reduced

  subroutine AP_reduced_prime(col1, col2, t, z, gs, apprime)
! Returns (1-z)*P^\prime * gS^2/t, with the same conventions as AP_reduced
    implicit none

    integer col1, col2
    double precision z, apprime(2), t, gs

    double precision TR, CF
    parameter(TR=1d0/2d0, CF=4d0/3d0)
    apprime(2) = 0d0
    if (col1 .eq. 8 .and. col2 .eq. 8) then
! g->gg splitting
      apprime(1) = 0d0
      apprime(2) = 0d0

    elseif (abs(col1) .eq. 3 .and. abs(col2) .eq. 3) then
! g->qqbar splitting
      apprime(1) = -2*TR*z*(1d0 - z)**2

    elseif (abs(col1) .eq. 3 .and. col2 .eq. 8) then
! q->qg splitting
      apprime(1) = -CF*(1d0 - z)**2

    elseif (col1 .eq. 8 .and. abs(col2) .eq. 3) then
! q->gq splitting
      apprime(1) = -CF*z*(1d0 - z)
    else
      write (*, *) 'Fatal error in AP_reduced_prime', col1, col2
      stop
    end if

    apprime(1) = apprime(1)*gs**2/t
    return
  end subroutine AP_reduced_prime

  subroutine Qterms_reduced_timelike(col1, col2, t, z, gs, Qterms)
! Eq's B.31 to B.34 of FKS paper, times (1-z)*g^2/t. The labeling
! conventions for particle IDs are the same as those in AP_reduced
    implicit none

    integer col1, col2
    double precision z, Qterms(2), t, gs

    double precision CA, TR
    parameter(CA=3d0, TR=1d0/2d0)
    Qterms(2) = 0d0
    if (col1 .eq. 8 .and. col2 .eq. 8) then
! g->gg splitting
      Qterms(1) = -4d0*CA*z*(1d0 - z)**2
      Qterms(2) = 0d0

    elseif (abs(col1) .eq. 3 .and. abs(col2) .eq. 3) then
! g->qqbar splitting
      Qterms(1) = 4d0*TR*z*(1d0 - z)**2

    elseif (abs(col1) .eq. 3 .and. col2 .eq. 8) then
! q->qg splitting
      Qterms(1) = 0d0

    elseif (col1 .eq. 8 .and. abs(col2) .eq. 3) then
! q->gq splitting
      Qterms(1) = 0d0
    else
      write (*, *) 'Fatal error in Qterms_reduced_timelike', col1, col2
      stop
    end if

    Qterms(1) = Qterms(1)*gs**2/t
    return
  end subroutine Qterms_reduced_timelike

  subroutine Qterms_reduced_spacelike(col1, col2, t, z, gs, Qterms)
! Eq's B.42 to B.45 of FKS paper, times (1-z)*gS^2/t. The labeling
! conventions for particle IDs are the same as those in AP_reduced.
! Thus, part1 has momentum fraction z, and it is the one off-shell
! (see (FKS.B.41))
    implicit none

    integer col1, col2
    double precision z, Qterms(2), t, gs

    double precision CA, CF
    parameter(CA=3d0, CF=4d0/3d0)
    Qterms(2) = 0d0
    if (col1 .eq. 8 .and. col2 .eq. 8) then
! g->gg splitting
      Qterms(1) = -4d0*CA*(1d0 - z)**2/z
      Qterms(2) = 0d0

    elseif (abs(col1) .eq. 3 .and. abs(col2) .eq. 3) then
! g->qqbar splitting
      Qterms(1) = 0d0

    elseif (abs(col1) .eq. 3 .and. col2 .eq. 8) then
! q->qg splitting
      Qterms(1) = 0d0

    elseif (col1 .eq. 8 .and. abs(col2) .eq. 3) then
! q->gq splitting
      Qterms(1) = -4d0*CF*(1d0 - z)**2/z
    else
      write (*, *) 'Fatal error in Qterms_reduced_spacelike', col1, col2
      stop
    end if

    Qterms(1) = Qterms(1)*gs**2/t
    return
  end subroutine Qterms_reduced_spacelike

end module fks_qcd_splitting
