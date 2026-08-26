module fastjet_timing_wrapper
  implicit none
  private

  public :: fastjet_timed
  public :: fastjet_etamax_timed

  interface
    subroutine amcatnlo_fastjetppgenkt(pqcd, nn, rfj, sycut, palg, &
                                       pjet, njet, jet)
      implicit none
      integer, intent(in) :: nn
      double precision, intent(in) :: pqcd(0:3, nn)
      double precision, intent(in) :: rfj
      double precision, intent(in) :: sycut
      double precision, intent(in) :: palg
      double precision, intent(out) :: pjet(0:3, nn)
      integer, intent(out) :: njet
      integer, intent(out) :: jet(nn)
    end subroutine amcatnlo_fastjetppgenkt

    subroutine amcatnlo_fastjetppgenkt_etamax(pqcd, nn, rfj, sycut, &
                                              etamax, palg, pjet, njet, jet)
      implicit none
      integer, intent(in) :: nn
      double precision, intent(in) :: pqcd(0:3, nn)
      double precision, intent(in) :: rfj
      double precision, intent(in) :: sycut
      double precision, intent(in) :: etamax
      double precision, intent(in) :: palg
      double precision, intent(out) :: pjet(0:3, nn)
      integer, intent(out) :: njet
      integer, intent(out) :: jet(nn)
    end subroutine amcatnlo_fastjetppgenkt_etamax
  end interface

contains

  subroutine fastjet_timed(pqcd, nn, rfj, sycut, palg, &
                           pjet, njet, jet)
    implicit none
    integer, intent(in) :: nn
    double precision, intent(in) :: pqcd(0:3, nn)
    double precision, intent(in) :: rfj
    double precision, intent(in) :: sycut
    double precision, intent(in) :: palg
    double precision, intent(out) :: pjet(0:3, nn)
    integer, intent(out) :: njet
    integer, intent(out) :: jet(nn)
    call amcatnlo_fastjetppgenkt(pqcd, nn, rfj, sycut, palg, &
                                 pjet, njet, jet)
  end subroutine fastjet_timed

  subroutine fastjet_etamax_timed(pqcd, nn, rfj, sycut, &
                                  etamax, palg, pjet, &
                                  njet, jet)
    implicit none
    integer, intent(in) :: nn
    double precision, intent(in) :: pqcd(0:3, nn)
    double precision, intent(in) :: rfj
    double precision, intent(in) :: sycut
    double precision, intent(in) :: etamax
    double precision, intent(in) :: palg
    double precision, intent(out) :: pjet(0:3, nn)
    integer, intent(out) :: njet
    integer, intent(out) :: jet(nn)
    call amcatnlo_fastjetppgenkt_etamax(pqcd, nn, rfj, sycut, etamax, &
                                        palg, pjet, njet, jet)
  end subroutine fastjet_etamax_timed

end module fastjet_timing_wrapper
