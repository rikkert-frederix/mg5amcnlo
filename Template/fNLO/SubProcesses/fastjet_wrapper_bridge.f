      subroutine amcatnlo_fastjetppgenkt_timed(pqcd,nn,rfj,sycut,
     &     palg,pjet,njet,jet)
      use fastjet_timing_wrapper, only: fastjet_timed
      implicit none
      integer nn,njet,jet(nn)
      double precision pqcd(0:3,nn),rfj,sycut,palg,pjet(0:3,nn)

      call fastjet_timed(pqcd,nn,rfj,sycut,palg,pjet,njet,jet)
      return
      end


      subroutine amcatnlo_fastjetppgenkt_etamax_timed(pqcd,nn,rfj,
     &     sycut,etamax,palg,pjet,njet,jet)
      use fastjet_timing_wrapper, only: fastjet_etamax_timed
      implicit none
      integer nn,njet,jet(nn)
      double precision pqcd(0:3,nn),rfj,sycut,etamax,palg
      double precision pjet(0:3,nn)

      call fastjet_etamax_timed(pqcd,nn,rfj,sycut,etamax,palg,pjet,
     &     njet,jet)
      return
      end
