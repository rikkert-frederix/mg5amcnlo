c Legacy momentum-reshuffling ABI for generated subprocess callers.

      subroutine reshuffle_momenta(p,q,iresh,pdg_old,pdg_new,pass)
      use process_dimensions, only: nexternal
      use momentum_reshuffling, only: reshuffle_momenta_impl =>
     &     reshuffle_momenta
      implicit none
      double precision p(0:3,nexternal-1),q(0:3,nexternal-1)
      integer iresh(2),pdg_old(2),pdg_new(2)
      logical pass

      call reshuffle_momenta_impl(p,q,iresh,pdg_old,pdg_new,pass)
      return
      end
