      subroutine HwU_write_file
      implicit none
      double precision xnorm
c     PineAPPL commons (this may not be the best place to put it)
      include "reweight_pineappl.inc"
      include "pineappl_common.inc"
      logical pineappl
      common /for_pineappl/ pineappl
      integer j
      if(pineappl)then
         do j=1,nh_obs
           appl_obs_num = j
           call APPL_term
         enddo
      endif
      open (unit=99,file='MADatNLO.HwU',status='unknown')
      xnorm=1d0
      call HwU_output(99,xnorm)
      close (99)
      return
      end
