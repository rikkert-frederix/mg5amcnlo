subroutine sborn_sf(p_born,m,n,wgt)
          implicit none
          include "nexternal.inc"
          include 'coupl.inc'
          double precision p_born(0:3,nexternal-1), wgt
          integer nsqampso
          parameter (nsqampso = %(nsqorders)d)

C Return the colour-linked Born and include the gs**2 factor.
          double precision wgt_col
          integer m,n
      include 'orders.inc'
      double precision amp_split_soft(amp_split_size)
      common /to_amp_split_soft/amp_split_soft

%(iflines_col)s
      wgt = wgt_col * g**2
      amp_split_soft(1:amp_split_size) =
     $     dble(amp_split_cnt(1:amp_split_size,1,qcd_pos)) * g**2

      return
      end
