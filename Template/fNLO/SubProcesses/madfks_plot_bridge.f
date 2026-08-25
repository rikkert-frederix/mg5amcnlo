c Process-generated declaration boundary for madfks_plot_module.

      subroutine initplot
      use madfks_plot_module, only: initplot_impl
      implicit none

      call initplot_impl()
      end


      subroutine topout
      use madfks_plot_module, only: topout_impl
      implicit none

      call topout_impl()
      end


      subroutine outfun(pp,ybst_til_tolab,www,ipdg,itype)
      use madfks_plot_module, only: outfun_impl
      implicit none
      include 'nexternal.inc'
      double precision pp(0:3,nexternal),ybst_til_tolab,www(*)
      integer ipdg(nexternal),itype
      double precision pmass(nexternal)
      common /to_mass/ pmass

      call outfun_impl(pp,ybst_til_tolab,www,ipdg,itype,pmass)
      end


      subroutine plot_pine_bridge(action,norm,ibody,itype,www)
      use run_state, only: pineappl
      implicit none
      integer action,ibody,itype
      double precision norm,www
      integer amp_pos_plot
      common /campposplot/ amp_pos_plot
      include 'pineappl_common.inc'

      if (.not.pineappl) return
      if (action.eq.1) then
         appl_Q2min=-1d0
         appl_Q2max=-1d0
         appl_xmin=-1d0
         appl_xmax=-1d0
         appl_nQ2=-1
         appl_Q2order=-1
         appl_nx=-1
         appl_xorder=-1
      elseif (action.eq.2) then
         appl_norm_histo=norm
      elseif (action.eq.3) then
         appl_itype=ibody
c Special treatment for collinear and soft-collinear counterterms.
         if (ibody.eq.2) then
            if (itype.eq.13) then
               appl_itype=4
            elseif (itype.eq.14) then
               appl_itype=5
            endif
         endif
         appl_amp_pos=amp_pos_plot
         appl_www_histo=www
      else
         write(*,*) 'Unknown madfks plot PineAPPL action:',action
         stop 1
      endif
      end
