c Process-generated declaration boundary for madfks_plot_module.

      subroutine initplot
      use madfks_plot_module, only: initplot_impl
      implicit none

      call initplot_impl()
      end
      subroutine plot_pine_bridge(action,norm,ibody,itype,www)
      use run_state, only: pineappl
      use fnlo_process_common, only: amp_pos_plot,appl_q2min,
     $     appl_q2max,appl_xmin,appl_xmax,appl_nq2,appl_q2order,
     $     appl_nx,appl_xorder,appl_norm_histo,appl_itype,
     $     appl_amp_pos,appl_www_histo
      implicit none
      integer action,ibody,itype
      double precision norm,www
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
