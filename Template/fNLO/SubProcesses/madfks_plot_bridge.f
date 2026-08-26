c Process-generated declaration boundary for madfks_plot_module.

      subroutine initplot
      use madfks_plot_module, only: initplot_impl
      implicit none

      call initplot_impl()
      end
