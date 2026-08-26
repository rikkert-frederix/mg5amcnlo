c Process-generated PineAPPL declarations are isolated here so that the
c histogram implementation in HwU.f90 has no generated INCLUDE files.

      subroutine hwu_pineappl_inithist
      use run_state, only: pineappl
      use fnlo_process_common, only: appl_obs_nbins
      implicit none

      if (pineappl) appl_obs_nbins = 0
      end


      subroutine hwu_pineappl_book(label,title_l,nbin_l,xmin,xmax)
      use run_state, only: pineappl
      use fnlo_process_common, only: appl_obs_nbins,appl_obs_bins,
     $     appl_obs_min,appl_obs_max,nh_obs,ih_obs
      implicit none
      integer label,nbin_l,i
      character(len=*) title_l
      double precision xmin,xmax,del

      if (pineappl.and.index(title_l,'Born').eq.0) then
         if (appl_obs_nbins.eq.0) then
            appl_obs_nbins=nbin_l
            del=(xmax-xmin)/nbin_l
            do i=0,appl_obs_nbins
               appl_obs_bins(i)=xmin+i*del
            enddo
         endif
         appl_obs_min=appl_obs_bins(0)
         appl_obs_max=appl_obs_bins(appl_obs_nbins)
         if (abs(appl_obs_max-xmax).gt.0.00000001d0) then
            write (*,*) 'PineAPPL Histogram: ',
     1           'Change of the upper limit:',xmax,'-->',appl_obs_max
         endif
         call APPL_init
         nh_obs=nh_obs+1
         ih_obs(nh_obs)=label
         appl_obs_nbins=0
      endif
      end


      subroutine hwu_pineappl_fill(label,x)
      use run_state, only: pineappl
      use fnlo_process_common, only: nh_obs,ih_obs,appl_obs_num,
     $     appl_obs_histo
      implicit none
      integer label,j
      double precision x

      if (pineappl) then
         do j=1,nh_obs
            if (label.eq.ih_obs(j)) then
               appl_obs_num=j
               appl_obs_histo=x
               call APPL_fill
            endif
         enddo
      endif
      end
      subroutine HwU_add_points
      use HwU_module,
     &     only: module_HwU_add_points => HwU_add_points
      implicit none
      call module_HwU_add_points
      end


      subroutine HwU_accum_iter(inclde,nPSpoints,values)
      use HwU_module,
     &     only: module_HwU_accum_iter => HwU_accum_iter
      implicit none
      logical inclde
      integer nPSpoints
      double precision values(2)
      call module_HwU_accum_iter(inclde,nPSpoints,values)
      end
      subroutine HwU_output(unit,xnorm)
      use HwU_module, only: module_HwU_output => HwU_output
      implicit none
      integer unit
      double precision xnorm
      call module_HwU_output(unit,xnorm)
      end
