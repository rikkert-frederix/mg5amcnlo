CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                          C
C                  HwU: Histograms with Uncertainties                      C
C                 By Rikkert Frederix, 12-2014--05-2017                    C
C                                                                          C
C     Book, fill and write out histograms. Particularly suited for NLO     C
C     computations with correlations between points (ie. event and         C
C     counter-event) and multiple weights for given points (e.g. scale     C
C     and PDF uncertainties through reweighting).                          C
C                                                                          C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      module HwU_wgts_info_len
         use iso_c_binding
         integer, parameter :: wgts_info_len=80 
         contains 
             integer function get_wgts_info_len() bind(c,name="get_wgts_info_len")
                 get_wgts_info_len = wgts_info_len
                 return
             end function get_wgts_info_len

      end module HwU_wgts_info_len


c The module contains effectively the common block with allocatable
c variables (something not possible in old fortran version)
      module HwU_variables
         use HwU_wgts_info_len
         implicit none
         integer :: max_plots,max_points,max_bins,nwgts,np
         integer :: error_estimation=3
         logical :: profile_errors=.false.
         logical, allocatable :: booked(:)
         integer, allocatable :: nbin(:),histi(:,:),hist_iter(:,:)
     $        ,p_bin(:),p_label(:)
         character(len=50), allocatable :: title(:)
         character(len=wgts_info_len), allocatable :: wgts_info(:)
         double precision, allocatable :: histy(:,:,:),histy_acc(:,:,:)
     $        ,histy2(:,:,:),histy_err(:,:,:),histxl(:,:),histxm(:,:)
     $        ,step(:),p_wgts(:,:)
         save
      end module HwU_variables

      
c To be called once at the start of each run. Initialises the packages
c and sets the number of weights that need to be included for each point.
      subroutine HwU_inithist(nweights,wgt_info)
      use HwU_variables
      implicit none
      integer i,nweights
      character*(*) wgt_info(*)

c     PineAPPL commons
      include "reweight_pineappl.inc"
      include "pineappl_common.inc"
      logical pineappl
      common /for_pineappl/ pineappl

C     Initialize the number of bins of the aMCfast grids
      if(pineappl) appl_obs_nbins = 0

      call HwU_deallocate_all
      max_plots=0
      max_points=0
      max_bins=0
      np=0
c     Number of weights associated to each point. The first weight is the
c     'central value'. Other weights are allowed to be non-zero when the
c     central one vanishes, as is needed for correlated migration weights.
      nwgts=nweights
      allocate(wgts_info(nwgts))
      profile_errors=.false.
      do i=1,nwgts
         wgts_info(i)=wgt_info(i)
         if (index(adjustl(wgts_info(i)),'FKSdamp ').eq.1)
     $        profile_errors=.true.
      enddo
      return
      end

      subroutine set_error_estimation(input)
c Error estimation
c     input 0: Simply sum all PS points; error is estimated by assuming
c        Binomial statistics: use for unweighted events. (Assumes that
c        all events have equal weight (up to a sign)).
c     input 1: Simply sum all PS points; error is estimated by the
c        variance between the separate iterations. This assumes that
c        each iteration has the same number of PS points (This is the
c        topdrawer default in MG5_aMC).
c     input 2: Sum PS points for a given iteration and error estimate by
c       square root of the sum of the squares. Perform a weighted average
c       iteration-by-iteration
c     input 3: Same as input 2, but weighted average is same as from MINT
      use HwU_variables
      implicit none
      integer input
      if (input.ge.0 .and. input.le.3) then
         error_estimation=input
      else
         write (*,*) 'unknown error estimation',input
         stop 1
      endif
      return
      end
            
c Book the histograms at the start of the run. Give a 'label' (an
c integer) that identifies the plot when filling it and a title
c ('title_l') for each plot. Also the number of bins ('nbin_l') and the
c plot range (from 'xmin' to 'xmax') should be given.
      subroutine HwU_book(label,title_l,nbin_l,xmin,xmax)
      use HwU_variables
      implicit none
      integer label,nbin_l,i,j
      character*(*) title_l
      double precision xmin,xmax
c     PineAPPL commons
      include "reweight_pineappl.inc"
      include "pineappl_common.inc"
      logical pineappl
      common /for_pineappl/ pineappl
      double precision del

c     Initialize the grids only if the switch "pineappl" is set to True
c     and if the title  does not contain the word "Born". 
      if(pineappl.and.index(title_l,"Born").eq.0)then
c     Observable parameters
c     Compute number of bins and edges only if they have not been given by the user.
         if(appl_obs_nbins.eq.0)then
            appl_obs_nbins = nbin_l
            ! bin width
            del = (xmax - xmin) / nbin_l
c     compute bin edges
            do i=0,appl_obs_nbins
               appl_obs_bins(i) = xmin + i * del
            enddo
         endif
         appl_obs_min = appl_obs_bins(0)
         appl_obs_max = appl_obs_bins(appl_obs_nbins)
         if(abs(appl_obs_max-xmax).gt.0.00000001d0)then
            write(*,*) 'PineAPPL Histogram: ', 
     1                 'Change of the upper limit:',xmax,'-->',
     2                  appl_obs_max
         endif
c     Initialize PineAPPL routines
         call APPL_init
c     Keep track of the position of this histogram
         nh_obs = nh_obs + 1
         ih_obs(nh_obs) = label
c     Reset number of bins to zero
         appl_obs_nbins = 0
      endif

c     Allocate space for new histograms if needed      
      call HwU_allocate_histo(label,nbin_l)
c     Setup the histogram
      booked(label)=.true.
      title(label)=title_l
      nbin(label)=nbin_l
c     Compute the bin width
      step(label)=(xmax-xmin)/dble(nbin(label))
      do i=1,nbin(label)
c     Compute the lower and upper bin edges
         histxl(label,i)=xmin+step(label)*dble(i-1)
         histxm(label,i)=xmin+step(label)*dble(i)
c     Set all the bins to zero.
         do j=1,nwgts
            histy(j,label,i)=0d0
            histy_acc(j,label,i)=0d0
            histy2(j,label,i)=0d0
            histy_err(j,label,i)=0d0
         enddo
         histi(label,i)=0
         hist_iter(label,i)=0
      enddo
      return
      end
      
c Fill the histograms identified with 'label' with a point at x with
c weights giving by 'wgts'. The dimension of the 'wgts' array should
c always be as large as the number of weights ('nweights') specified in
c the 'HwU_inithist' subroutine. That means that each point should have
c the same number of weights.
      subroutine HwU_fill(label,x,wgts)
      use HwU_variables
      implicit none
      integer label,i,j,bin
      logical nonzero_weight
      double precision x, wgts(*)

c     PineAPPL commons
      include "reweight_pineappl.inc"
      include "pineappl_common.inc"
      logical pineappl
      common /for_pineappl/ pineappl
      if(pineappl)then
         do j=1,nh_obs
            if(label.eq.ih_obs(j))then
               appl_obs_num   = j
               appl_obs_histo = x
c     Fill the reference PineAPPL histograms
c     Fill the PineAPPL files
               call APPL_fill
            endif
         enddo
      endif

c     Preserve the legacy central-weight condition unless profile errors
c     are requested. A migration counterevent can have zero central weight
c     and a non-zero damping-profile weight.
      if (.not.profile_errors .and. wgts(1).eq.0d0) return
      nonzero_weight=.false.
      do j=1,nwgts
         if (wgts(j).ne.0d0) then
            nonzero_weight=.true.
            exit
         endif
      enddo
      if (.not.nonzero_weight) return
c     Check if point is within plotting range
      if (x.lt.histxl(label,1) .or.
     $     x.gt.histxm(label,nbin(label))) return
c     Compute the bin to fill
      bin=int((x-histxl(label,1))/step(label)) +1
c     To prevent numerical inaccuracy, double check that bin is
c     reasonable
      if (bin.lt.1 .or. bin.gt.nbin(label)) return
c     Check if we already had another (correlated) point in this bin. If
c     so, add the current weight to that point.
      do i=1,np
         if ( p_label(i).eq.label .and. p_bin(i).eq.bin) then
            do j=1,nwgts
               p_wgts(j,i)=p_wgts(j,i)+wgts(j)
            enddo
            return
         endif
      enddo
c     If a new bin, add it to the list of points
      np=np+1
      call HwU_allocate_p
      p_label(np)=label
      p_bin(np)=bin
      do j=1,nwgts
         p_wgts(j,np)=wgts(j)
      enddo
      return
      end
      
c Call after all correlated contributions for a give phase-space
c point. I.e., every time you get a new set of random numbers from
c MINT/VEGAS. It adds the current list of points to the histograms and
c adds the square of every correlated weight to compute its statistical
c uncertainty. Squaring here, rather than in HwU_fill, is essential: event
c and counterevent contributions have already been summed in p_wgts.
      subroutine HwU_add_points
      use HwU_variables
      implicit none
      integer i,j
      do i=1,np
         do j=1,nwgts
            histy(j,p_label(i),p_bin(i))=
     $           histy(j,p_label(i),p_bin(i))+p_wgts(j,i)
            if (profile_errors) histy2(j,p_label(i),p_bin(i))=
     $           histy2(j,p_label(i),p_bin(i))+p_wgts(j,i)**2
         enddo
         if (.not.profile_errors) histy2(1,p_label(i),p_bin(i))=
     $        histy2(1,p_label(i),p_bin(i))+p_wgts(1,i)**2
         histi(p_label(i),p_bin(i))=histi(p_label(i),p_bin(i))+1
      enddo
      np=0
      return
      end


c *****For plotting during fixed order f(N)LO runs*****
c Call after every iteration. This adds the histograms of the current
c iteration ('histy') to the accumulated results ('histy_acc'), with the
c uncertainty estimate given in 'histy_err' and empties the arrays for
c the current iteration so that they can be filled with the next
c iteration.
      subroutine HwU_accum_iter(inclde,nPSpoints,values)
      use HwU_variables
      implicit none
      logical inclde
      integer nPSpoints,label,i,j
      double precision nPSinv,etot,niter,y_squared,values(2)
      data niter /0d0/
      nPSinv = 1d0/dble(nPSpoints)
      if (inclde) niter = niter+1d0
      do label=1,max_plots
         if (.not. booked(label)) cycle
         if (inclde) then
            call accumulate_results(label,nPSinv,niter,values)
         endif
c     Reset the histo of the current iteration to zero so that they are
c     ready for the next iteration.
         do i=1,nbin(label)
            do j=1,nwgts
               histy(j,label,i)=0d0
               histy2(j,label,i)=0d0
            enddo
            histi(label,i)=0
         enddo
      enddo
      return
      end

c *****For plotting unweighted events****
c Overwrites the accumulated results ('histy_acc') with the current
c results ('histy') and provides an uncertainty estimate in
c 'histy_err'. This means that during the plotting of events, at
c intermediate stages this function can be called (together with
c HwU_output) to write intermediate plots to disk.
      subroutine finalize_histograms(nPSpoints)
      use HwU_variables
      implicit none
      integer label,nPSpoints,i,j
      double precision nPSinv,niter,dummy(2)
      nPSinv=1d0/dble(nPSpoints)
      niter=1d0
      do label=1,max_plots
         if (.not. booked(label)) cycle
         do i=1,nbin(label)
c     Set all the bins to zero.
            do j=1,nwgts
               histy_acc(j,label,i)=0d0
               histy_err(j,label,i)=0d0
            enddo
            hist_iter(label,i)=0
         enddo
         call accumulate_results(label,nPSinv,niter,dummy)
      enddo
      return
      end

c This adds the histograms of the current iteration ('histy') to the
c accumulated results ('histy_acc'), with one uncertainty estimate in
c 'histy_err' for each weight. When error_estimation==2 each weight uses
c its own bin uncertainty for the iteration average. With
c error_estimation==3 the means retain the common MINT iteration weights,
c while each uncertainty is propagated separately. With
c error_estimation==1 the iteration variance is evaluated weight by weight;
c error_estimation==0 uses the corresponding sum of squared weights.
      subroutine accumulate_results(label,nPSinv,niter,values)
      use HwU_variables
      implicit none
      integer label,i,j
      double precision nPSinv,niter,values(2),a1,a2,denom
      double precision,allocatable :: vtot(:),etot(:),y_squared(:)
      if (.not. allocated(vtot)) then
         allocate(vtot(nwgts))
         allocate(etot(nwgts))
         allocate(y_squared(nwgts))
      endif
      if (.not.profile_errors) then
         call accumulate_results_legacy(label,nPSinv,niter,values)
         return
      endif
      if (error_estimation.eq.2) then
c     Use the weighted average bin-by-bin. This is not really justified
c     for fNLO computations, because for bins with low statistics, the
c     variance computation cannot be trusted.
         do i=1,nbin(label)
c     Skip bin if no entries
            if (histi(label,i).eq.0) cycle
c     Divide weights by the number of PS points. This means that this is
c     now normalised to the total cross section in that bin
            do j=1,nwgts
               vtot(j)=histy(j,label,i)*nPSinv
            enddo
            do j=1,nwgts
c     Error estimation of the current bin
               etot(j)=sqrt(abs(histy2(j,label,i)*nPSinv-vtot(j)**2)
     $              *nPSinv)
c     Include "Bessel's correction" to have a corrected (even though
c     still biased) estimator of the standard deviation.
               if (histi(label,i).gt.1) then
                  etot(j)=etot(j)*sqrt(dble(histi(label,i))
     &                 /(dble(histi(label,i))-1.5d0))
               else
                  etot(j)=abs(vtot(j))*10d0 ! make one-entry error large
               endif
c     Initialise this bin on its first contributing iteration.
               if (hist_iter(label,i).eq.0) then
                  histy_acc(j,label,i)=vtot(j)
                  histy_err(j,label,i)=etot(j)
               elseif (histy_err(j,label,i).ne.0d0 .and.
     $                 etot(j).ne.0d0) then
c     Add the results of the current iteration to the accumulated results
                  denom=1d0/histy_err(j,label,i)+1d0/etot(j)
                  histy_acc(j,label,i)=(histy_acc(j,label,i)
     &                 /histy_err(j,label,i)+vtot(j)/etot(j))/denom
                  histy_err(j,label,i)=1d0/sqrt(1d0/
     $                 histy_err(j,label,i)**2+1d0/etot(j)**2)
               else
c     A zero sample variance is common for sparse profile weights and must
c     not be used as an infinite inverse-variance weight. Fall back to an
c     iteration-count average for this update.
                  a1=1d0/dble(hist_iter(label,i)+1)
                  a2=1d0-a1
                  histy_acc(j,label,i)=a2*histy_acc(j,label,i)
     $                 +a1*vtot(j)
                  histy_err(j,label,i)=sqrt((a2*
     $                 histy_err(j,label,i))**2+(a1*etot(j))**2)
               endif
            enddo
            hist_iter(label,i)=hist_iter(label,i)+1
         enddo
      elseif (error_estimation.eq.3) then
         do i=1,nbin(label)
c     Skip bin if no entries
            if (histi(label,i).eq.0) cycle
c     Divide weights by the number of PS points. This means that this is
c     now normalised to the total cross section in that bin
            do j=1,nwgts
               vtot(j)=histy(j,label,i)*nPSinv
            enddo
            do j=1,nwgts
c     Error estimation of the current bin
               etot(j)=sqrt(abs(histy2(j,label,i)*nPSinv-vtot(j)**2)
     $              *nPSinv)
c     Include "Bessel's correction" to have a corrected (even though
c     still biased) estimator of the standard deviation.
               if (histi(label,i).gt.1) then
                  etot(j)=etot(j)*sqrt(dble(histi(label,i))
     &                 /(dble(histi(label,i))-1.5d0))
               else
                  etot(j)=abs(vtot(j))*10d0 ! make one-entry error large
               endif
            enddo
c     Initialise this bin on its first contributing iteration.
            if (hist_iter(label,i).eq.0) then
               do j=1,nwgts
                  histy_acc(j,label,i)=vtot(j)
                  histy_err(j,label,i)=etot(j)
               enddo
            elseif (values(1).ne.0d0 .and. values(2).ne.0d0) then
c     Add the results of the current iteration to the accumulated results
               do j=1,nwgts
                  histy_acc(j,label,i)=(histy_acc(j,label,i)
     &                 /values(2)+vtot(j)/values(1))/(1d0
     &                 /values(2) + 1d0/values(1))
               enddo
               a1=((1d0/values(1))/((1d0/values(1))+1d0/values(2)))**2
               a2=((1d0/values(2))/((1d0/values(1))+1d0/values(2)))**2
               do j=1,nwgts
                  histy_err(j,label,i)=sqrt(a2*
     $                 histy_err(j,label,i)**2+a1*etot(j)**2)
               enddo
            else
c     Guard the degenerate MINT-error case with an iteration-count average.
               a1=1d0/dble(hist_iter(label,i)+1)
               a2=1d0-a1
               do j=1,nwgts
                  histy_acc(j,label,i)=a2*histy_acc(j,label,i)
     $                 +a1*vtot(j)
                  histy_err(j,label,i)=sqrt((a2*
     $                 histy_err(j,label,i))**2+(a1*etot(j))**2)
               enddo
            endif
            hist_iter(label,i)=hist_iter(label,i)+1
         enddo 
      elseif(error_estimation.eq.1) then
c     simply sum the weights in the bins
         do i=1,nbin(label)
            do j=1,nwgts
               if (niter.ne.1d0) y_squared(j)=((niter-1)*
     $              histy_err(j,label,i))**2+(niter-1)*
     $              histy_acc(j,label,i)**2
               vtot(j)=histy(j,label,i)*nPSinv
               histy_acc(j,label,i)=(histy_acc(j,label,i)*(niter-1d0)
     &              +vtot(j))/niter
c     base the error on the variance in the results per iteration. For a
c     small number of iterations, this underestimates the actual
c     uncertainty.
               if (niter.eq.1d0) then
                  histy_err(j,label,i)=0d0
               else
                  histy_err(j,label,i)=sqrt(max(0d0,((y_squared(j)
     $                 +vtot(j)**2)/niter-histy_acc(j,label,i)**2)
     $                 /niter))
               endif
            enddo
         enddo
      elseif(error_estimation.eq.0) then
c     simply sum the weights in the bins
         do i=1,nbin(label)
            do j=1,nwgts
               vtot(j)=histy(j,label,i)*nPSinv
               histy_acc(j,label,i)=(histy_acc(j,label,i)*(niter-1d0)
     &              +vtot(j))/niter
            enddo
c     Error estimation of the current bin using Poisson statistics,
c     making sure we normalise with the number of iterations.
            if (histi(label,i).ne.0) then
               do j=1,nwgts
                  etot(j)=sqrt(histy2(j,label,i))*nPSinv
                  histy_err(j,label,i)=sqrt(((niter-1d0)*
     $                 histy_err(j,label,i))**2+etot(j)**2)/niter
               enddo
            else
               do j=1,nwgts
                  histy_err(j,label,i)=(niter-1d0)/niter*
     $                 histy_err(j,label,i)
               enddo
            endif
         enddo
      endif
      end


c Preserve the established HwU estimator exactly when no FKS migration
c profiles are present.  This keeps zero-profile runs numerically and
c format compatible with the 3.x branch.
      subroutine accumulate_results_legacy(label,nPSinv,niter,values)
      use HwU_variables
      implicit none
      integer label,i,j
      double precision nPSinv,etot,niter,y_squared,values(2),a1,a2
      double precision,allocatable :: vtot(:)
      if (.not.allocated(vtot)) allocate(vtot(nwgts))
      if (error_estimation.eq.2) then
         do i=1,nbin(label)
            if (histi(label,i).eq.0) cycle
            do j=1,nwgts
               vtot(j)=histy(j,label,i)*nPSinv
            enddo
            etot=sqrt(abs(histy2(1,label,i)*nPSinv-vtot(1)**2)
     $           *nPSinv)
            if (histi(label,i).gt.1) then
               etot=etot*sqrt(dble(histi(label,i))
     $              /(dble(histi(label,i))-1.5d0))
            else
               etot=abs(vtot(1))*10d0
            endif
            if (histy_err(1,label,i).eq.0d0) then
               do j=1,nwgts
                  histy_acc(j,label,i)=vtot(j)
               enddo
               histy_err(1,label,i)=etot
            else
               do j=1,nwgts
                  histy_acc(j,label,i)=(histy_acc(j,label,i)
     $                 /histy_err(1,label,i)+vtot(j)/etot)/(1d0
     $                 /histy_err(1,label,i)+1d0/etot)
               enddo
               histy_err(1,label,i)=1d0/sqrt(1d0/
     $              histy_err(1,label,i)**2+1d0/etot**2)
            endif
         enddo
      elseif (error_estimation.eq.3) then
         do i=1,nbin(label)
            if (histi(label,i).eq.0) cycle
            do j=1,nwgts
               vtot(j)=histy(j,label,i)*nPSinv
            enddo
            etot=sqrt(abs(histy2(1,label,i)*nPSinv-vtot(1)**2)
     $           *nPSinv)
            if (histi(label,i).gt.1) then
               etot=etot*sqrt(dble(histi(label,i))
     $              /(dble(histi(label,i))-1.5d0))
            else
               etot=abs(vtot(1))*10d0
            endif
            if (histy_err(1,label,i).eq.0d0) then
               do j=1,nwgts
                  histy_acc(j,label,i)=vtot(j)
               enddo
               histy_err(1,label,i)=etot
            else
               do j=1,nwgts
                  histy_acc(j,label,i)=(histy_acc(j,label,i)
     $                 /values(2)+vtot(j)/values(1))/(1d0
     $                 /values(2)+1d0/values(1))
               enddo
               a1=((1d0/values(1))/((1d0/values(1))
     $              +1d0/values(2)))**2
               a2=((1d0/values(2))/((1d0/values(1))
     $              +1d0/values(2)))**2
               histy_err(1,label,i)=sqrt(a2*
     $              histy_err(1,label,i)**2+a1*etot**2)
            endif
         enddo
      elseif (error_estimation.eq.1) then
         do i=1,nbin(label)
            if (histi(label,i).eq.0 .and.
     $           histy_acc(1,label,i).eq.0d0) cycle
            if (niter.ne.1d0) y_squared=((niter-1d0)*
     $           histy_err(1,label,i))**2+(niter-1d0)*
     $           histy_acc(1,label,i)**2
            do j=1,nwgts
               vtot(j)=histy(j,label,i)*nPSinv
               histy_acc(j,label,i)=(histy_acc(j,label,i)
     $              *(niter-1d0)+vtot(j))/niter
            enddo
            if (niter.eq.1d0) then
               histy_err(1,label,i)=0d0
            else
               histy_err(1,label,i)=sqrt(((y_squared+vtot(1)**2)
     $              /niter-histy_acc(1,label,i)**2)/niter)
            endif
         enddo
      elseif (error_estimation.eq.0) then
         do i=1,nbin(label)
            if (histi(label,i).eq.0 .and.
     $           histy_acc(1,label,i).eq.0d0) cycle
            do j=1,nwgts
               vtot(j)=histy(j,label,i)*nPSinv
               histy_acc(j,label,i)=(histy_acc(j,label,i)
     $              *(niter-1d0)+vtot(j))/niter
            enddo
            if (histi(label,i).ne.0) then
               etot=sqrt(histy2(1,label,i))*nPSinv
               histy_err(1,label,i)=sqrt(((niter-1d0)*
     $              histy_err(1,label,i))**2+etot**2)/niter
            else
               histy_err(1,label,i)=(niter-1d0)/niter*
     $              histy_err(1,label,i)
            endif
         enddo
      endif
      return
      end
      
c Write the histograms to disk at the end of the run, multiplying the
c output by 'xnorm'
      subroutine HwU_output(unit,xnorm)
      use HwU_variables
      implicit none
      integer unit,i,j,label
      integer max_length
      character(len=:), allocatable :: buffer
      character*4 str_nbin
      double precision xnorm
      if (.not. allocated(buffer))
     &     allocate(character(len=(2*nwgts+2)*17) :: buffer)
c     column info: x_min, x_max, y (central value), dy, {extra
c     weights}, {dy for each migration-profile weight}.
      write (unit,'(a)', advance='no') '##& xmin'
      write (unit,'(a)', advance='no') ' & xmax'
      write (unit,'(a)', advance='no') ' & '//trim(adjustl(wgts_info(1)))
      write (unit,'(a)', advance='no') ' & dy'
      do j=2,nwgts
         write (unit,'(a)', advance='no') ' & '//trim(adjustl(wgts_info(j)))
      enddo
      do j=2,nwgts
         if (index(adjustl(wgts_info(j)),'FKSdamp ').eq.1)
     $        write (unit,'(a)', advance='no') ' & dy['//
     $        trim(adjustl(wgts_info(j)))//']'
      enddo
      write (unit,'(a)') ''
      write (unit,'(a)') ''
      do label=1,max_plots
         if (.not. booked(label)) cycle
c     title
c        For some weird reason, it is no possible to include directly
c        the integer in two line below with the format '(12a,i4.4,3a,a,2a)'
c        this is why I have to do it in two steps instead.
         write (str_nbin,'(i4.4)') nbin(label)
         write (unit, '(12a,4a,3a,a,2a)') '<histogram> ',str_nbin,' " ',trim(title(label)), ' "'
c         write (*, '(12a,4a,3a,a,2a)') '<histogram> ',str_nbin,' " ',trim(title(label)), ' "'
c     data
         do i=1,nbin(label)
           write (unit,'(2x,e14.7)', advance='no') histxl(label,i)
           write (unit,'(2x,e14.7)', advance='no') histxm(label,i)
           write (unit,'(2x,e14.7)', advance='no') histy_acc(1,label,i)*xnorm
           write (unit,'(2x,e14.7)', advance='no')
     $          histy_err(1,label,i)*xnorm
           do j=2,nwgts
              write (unit, '(2x,e14.7)', advance='no') histy_acc(j,label,i)*xnorm
           enddo
           do j=2,nwgts
              if (index(adjustl(wgts_info(j)),'FKSdamp ').eq.1)
     $             write (unit, '(2x,e14.7)', advance='no')
     $             histy_err(j,label,i)*xnorm
           enddo
	   write(unit, *) ''
         enddo
c     2 empty lines after each plot
         write (unit,'(12a)') '<\histogram>'
         write (unit,'(a)') ''
         write (unit,'(a)') ''         
      enddo
      deallocate(buffer)
      return
      end

c Clean all the allocatable variables:
      subroutine HwU_deallocate_all
      use HwU_variables
      implicit none
      if (allocated(wgts_info)) deallocate(wgts_info)
      if (allocated(booked)) deallocate(booked)
      if (allocated(title)) deallocate(title)
      if (allocated(nbin)) deallocate(nbin)
      if (allocated(step)) deallocate(step)
      if (allocated(histxl)) deallocate(histxl)
      if (allocated(histxm)) deallocate(histxm)
      if (allocated(histy)) deallocate(histy)
      if (allocated(histy_acc)) deallocate(histy_acc)
      if (allocated(histi)) deallocate(histi)
      if (allocated(hist_iter)) deallocate(hist_iter)
      if (allocated(histy2)) deallocate(histy2)
      if (allocated(histy_err)) deallocate(histy_err)
      if (allocated(p_bin)) deallocate(p_bin)
      if (allocated(p_label)) deallocate(p_label)
      if (allocated(p_wgts)) deallocate(p_wgts)
      return
      end


      subroutine HwU_allocate_p
      use HwU_variables
      implicit none
      integer,allocatable :: itemp1(:)
      double precision, allocatable :: temp2(:,:)
      if (.not. allocated(p_bin)) then
         allocate(p_bin(max_plots))
         allocate(p_label(max_plots))
         allocate(p_wgts(nwgts,max_plots))
         max_points=max_plots
      else
         if (np.gt.max_points) then
c p_bin
            allocate(itemp1(np+max_plots))
            itemp1(1:max_points)=p_bin
            call move_alloc(itemp1,p_bin)
      
c p_label
            allocate(itemp1(np+max_plots))
            itemp1(1:max_points)=p_label
            call move_alloc(itemp1,p_label)
c p_wgts
            allocate(temp2(nwgts,np+max_plots))
            temp2(1:nwgts,1:max_points)=p_wgts
            call move_alloc(temp2,p_wgts)
            max_points=np+max_plots
         endif
      endif
      return
      end

      subroutine HwU_allocate_histo(label,nbin_l)
      use HwU_variables
      implicit none
      logical,allocatable :: ltemp(:)
      integer,allocatable :: itemp1(:),itemp2(:,:)
      character(len=50),allocatable :: ctemp(:)
      double precision, allocatable :: temp1(:),temp2(:,:),temp3(:,:,:)
      integer label,i,nbin_l,label_max,nbin_max
      logical debug
      parameter (debug=.true.)
c Check if variables are already allocated. If not, simply allocate a
c single histogram
      if (.not. allocated(booked)) then
         allocate(booked(1))
         booked(1)=.false.
         allocate(title(1))
         allocate(nbin(1))
         allocate(step(1))
         allocate(histxl(1,nbin_l))
         allocate(histxm(1,nbin_l))
         allocate(histy(nwgts,1,nbin_l))
         allocate(histy_acc(nwgts,1,nbin_l))
         allocate(histi(1,nbin_l))
         allocate(hist_iter(1,nbin_l))
         allocate(histy2(nwgts,1,nbin_l))
         allocate(histy_err(nwgts,1,nbin_l))
         max_plots=1
         max_bins=nbin_l
      endif
c If current label is greater than the plots already allocated, increase
c the size of the allocated arrays. This is kind of slow, but shouldn't
c really matter since it's only done at the start of a run.
      if (label.gt.max_plots .or. nbin_l.gt.max_bins) then
         label_max=max(label,max_plots)
         nbin_max=max(nbin_l,max_bins)
c booked
         allocate(ltemp(label_max))
         ltemp(1:max_plots)=booked
         call move_alloc(ltemp,booked)
         do i=max_plots+1,label_max
            booked(i)=.false. ! histos have not yet been setup
         enddo
c title
         allocate(ctemp(label_max))
         ctemp(1:max_plots)=title
         call move_alloc(ctemp,title)
c nbin
         allocate(itemp1(label_max))
         itemp1(1:max_plots)=nbin
         call move_alloc(itemp1,nbin)
c step
         allocate(temp1(label_max))
         temp1(1:max_plots)=step
         call move_alloc(temp1,step)
c histxl
         allocate(temp2(label_max,nbin_max))
         temp2(1:max_plots,1:max_bins)=histxl
         call move_alloc(temp2,histxl)
c histxm
         allocate(temp2(label_max,nbin_max))
         temp2(1:max_plots,1:max_bins)=histxm
         call move_alloc(temp2,histxm)
c histy
         allocate(temp3(nwgts,label_max,nbin_max))
         temp3(1:nwgts,1:max_plots,1:max_bins)=histy
         call move_alloc(temp3,histy)
c histy_acc
         allocate(temp3(nwgts,label_max,nbin_max))
         temp3(1:nwgts,1:max_plots,1:max_bins)=histy_acc
         call move_alloc(temp3,histy_acc)
c histi
         allocate(itemp2(label_max,nbin_max))
         itemp2(1:max_plots,1:max_bins)=histi
         call move_alloc(itemp2,histi)
c hist_iter
         allocate(itemp2(label_max,nbin_max))
         itemp2(1:max_plots,1:max_bins)=hist_iter
         call move_alloc(itemp2,hist_iter)
c histy2
         allocate(temp3(nwgts,label_max,nbin_max))
         temp3(1:nwgts,1:max_plots,1:max_bins)=histy2
         call move_alloc(temp3,histy2)
c histy_err
         allocate(temp3(nwgts,label_max,nbin_max))
         temp3(1:nwgts,1:max_plots,1:max_bins)=histy_err
         call move_alloc(temp3,histy_err)
c Update maximums
         max_plots=label_max
         max_bins=nbin_max
      elseif (booked(label)) then
         write (*,*) 'ERROR in HwU.f: histogram already booked',label
         stop
      endif
      return
      end
      
      
c dummy subroutine
      subroutine accum(idummy)
      integer idummy
      end
c dummy subroutine
      subroutine addfil(string)
      character*(*) string
      end
