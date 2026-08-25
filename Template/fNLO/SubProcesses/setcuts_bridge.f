      subroutine setcuts
c-----------------------------------------------------------------------
c     Keep Python-generated model declarations and assignments out of
c     the free-form setcuts module, and mirror its particle classes in
c     the historical COMMON block for fixed-form callers.
c-----------------------------------------------------------------------
      use setcuts_module, only: setcuts_impl
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'coupl.inc'
      double precision zero
      parameter (zero=0d0)
      double precision pmass(nexternal)
      common /to_mass/pmass
      integer idup(nexternal,maxproc)
      integer mothup(2,nexternal,maxproc)
      integer icolup(2,nexternal,maxflow),niprocs
      common /c_leshouche_inc/idup,mothup,icolup,niprocs
      logical is_a_j(nexternal),is_a_lp(nexternal)
      logical is_a_lm(nexternal),is_a_ph(nexternal)
      common /to_specisa/is_a_j,is_a_lp,is_a_lm,is_a_ph
      double precision etmin(nincoming+1:nexternal-1)
      double precision etmax(nincoming+1:nexternal-1)
      double precision mxxmin(nincoming+1:nexternal-1,
     $     nincoming+1:nexternal-1)
      common /to_cuts/etmin,etmax,mxxmin
      include 'pmass.inc'

      call setcuts_impl(nf,pmass,idup,etmin,etmax,mxxmin,
     $     is_a_j,is_a_lp,is_a_lm,is_a_ph)
      return
      end


      subroutine set_tau_min
c-----------------------------------------------------------------------
c     born_props.inc is generated for each subprocess.  Materialize its
c     assignments here, then pass ordinary arrays to the F90 module.
c-----------------------------------------------------------------------
      use setcuts_module, only: set_tau_min_impl
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'coupl.inc'
      double precision zero
      parameter (zero=0d0)
      integer i,j
      integer pow(-nexternal:0,lmaxconfigs)
      double precision pmass(-nexternal:0,lmaxconfigs)
      double precision pwidth(-nexternal:0,lmaxconfigs)
      integer itree(2,-max_branch:-1),iconf
      common /to_itree/itree,iconf
      integer nfksprocess
      common /c_nfksprocess/nfksprocess
      double precision tau_born_lower_bound
      double precision tau_lower_bound_resonance,tau_lower_bound
      common /ctau_lower_bound/tau_born_lower_bound,
     $     tau_lower_bound_resonance,tau_lower_bound
      double precision cbw_mass(-1:1,-nexternal:-1)
      double precision cbw_width(-1:1,-nexternal:-1)
      integer cbw_level_max,cbw(-nexternal:-1)
      integer cbw_level(-nexternal:-1)
      common /c_conflictingbw/cbw_mass,cbw_width,cbw_level_max,
     $     cbw,cbw_level
      double precision s_mass(-nexternal:nexternal)
      common /to_phase_space_s_channel/s_mass
      integer idup(nexternal,maxproc)
      integer mothup(2,nexternal,maxproc)
      integer icolup(2,nexternal,maxflow),niprocs
      common /c_leshouche_inc/idup,mothup,icolup,niprocs
      double precision emass(nexternal)
      common /to_mass/emass
      logical is_a_j(nexternal),is_a_lp(nexternal)
      logical is_a_lm(nexternal),is_a_ph(nexternal)
      common /to_specisa/is_a_j,is_a_lp,is_a_lm,is_a_ph
      double precision etmin(nincoming+1:nexternal-1)
      double precision etmax(nincoming+1:nexternal-1)
      double precision mxxmin(nincoming+1:nexternal-1,
     $     nincoming+1:nexternal-1)
      common /to_cuts/etmin,etmax,mxxmin

      do i=1,lmaxconfigs
         do j=-nexternal,0
            pmass(j,i)=0d0
            pwidth(j,i)=0d0
            pow(j,i)=0
         enddo
      enddo
      include 'born_props.inc'

      call set_tau_min_impl(pmass,pwidth,itree,iconf,nfksprocess,
     $     idup,emass,etmin,etmax,mxxmin,is_a_j,is_a_lp,is_a_lm,
     $     is_a_ph,tau_born_lower_bound,tau_lower_bound_resonance,
     $     tau_lower_bound,cbw_mass,cbw_width,cbw_level_max,cbw,
     $     cbw_level,s_mass)
      return
      end


      subroutine schan_order(ns_channel,order)
c-----------------------------------------------------------------------
c     Preserve the historical external ABI while the ordering cache and
c     algorithm live in setcuts_module.
c-----------------------------------------------------------------------
      use setcuts_module, only: schan_order_impl
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      integer ns_channel,order(-nexternal:0)
      integer itree(2,-max_branch:-1),iconf
      common /to_itree/itree,iconf

      call schan_order_impl(ns_channel,order,itree)
      return
      end
