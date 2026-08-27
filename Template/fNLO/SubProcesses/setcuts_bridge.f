      subroutine setcuts
c-----------------------------------------------------------------------
c     Keep Python-generated model declarations and assignments out of
c     the free-form setcuts module, and mirror its particle classes in
c     the module-owned shared state for fixed-form callers.
c-----------------------------------------------------------------------
      use setcuts_module, only: setcuts_impl
      use fnlo_process_common, only: pmass=>particle_masses,idup,
     $     etmin,etmax,mxxmin,is_a_j,is_a_lp,is_a_lm,is_a_ph
     $     ,from_decay
      implicit none
      include 'coupl.inc'
      double precision zero
      parameter (zero=0d0)
      include 'pmass.inc'

      from_decay=.false.

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
      use fnlo_process_common, only: nexternal,lmaxconfigs,
     $     itree=>config_tree,iconf=>config_index,nfksprocess,
     $     tau_born_lower_bound,tau_lower_bound_resonance,
     $     tau_lower_bound,cbw_mass,cbw_width,cbw_level_max,cbw,
     $     cbw_level,s_mass=>schannel_masses,idup,
     $     emass=>particle_masses,is_a_j,is_a_lp,is_a_lm,
     $     is_a_ph,etmin,etmax,mxxmin
      implicit none
      include 'coupl.inc'
      double precision zero
      parameter (zero=0d0)
      integer i,j
      double precision pmass(-nexternal:0,lmaxconfigs)
      double precision pwidth(-nexternal:0,lmaxconfigs)

      do i=1,lmaxconfigs
         do j=-nexternal,0
            pmass(j,i)=0d0
            pwidth(j,i)=0d0
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
      use fnlo_process_common, only: nexternal,itree=>config_tree
      implicit none
      integer ns_channel,order(-nexternal:0)

      call schan_order_impl(ns_channel,order,itree)
      return
      end
