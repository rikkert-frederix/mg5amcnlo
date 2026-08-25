c Python-generated COMMON boundary for the allocatable cuts state.
c Call after SETCUTS, including whenever the selected FKS channel changes.

      subroutine sync_cuts_bridge_state()
      use cuts_module, only: initialize_cuts_runtime_state,
     $     initialize_cuts_event_state
      implicit none
      include 'nexternal.inc'
      include 'genps.inc'
      include 'orders.inc'
      logical split_type_used(nsplitorders)
      common /to_split_type_used/split_type_used
      integer need_matching_s(nexternal),need_matching_h(nexternal)
      integer need_matching_cuts(nexternal)
      common /c_need_matching/need_matching_s,need_matching_h,
     $     need_matching_cuts
      logical particle_tag(nexternal)
      common /c_particle_tag/particle_tag
      double precision etmin(nincoming+1:nexternal-1)
      double precision etmax(nincoming+1:nexternal-1)
      double precision mxxmin(nincoming+1:nexternal-1,
     $     nincoming+1:nexternal-1)
      common /to_cuts/etmin,etmax,mxxmin
      logical is_a_j(nexternal),is_a_lp(nexternal)
      logical is_a_lm(nexternal),is_a_ph(nexternal)
      common /to_specisa/is_a_j,is_a_lp,is_a_lm,is_a_ph
      double precision fxfx_ren_scales(0:nexternal)
      double precision fxfx_fac_scale(2)
      integer nfxfx_ren_scales
      common /c_fxfx_scales/fxfx_ren_scales,nfxfx_ren_scales,
     $     fxfx_fac_scale
      double precision pmass(nexternal)
      common /to_mass/pmass
      integer idup(nexternal,maxproc)
      integer mothup(2,nexternal,maxproc)
      integer icolup(2,nexternal,maxflow),niprocs
      common /c_leshouche_inc/idup,mothup,icolup,niprocs
      double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
      common /parton_cms_stuff/ybst_til_tolab,ybst_til_tocm,
     $     sqrtshat,shat

      call initialize_cuts_runtime_state(split_type_used,
     $     need_matching_s,need_matching_h,need_matching_cuts,
     $     particle_tag,etmin,etmax,mxxmin,is_a_j,is_a_lp,is_a_lm,
     $     is_a_ph,fxfx_ren_scales,nfxfx_ren_scales,
     $     fxfx_fac_scale)
      call initialize_cuts_event_state(pmass,idup,ybst_til_tolab)
      return
      end


c The generated FastJet wrapper name exceeds the Fortran 95 identifier
c limit, so keep this short standards-compliant adapter at the ABI edge.
      subroutine cuts_fastjet_etamax(pqcd,nn,rfj,sycut,etamax,
     $     palg,pjet,njet,jet)
      implicit none
      integer nn,njet,jet(nn)
      double precision pqcd(0:3,nn),rfj,sycut,etamax,palg
      double precision pjet(0:3,nn)
      call amcatnlo_fastjetppgenkt_etamax_timed(pqcd,nn,rfj,
     $     sycut,etamax,palg,pjet,njet,jet)
      return
      end
