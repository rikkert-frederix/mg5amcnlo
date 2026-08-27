module setcuts_module
  use process_dimensions, only: nexternal, nincoming, max_branch, &
       lmaxconfigs, maxproc, fks_configs, validate_process_dimensions
  use fks_metadata, only: validate_fks_metadata, fks_i_d, fks_j_d
  use decay_chain_metadata, only: initialize_decay_chain_metadata, &
       has_decay_chains, born_context, context_core_count, &
       core_target_kind, core_target_id, direct_leg_target, &
       decay_node_target
  use decay_chain_kinematics, only: initialize_decay_chain_kinematics, &
       core_mass
  use run_state, only: lpp, ebeam, maxjetflavor, gamma_is_j, pdg_cut, &
       ptmin4pdg, ptmax4pdg, mxxmin4pdg, mxxpart_antipart, ptj, &
       ptgmin, ptl, mll, mll_sf
  use mint_module, only: maxchannels, ichan, iconfig, new_point
  use fks_random_module, only: random_unit_interval
  implicit none
  private

  double precision, parameter :: vtiny = 1d-8
  logical, save :: setcuts_state_initialized = .false.
  double precision, save :: stot = 0d0
  double precision, allocatable, save :: taumin(:, :)
  double precision, allocatable, save :: taumin_s(:, :)
  double precision, allocatable, save :: taumin_j(:, :)
  double precision, allocatable, save :: mass_min(:, :)
  integer, allocatable, save :: cbw_fks_level_max(:, :)
  integer, allocatable, save :: cbw_fks(:, :, :)
  integer, allocatable, save :: cbw_fks_level(:, :, :)
  double precision, allocatable, save :: cbw_fks_mass(:, :, :, :)
  double precision, allocatable, save :: cbw_fks_width(:, :, :, :)
  double precision, allocatable, save :: s_mass_fks(:, :, :)
  logical, allocatable, save :: firsttime_chans(:)
  integer, allocatable, save :: saved_schan_order(:)

  public :: setcuts_impl, set_tau_min_impl, set_decay_tau_min_impl
  public :: schan_order_impl

contains

  subroutine setcuts_impl(nf, pmass, idup, etmin, etmax, mxxmin, &
       is_a_j_compat, is_a_lp_compat, is_a_lm_compat, is_a_ph_compat)
    implicit none
    double precision, intent(in) :: nf
    double precision, intent(in) :: pmass(:)
    integer, intent(in) :: idup(:, :)
    double precision, intent(inout) :: etmin(nincoming + 1:)
    double precision, intent(inout) :: etmax(nincoming + 1:)
    double precision, intent(inout) :: mxxmin(nincoming + 1:, &
         nincoming + 1:)
    logical, intent(inout) :: is_a_j_compat(:), is_a_lp_compat(:)
    logical, intent(inout) :: is_a_lm_compat(:), is_a_ph_compat(:)
    integer :: i, j, k

    call validate_process_dimensions()
    call validate_setcuts_inputs(pmass, idup, etmin, etmax, mxxmin, &
         is_a_j_compat, is_a_lp_compat, is_a_lm_compat, is_a_ph_compat)

    if (nincoming == 1) then
      lpp(1) = 0
      lpp(2) = 0
      ebeam(1) = pmass(1) / 2d0
      ebeam(2) = pmass(1) / 2d0
    end if

    if (maxjetflavor < int(nf)) then
      write(*, '(a,i3,a,i3)') &
           'WARNING: the value of maxjetflavorspecified in the run_card (', &
           maxjetflavor, ') is inconsistent with the number of light ' // &
           'flavours inthe model. Hence it will be set to:', int(nf)
      maxjetflavor = int(nf)
    end if

    do i = nincoming + 1, nexternal
      is_a_j_compat(i) = .false.
      is_a_lp_compat(i) = .false.
      is_a_lm_compat(i) = .false.
      is_a_ph_compat(i) = .false.
      if (abs(idup(i, 1)) <= maxjetflavor) is_a_j_compat(i) = .true.
      if (abs(idup(i, 1)) == 21) is_a_j_compat(i) = .true.

      if (idup(i, 1) == 11) is_a_lm_compat(i) = .true.
      if (idup(i, 1) == 13) is_a_lm_compat(i) = .true.
      if (idup(i, 1) == 15) is_a_lm_compat(i) = .true.
      if (idup(i, 1) == -11) is_a_lp_compat(i) = .true.
      if (idup(i, 1) == -13) is_a_lp_compat(i) = .true.
      if (idup(i, 1) == -15) is_a_lp_compat(i) = .true.

      if (idup(i, 1) == 22 .and. .not. gamma_is_j) then
        is_a_ph_compat(i) = .true.
      end if
      if (idup(i, 1) == 22 .and. gamma_is_j) is_a_j_compat(i) = .true.
    end do

    do i = nincoming + 1, nexternal - 1
      etmin(i) = 0d0
      etmax(i) = -1d0
      do j = i, nexternal - 1
        mxxmin(i, j) = 0d0
      end do
    end do

    if (pdg_cut(1) /= 0) then
      do j = 1, pdg_cut(0)
        do i = nincoming + 1, nexternal - 1
          if (abs(idup(i, 1)) /= pdg_cut(j)) cycle
          if (pmass(i) == 0d0) then
            write(*, *) 'Illegal use of pdg specific cut.'
            write(*, *) 'For NLO process, only massive particle can be included'
            stop 1
          end if
          if (is_a_lp_compat(i) .or. is_a_lm_compat(i) .or. &
              is_a_j_compat(i) .or. is_a_ph_compat(i)) then
            write(*, *) 'Illegal use of pdg specific cut.'
            write(*, *) 'This can not be used for jet/lepton/photon/gluon'
            stop 1
          end if
          etmin(i) = ptmin4pdg(j)
          etmax(i) = ptmax4pdg(j)
          if (mxxmin4pdg(j) /= 0d0) then
            do k = i + 1, nexternal - 1
              if (mxxpart_antipart(j)) then
                if (idup(k, 1) == -idup(i, 1)) then
                  mxxmin(i, k) = mxxmin4pdg(j)
                end if
              else
                if (abs(idup(k, 1)) == pdg_cut(j)) then
                  mxxmin(i, k) = mxxmin4pdg(j)
                end if
              end if
            end do
          end if
        end do
      end do
    end if

  end subroutine setcuts_impl


subroutine set_decay_tau_min_impl(nfksprocess, etmin, is_a_j_compat, &
     is_a_lp_compat, is_a_lm_compat, is_a_ph_compat, &
     tau_born_lower_bound, tau_lower_bound_resonance, tau_lower_bound, &
     cbw_mass, cbw_width, cbw_level_max, cbw, cbw_level, s_mass)
! Set the tau bounds for an exact-NWA decay-chain process from the
! undecayed production core.  A forced decay contributes its on-shell
! parent mass exactly once; any enabled cuts on its visible daughters
! are enforced only after the decay momenta have been generated.
implicit none
integer, intent(in) :: nfksprocess
double precision, intent(in) :: etmin(nincoming + 1:)
logical, intent(in) :: is_a_j_compat(:), is_a_lp_compat(:)
logical, intent(in) :: is_a_lm_compat(:), is_a_ph_compat(:)
double precision, intent(out) :: tau_born_lower_bound
double precision, intent(out) :: tau_lower_bound_resonance
double precision, intent(out) :: tau_lower_bound
double precision, intent(out) :: cbw_mass(-1:, -nexternal:)
double precision, intent(out) :: cbw_width(-1:, -nexternal:)
integer, intent(out) :: cbw_level_max
integer, intent(out) :: cbw(-nexternal:), cbw_level(-nexternal:)
double precision, intent(out) :: s_mass(-nexternal:)
integer :: context, leg, target, j_fks
double precision :: production_mass, cut_mass, leg_mass, leg_bound
double precision :: collider_energy_squared

call validate_process_dimensions()
call validate_fks_metadata()
call initialize_decay_chain_metadata()
if (.not. has_decay_chains()) then
   call fail_setcuts('decay tau bounds requested without decay metadata')
end if
call initialize_decay_chain_kinematics()
if (nfksprocess < 1 .or. nfksprocess > fks_configs) then
   call fail_setcuts('decay tau bounds received an invalid FKS configuration')
end if
if (size(etmin) /= nexternal - nincoming - 1 .or. &
    size(is_a_j_compat) /= nexternal .or. &
    size(is_a_lp_compat) /= nexternal .or. &
    size(is_a_lm_compat) /= nexternal .or. &
    size(is_a_ph_compat) /= nexternal) then
   call fail_setcuts('decay tau bounds received inconsistent cut data')
end if
if (size(cbw_mass, 1) /= 3 .or. size(cbw_mass, 2) /= nexternal .or. &
    any(shape(cbw_width) /= shape(cbw_mass)) .or. &
    size(cbw) /= nexternal .or. size(cbw_level) /= nexternal .or. &
    size(s_mass) /= 2*nexternal + 1) then
   call fail_setcuts('decay tau bounds received inconsistent output storage')
end if

context = born_context()
production_mass = 0d0
cut_mass = 0d0
do leg = nincoming + 1, context_core_count(context)
   leg_mass = core_mass(context, leg)
   production_mass = production_mass + leg_mass
   leg_bound = leg_mass

   select case (core_target_kind(context, leg))
   case (decay_node_target)
! A decay node is an external on-shell particle of the production core.
! In particular, do not add ptj/ptl/ptgmin once for every decay leaf.
      continue
   case (direct_leg_target)
      target = core_target_id(context, leg)
      if (target <= nincoming .or. target >= nexternal) then
         call fail_setcuts('production leg has an invalid visible target')
      end if
      if (is_a_j_compat(target)) then
         leg_bound = dsqrt(ptj**2 + leg_mass**2)
      elseif (is_a_ph_compat(target)) then
         leg_bound = dsqrt(ptgmin**2 + leg_mass**2)
      elseif (is_a_lp_compat(target) .or. is_a_lm_compat(target)) then
         leg_bound = dsqrt(ptl**2 + leg_mass**2)
      elseif (etmin(target) > 0d0) then
         leg_bound = dsqrt(etmin(target)**2 + leg_mass**2)
      end if
   case default
      call fail_setcuts('production leg has an unknown decay-map target')
   end select
   cut_mass = cut_mass + leg_bound
end do

collider_energy_squared = 4d0*ebeam(1)*ebeam(2)
if (collider_energy_squared <= 0d0) then
   call fail_setcuts('decay tau bounds received invalid beam energies')
end if
j_fks = fks_j_d(nfksprocess)
if (j_fks > nincoming) then
   tau_born_lower_bound = cut_mass**2/collider_energy_squared
else
! For initial-state radiation the real parton can make an otherwise soft
! underlying Born pass the production cuts, so its hard Born bound is mass-only.
   tau_born_lower_bound = production_mass**2/collider_energy_squared
end if
tau_lower_bound = cut_mass**2/collider_energy_squared
tau_lower_bound_resonance = tau_lower_bound

! The decay-chain phase space generates the production core directly and
! therefore does not use the visible-process propagator tree.
cbw_mass = 0d0
cbw_width = 0d0
cbw_level_max = 0
cbw = 0
cbw_level = 0
s_mass = 0d0
end subroutine set_decay_tau_min_impl


subroutine set_tau_min_impl(pmass, pwidth, itree, iconf, nfksprocess, &
     idup, emass, etmin, etmax, mxxmin, is_a_j_compat, is_a_lp_compat, &
     is_a_lm_compat, is_a_ph_compat, tau_born_lower_bound, &
     tau_lower_bound_resonance, tau_lower_bound, cbw_mass, cbw_width, &
     cbw_level_max, cbw, cbw_level, s_mass)
! Sets the lower bound for tau=x1*x2, using information on particle
! masses and on the jet minimum pt, as entered in run_card.dat,
! variable ptj
implicit none
double precision, intent(in) :: pmass(-nexternal:, :)
double precision, intent(in) :: pwidth(-nexternal:, :)
integer, intent(in) :: itree(2, -max_branch:-1), iconf, nfksprocess
integer, intent(in) :: idup(:, :)
double precision, intent(in) :: emass(:)
double precision, intent(in) :: etmin(nincoming + 1:)
double precision, intent(in) :: etmax(nincoming + 1:)
double precision, intent(in) :: mxxmin(nincoming + 1:, nincoming + 1:)
logical, intent(in) :: is_a_j_compat(:), is_a_lp_compat(:)
logical, intent(in) :: is_a_lm_compat(:), is_a_ph_compat(:)
double precision, intent(out) :: tau_born_lower_bound
double precision, intent(out) :: tau_lower_bound_resonance
double precision, intent(out) :: tau_lower_bound
double precision, intent(inout) :: cbw_mass(-1:, -nexternal:)
double precision, intent(inout) :: cbw_width(-1:, -nexternal:)
integer, intent(out) :: cbw_level_max
integer, intent(out) :: cbw(-nexternal:), cbw_level(-nexternal:)
double precision, intent(out) :: s_mass(-nexternal:)
double precision :: xk(-nexternal:nexternal)
integer :: i,j,k,d1,d2,iFKS,nt
double precision xm(-nexternal:nexternal),xm1,xm2,xmi
double precision xw(-nexternal:nexternal),xw1,xw2,xwi
integer tsign,i_fks,j_fks
! BW stuff
double precision :: masslow(-nexternal:-1),widthlow(-nexternal:-1)
double precision :: sum_all_s
integer t_channel
double precision smin_update , mxx
integer nb_iden_pdg

call validate_process_dimensions()
call validate_fks_metadata()
call initialize_setcuts_state()
call validate_tau_inputs(pmass, pwidth, itree, iconf, nfksprocess, &
     idup, emass, etmin, etmax, mxxmin, is_a_j_compat, is_a_lp_compat, &
     is_a_lm_compat, is_a_ph_compat, cbw_mass, cbw_width, cbw, &
     cbw_level, s_mass)

! The following assumes that light QCD particles are at the end of the
! list. Exclude one of them (i_fks) to set tau bound at the Born level
! This sets a hard cut in the minimal shat of the Born phase-space
! generation.
!
! The contribution from ptj should be treated only as a 'soft lower
! bound' if j_fks is initial state: the real-emission i_fks parton is
! not necessarily the softest.  Therefore, it could be that even though
! the Born does not have enough energy to pass the cuts set by ptj, the
! event could.
if (firsttime_chans(ichan)) then
   do i=-nexternal,nexternal
      xm(i)=0d0
      xw(i)=0d0
      mass_min(i,ichan)=0d0
   end do
   firsttime_chans(ichan)=.false.
   do iFKS=1,fks_configs
      j_fks=fks_j_d(iFKS)
      i_fks=fks_i_d(iFKS)
      taumin(iFKS,ichan)=0.d0
      taumin_s(iFKS,ichan)=0.d0
      taumin_j(iFKS,ichan)=0.d0
      do i=nincoming+1,nexternal
! Skip i_fks
         if (i.eq.i_fks) cycle
! Add the minimal jet pTs to tau
         if (is_a_j_compat(i)) then
            if  (j_fks.gt.nincoming .and. j_fks.lt.nexternal) then
               taumin(iFKS,ichan)=taumin(iFKS,ichan)+dsqrt(ptj**2 +emass(i)**2)
               taumin_s(iFKS,ichan)=taumin_s(iFKS,ichan)+dsqrt(ptj**2 +emass(i)**2)
               taumin_j(iFKS,ichan)=taumin_j(iFKS,ichan)+dsqrt(ptj**2 +emass(i)**2)
            elseif (j_fks.ge.1 .and. j_fks.le.nincoming) then
               taumin(iFKS,ichan)=taumin(iFKS,ichan)+emass(i)
               taumin_s(iFKS,ichan)=taumin_s(iFKS,ichan)+dsqrt(ptj**2 +emass(i)**2)
               taumin_j(iFKS,ichan)=taumin_j(iFKS,ichan)+dsqrt(ptj**2 +emass(i)**2)
            elseif (j_fks.eq.nexternal) then
               write (*,*) &
  &                     'ERROR, j_fks cannot be the final parton' &
  &                     ,j_fks
               stop
            else
               write (*,*) 'ERROR, j_fks not correctly defined' &
  &                     ,j_fks
               stop
            endif
            xm(i)=emass(i)+ptj
! Add the minimal photon pTs to tau
         elseif (is_a_ph_compat(i)) then
            if (abs(emass(i)).gt.vtiny) then
               write (*,*) 'Error in set_tau_min in setcuts.f:'
               write (*,*) 'mass of a photon should be zero',i &
  &                     ,emass(i)
               stop
            endif
            if  (j_fks.gt.nincoming) &
  &                  taumin(iFKS,ichan)=taumin(iFKS,ichan)+ptgmin
            taumin_s(iFKS,ichan)=taumin_s(iFKS,ichan)+ptgmin
            taumin_j(iFKS,ichan)=taumin_j(iFKS,ichan)+ptgmin
            xm(i)=emass(i)+ptgmin
         elseif (is_a_lp_compat(i)) then
! Add the postively charged lepton pTs to tau
            if (j_fks.gt.nincoming) then
               taumin(iFKS,ichan)=taumin(iFKS,ichan)+dsqrt(ptl**2+emass(i)**2)
            else
               taumin(iFKS,ichan)=taumin(iFKS,ichan)+emass(i)
            endif
            taumin_s(iFKS,ichan)=taumin_s(iFKS,ichan)+dsqrt(emass(i)**2+ptl**2)
            taumin_j(iFKS,ichan)=taumin_j(iFKS,ichan)+dsqrt(emass(i)**2+ptl**2)
            xm(i)=emass(i)+ptl
! Add the lepton invariant mass to tau if there is at least another
! lepton of opposite charge. (Only add half of it, i.e. 'the part
! contributing from this lepton'). Remove possible overcounting with the
! lepton pT
            do j=nincoming+1,nexternal
               if (is_a_lm_compat(j) .and. idup(i,1).eq.-idup(j,1) .and. &
  &                     (mll_sf.ne.0d0 .or. mll.ne.0d0) ) then
                  if (j_fks.gt.nincoming) &
  &                        taumin(iFKS,ichan) = taumin(iFKS,ichan)-dsqrt(ptl**2+emass(i)**2) + &
  &                               max(mll/2d0,mll_sf/2d0,dsqrt(ptl**2+emass(i)**2))
                  taumin_s(iFKS,ichan) = taumin_s(iFKS,ichan)-dsqrt(ptl**2+emass(i)**2) &
  &                        + max(mll/2d0,mll_sf/2d0,dsqrt(ptl**2+emass(i)**2))
                  taumin_j(iFKS,ichan) = taumin_j(iFKS,ichan)-dsqrt(ptl**2+emass(i)**2) &
  &                        + max(mll/2d0,mll_sf/2d0,dsqrt(ptl**2+emass(i)**2))
                  xm(i)=xm(i)-ptl-emass(i)+max(mll/2d0,mll_sf/2d0 &
  &                        ,ptl+emass(i))
                  exit
               elseif (is_a_lm_compat(j) .and. mll.ne.0d0) then
                  if (j_fks.gt.nincoming) &
  &                        taumin(iFKS,ichan)= taumin(iFKS,ichan)-dsqrt(ptl**2+emass(i)**2) + &
  &                                      max(mll/2d0,dsqrt(ptl**2+emass(i)**2))
                  taumin_s(iFKS,ichan) = taumin_s(iFKS,ichan)-dsqrt(ptl**2+emass(i)**2) &
  &                        + max(mll/2d0, dsqrt(ptl**2+emass(i)**2))
                  taumin_j(iFKS,ichan) = taumin_j(iFKS,ichan)-dsqrt(ptl**2+emass(i)**2) &
  &                        + max(mll/2d0,dsqrt(ptl**2+emass(i)**2))
                  xm(i)=xm(i)-ptl-emass(i)+max(mll/2d0,ptl &
  &                        +emass(i))
                  exit
               endif
            enddo
         elseif (is_a_lm_compat(i)) then
! Add the negatively charged lepton pTs to tau
            if (j_fks.gt.nincoming) then
               taumin(iFKS,ichan)=taumin(iFKS,ichan)+dsqrt(ptl**2+emass(i)**2)
            else
               taumin(iFKS,ichan)=taumin(iFKS,ichan)+emass(i)
            endif
            taumin_s(iFKS,ichan)=taumin_s(iFKS,ichan)+dsqrt(ptl**2+emass(i)**2)
            taumin_j(iFKS,ichan)=taumin_j(iFKS,ichan)+dsqrt(ptl**2+emass(i)**2)
            xm(i)=emass(i)+ptl
! Add the lepton invariant mass to tau if there is at least another
! lepton of opposite charge. (Only add half of it, i.e. 'the part
! contributing from this lepton'). Remove possible overcounting with the
! lepton pT
            do j=nincoming+1,nexternal
               if (is_a_lp_compat(j) .and. idup(i,1).eq.-idup(j,1) .and. &
  &                     (mll_sf.ne.0d0 .or. mll.ne.0d0) ) then
                  if (j_fks.gt.nincoming) &
  &                        taumin(iFKS,ichan) = taumin(iFKS,ichan)-dsqrt(ptl**2+emass(i)**2) + &
  &                               max(mll/2d0,mll_sf/2d0,dsqrt(ptl**2+emass(i)**2))
                  taumin_s(iFKS,ichan) = taumin_s(iFKS,ichan)-dsqrt(ptl**2+emass(i)**2) &
  &                        + max(mll/2d0,mll_sf/2d0,dsqrt(ptl**2+emass(i)**2))
                  taumin_j(iFKS,ichan) = taumin_j(iFKS,ichan)-dsqrt(ptl**2+emass(i)**2) &
  &                        + max(mll/2d0,mll_sf/2d0,dsqrt(ptl**2+emass(i)**2))
                  xm(i)=xm(i)-ptl-emass(i)+max(mll/2d0,mll_sf/2d0 &
  &                        ,ptl+emass(i))
                  exit
               elseif (is_a_lp_compat(j) .and. mll.ne.0d0) then
                  if (j_fks.gt.nincoming) &
  &                        taumin(iFKS,ichan) = taumin(iFKS,ichan)-dsqrt(ptl**2+emass(i)**2) + &
  &                                       max(mll/2d0,dsqrt(ptl**2+emass(i)**2))
                  taumin_s(iFKS,ichan) = taumin_s(iFKS,ichan)-dsqrt(ptl**2+emass(i)**2) &
  &                        + max(mll/2d0,dsqrt(ptl**2+emass(i)**2))
                  taumin_j(iFKS,ichan) = taumin_j(iFKS,ichan)-dsqrt(ptl**2+emass(i)**2) &
  &                        + max(mll/2d0,dsqrt(ptl**2+emass(i)**2))
                  xm(i)=xm(i)-ptl-emass(i)+max(mll/2d0,ptl &
  &                        +emass(i))
                  exit
               endif
            enddo
         else
            if (i.eq.nexternal)then
                  taumin(iFKS,ichan)=taumin(iFKS,ichan) + emass(i)
                  taumin_s(iFKS,ichan)=taumin_s(iFKS,ichan) +  emass(i)
                  taumin_j(iFKS,ichan)=taumin_j(iFKS,ichan) + emass(i)
                  xm(i) = emass(i)
            else
               smin_update = 0
               nb_iden_pdg = 1
               mxx = 0d0
!                    assume smin apply always on the same set of particle
               do j=nincoming+1,nexternal-1
                  if (mxxmin(i,j).ne.0d0.or.mxxmin(j,i).ne.0d0) then
                     nb_iden_pdg = nb_iden_pdg +1
                     if (mxx.eq.0d0) mxx = max(mxxmin(i,j), mxxmin(j,i))
                  endif
               enddo
               ! S >= (2*N-N^2)*M1^2 + (N^2-N)/2 * Mxx^2
               smin_update = nb_iden_pdg*((2-nb_iden_pdg)*emass(i)**2 + (nb_iden_pdg-1)/2.*mxx**2)
               ! compare with the update from pt cut
               if (smin_update.lt.nb_iden_pdg**2*(etmin(i)**2 + emass(i)**2))then
                  ! the pt is more restrictive
                  smin_update = dsqrt(etmin(i)**2 + emass(i)**2)
               else
                  smin_update = dsqrt(smin_update)/nb_iden_pdg ! share over N particle, and change dimension
               endif
               ! update in sqrt(s) so take the
               if  (j_fks.gt.nincoming) then
                  taumin(iFKS,ichan)=taumin(iFKS,ichan) + smin_update
               else
                  taumin(iFKS,ichan)=taumin(iFKS,ichan) + emass(i)
               endif
               taumin_s(iFKS,ichan)=taumin_s(iFKS,ichan) + smin_update
               taumin_j(iFKS,ichan)=taumin_j(iFKS,ichan) + smin_update
               xm(i) = smin_update
            endif
         endif
         xw(i)=0d0
      enddo
      stot = 4d0*ebeam(1)*ebeam(2)
      tau_Born_lower_bound=taumin(iFKS,ichan)**2/stot
      tau_lower_bound=taumin_j(iFKS,ichan)**2/stot
!
! Also find the minimum lower bound if all internal s-channel particles
! were on-shell
      tsign=-1
      nt=0
      do i=-1,-(nexternal-3),-1 ! All propagators
         if ( itree(1,i) .eq. 1 .or. itree(1,i) .eq. 2 ) tsign=1
         if (tsign.eq.-1) then ! s-channels
            d1=itree(1,i)
            d2=itree(2,i)
! If daughter is a jet, we should treat the ptj as a mass. Except if
! d1=nexternal, because we check the Born, so final parton should be
! skipped. [This is already done above; also for the leptons]
            xm1=xm(d1)
            xm2=xm(d2)
            xw1=xw(d1)
            xw2=xw(d2)
! On-shell mass of the intermediate resonance
            xmi=pmass(i,iconf)
! Width of the intermediate resonance
            xwi=pwidth(i,iconf)
! Set the intermediate mass equal to the max of its actual mass and
! the sum of the masses of the two daugters.
            if (xmi.gt.xm1+xm2) then
               xm(i)=xmi
               xw(i)=xwi
            else
               xm(i)=xm1+xm2
               xw(i)=xw1+xw2 ! just sum the widths
            endif
! Add the new mass to the bound. To avoid double counting, we should
! subtract the daughters, because they are already included above or in
! the previous iteration of the loop
            taumin_s(iFKS,ichan)=taumin_s(iFKS,ichan)+xm(i)-xm1-xm2
         else             ! t-channels
            if (i.eq.-(nexternal-3)) then
               xm(i)=0d0
               cycle
            endif
            nt=nt+1
            d1=itree(2,i) ! only use 2nd daughter (which is the outgoing one)
            xm1=xm(d1)
            if (nt.gt.1) xm1=max(xm1,xk(nt-1))
            xk(nt)=xm1
            j=i-1         ! this is the closest to p2
            d2=itree(2,j)
            xm2=xm(d2)
            xm(i)=min(xm1,xm2)
         endif
      enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Determine the "minimal" s-channel invariant masses
      do i=nincoming+1,nexternal-1
         s_mass_FKS(iFKS,i,ichan)=xm(i)**2
      enddo
      do i=-1,-(nexternal-3),-1 ! All propagators
         s_mass_FKS(iFKS,i,ichan)=xm(i)**2
      enddo
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Determine the conflicting Breit-Wigner's. Note that xm(i) contains the
! mass of the BW
      do i=nincoming+1,nexternal-1
         mass_min(i,ichan)=xm(i) ! minimal allowed resonance mass (including masses set by cuts)
      enddo
      cBW_FKS_level_max(iFKS,ichan)=0
      t_channel=0
      do i=-1,-(nexternal-3),-1 ! All propagators
         cBW_FKS_mass(iFKS,1,i,ichan)=0d0
         cBW_FKS_width(iFKS,1,i,ichan)=0d0
         cBW_FKS_mass(iFKS,-1,i,ichan)=0d0
         cBW_FKS_width(iFKS,-1,i,ichan)=0d0
         masslow(i)=9d99
         widthlow(i)=0d0
         if ( itree(1,i).eq.1 .or. itree(1,i).eq.2 ) t_channel=i
         if (t_channel.ne.0) exit ! only s-channels
         mass_min(i,ichan)=mass_min(itree(1,i),ichan) &
  &               +mass_min(itree(2,i),ichan)
         if (xm(i).lt.mass_min(i,ichan)-vtiny) then
            write (*,*) &
  &                  'ERROR in the determination of conflicting BW',i &
  &                  ,xm(i),mass_min(i,ichan)
            stop
         endif
         if (pmass(i,iconf).lt.xm(i) .and. &
  &               pwidth(i,iconf).gt.0d0) then
!     Possible conflict in BW
            if (pmass(i,iconf).lt.mass_min(i,ichan)) then
!     Resonance can never go on-shell due to the kinematics of the event
               cBW_FKS(iFKS,i,ichan)=2
               cBW_FKS_level(iFKS,i,ichan)=0
            elseif(pmass(i,iconf).lt.xm(i)) then
!     Conflicting Breit-Wigner
               cBW_FKS(iFKS,i,ichan)=1
               cBW_FKS_level(iFKS,i,ichan)=1
               cBW_FKS_level_max(iFKS,ichan)=max(cBW_FKS_level_max(iFKS,ichan) &
  &                     ,cBW_FKS_level(iFKS,i,ichan))
!     Set here the mass (and width) of the alternative mass; it's the
!     sum of daughter masses. (2nd argument is '1', because this
!     alternative mass is LARGER than the resonance mass).
               cBW_FKS_mass(iFKS,1,i,ichan)=xm(i)
               cBW_FKS_width(iFKS,1,i,ichan)=xw(i)
            endif
!     set the daughters also as conflicting (recursively)
            masslow(i)=pmass(i,iconf)
            widthlow(i)=pwidth(i,iconf)
            do j=i,-1
               if (cBW_FKS(iFKS,j,ichan).eq.0) cycle
               do k=1,2   ! loop over the 2 daughters
                  if (itree(k,j).ge.0) cycle
                  if (cBW_FKS(iFKS,itree(k,j),ichan).eq.2) cycle
                  cBW_FKS(iFKS,itree(k,j),ichan)=1
                  cBW_FKS_level(iFKS,itree(k,j),ichan)= &
  &                        cBW_FKS_level(iFKS,j,ichan)+1
                  cBW_FKS_level_max(iFKS,ichan)= &
  &                        max(cBW_FKS_level_max(iFKS,ichan) &
  &                        ,cBW_FKS_level(iFKS,itree(k,j),ichan))
!     Set here the mass (and width) of the alternative mass; it's the
!     difference between the mother and the sister masses. (3rd argument
!     is '-1', because this alternative mass is SMALLER than the
!     resonance mass).
                  masslow(itree(k,j))=min(masslow(itree(k,j)), &
  &                        max(masslow(j)-xm(itree(3-k,j)),0d0)) ! mass difference
                  widthlow(itree(k,j))=max(widthlow(itree(k,j)), &
  &                        widthlow(j)+xw(itree(3-k,j))) ! sum of widths
                  if (pwidth(itree(k,j),iconf).eq.0d0 .or. &
  &                        masslow(itree(k,j)).ge.pmass(itree(k,j) &
  &                        ,iconf)) cycle
                  cBW_FKS_mass(iFKS,-1,itree(k,j),ichan)= &
  &                        masslow(itree(k,j))
                  cBW_FKS_width(iFKS,-1,itree(k,j),ichan)= &
  &                        widthlow(itree(k,j))
               enddo
            enddo
         else
!     Normal Breit-Wigner
            cBW_FKS(iFKS,i,ichan)=0
         endif
      enddo
! loop over t-channel to make sure that s-hat is consistent with sum of
! s-channel masses
      if (t_channel.ne.0) then
         sum_all_s=0d0
         do i=t_channel,-(nexternal-3),-1
! Breit-wigner can never go on-shell:
            if (itree(2,i).gt.0) cycle
            if ( pmass(itree(2,i),iconf).gt.sqrt(stot) .and. &
  &                  pwidth(itree(2,i),iconf).gt.0d0) then
               cBW_FKS(iFKS,itree(2,i),ichan)=2
            endif
!     s-channel is always 2nd argument of itree, sum it to sum_all_s
            sum_all_s=sum_all_s+xm(itree(2,i))
         enddo
         if (sum_all_s.gt.sqrt(stot)) then
!     conflicting BWs: set all s-channels as conflicting
            do i=t_channel,-(nexternal-3),-1
               if (itree(2,i).gt.0) cycle
               if (cBW_FKS(iFKS,itree(2,i),ichan).ne.2) then
                  cBW_FKS(iFKS,itree(2,i),ichan)=1
                  cBW_FKS_mass(iFKS,-1,itree(2,i),ichan)=sqrt(stot)/2d0
                  cBW_FKS_width(iFKS,-1,itree(2,i),ichan)=xw(itree(2,i))
               endif
            enddo
         endif
      endif


! Conflicting BW's determined. They are saved in cBW_FKS
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! If the lower bound found here is smaller than the hard bound,
! simply set the soft bound equal to the hard bound.
      taumin_s(iFKS,ichan)= &
  &            max(taumin_j(iFKS,ichan),taumin_s(iFKS,ichan))
!
! For the bound, we have to square and divide by stot.
      tau_lower_bound_resonance=taumin_s(iFKS,ichan)**2/stot
!
      if (j_fks.gt.nincoming) then
         write (*,'(a7,1x,i3,1x,i5,1x,a1,3(e12.5,1x))') 'tau_min' &
  &               ,iFKS,ichan,':',taumin(iFKS,ichan),taumin_j(iFKS &
  &               ,ichan),taumin_s(iFKS,ichan)
      else
         write (*,'(a7,1x,i3,1x,i5,1x,a1,e12.5,1x,a13,e12.5,1x)') &
  &               'tau_min',iFKS,ichan,':',taumin(iFKS,ichan) &
  &               ,'     --      ',taumin_s(iFKS,ichan)
      endif
   enddo
endif
tau_Born_lower_bound=taumin(nFKSprocess,ichan)**2/stot
tau_lower_bound=taumin_j(nFKSprocess,ichan)**2/stot
tau_lower_bound_resonance=taumin_s(nFKSprocess,ichan)**2/stot
do i=-nexternal,-1
   cBW(i)=cBW_FKS(nFKSprocess,i,ichan)
   cBW_level(i)=cBW_FKS_level(nFKSprocess,i,ichan)
   do j=-1,1,2
      cBW_mass(j,i)=cBW_FKS_mass(nFKSprocess,j,i,ichan)
      cBW_width(j,i)=cBW_FKS_width(nFKSprocess,j,i,ichan)
   enddo
enddo
do i=-nexternal,nexternal
   s_mass(i)=s_mass_FKS(nFKSprocess,i,ichan)
enddo
cBW_level_max=cBW_FKS_level_max(nFKSprocess,ichan)
end subroutine set_tau_min_impl


subroutine schan_order_impl(ns_channel, order, itree)
implicit none
integer, intent(in) :: ns_channel
integer, intent(inout) :: order(-nexternal:0)
integer, intent(in) :: itree(2, -max_branch:-1)
double precision :: rnd
integer :: i, j, ipos, npos, pos(nexternal)
logical :: done(-nexternal:nexternal)

call validate_process_dimensions()
call initialize_setcuts_state()
if (ns_channel < 0 .or. ns_channel > nexternal) then
   call fail_setcuts('sChan_order received an invalid channel count')
end if
if (.not. new_point) then
   do j=-1,-ns_channel,-1
      order(j)=saved_schan_order(j)
   enddo
   return
endif
do i=-ns_channel,0
   done(i)=.false.
enddo
do i=1,nexternal
   done(i)=.true.
enddo
do j=-1,-ns_channel,-1
   npos=0
   do i=-1,-ns_channel,-1
      if((.not. done(i)) .and. &
  &            done(itree(1,i))  .and. done(itree(2,i))) then
         npos=npos+1
         pos(npos)=i
      endif
   enddo
   if (npos.gt.1) then
      rnd = random_unit_interval(iconfig)
      ipos=min(int(rnd*npos)+1,npos)
      saved_schan_order(j)=pos(ipos)
      done(pos(ipos))=.true.
   elseif (npos.eq.1) then
      saved_schan_order(j)=pos(npos)
      done(pos(npos))=.true.
   else
      write (*,*) 'ERROR in sChan_order',npos
   endif
   order(j)=saved_schan_order(j)
enddo
new_point=.false.
end subroutine schan_order_impl


subroutine initialize_setcuts_state()
  implicit none

  call validate_process_dimensions()
  if (setcuts_state_initialized) then
    call validate_setcuts_state()
    return
  end if

  allocate(taumin(fks_configs, maxchannels))
  allocate(taumin_s(fks_configs, maxchannels))
  allocate(taumin_j(fks_configs, maxchannels))
  allocate(mass_min(-nexternal:nexternal, maxchannels))
  allocate(cbw_fks_level_max(fks_configs, maxchannels))
  allocate(cbw_fks(fks_configs, -nexternal:-1, maxchannels))
  allocate(cbw_fks_level(fks_configs, -nexternal:-1, maxchannels))
  allocate(cbw_fks_mass(fks_configs, -1:1, -nexternal:-1, &
       maxchannels))
  allocate(cbw_fks_width(fks_configs, -1:1, -nexternal:-1, &
       maxchannels))
  allocate(s_mass_fks(fks_configs, -nexternal:nexternal, maxchannels))
  allocate(firsttime_chans(maxchannels))
  allocate(saved_schan_order(-nexternal:0))

  taumin = 0d0
  taumin_s = 0d0
  taumin_j = 0d0
  mass_min = 0d0
  cbw_fks_level_max = 0
  cbw_fks = 0
  cbw_fks_level = 0
  cbw_fks_mass = 0d0
  cbw_fks_width = 0d0
  s_mass_fks = 0d0
  firsttime_chans = .true.
  saved_schan_order = 0
  stot = 0d0
  setcuts_state_initialized = .true.
  call validate_setcuts_state()
end subroutine initialize_setcuts_state




subroutine validate_setcuts_state()
  implicit none

  if (.not. setcuts_state_initialized) then
    call fail_setcuts('module state is not initialized')
  end if
  if (.not. allocated(taumin) .or. .not. allocated(taumin_s) .or. &
      .not. allocated(taumin_j) .or. .not. allocated(mass_min) .or. &
      .not. allocated(cbw_fks_level_max) .or. &
      .not. allocated(cbw_fks) .or. .not. allocated(cbw_fks_level) .or. &
      .not. allocated(cbw_fks_mass) .or. &
      .not. allocated(cbw_fks_width) .or. &
      .not. allocated(s_mass_fks) .or. &
      .not. allocated(firsttime_chans) .or. &
      .not. allocated(saved_schan_order)) then
    call fail_setcuts('module storage is incomplete')
  end if
  if (size(taumin, 1) /= fks_configs .or. &
      size(taumin, 2) /= maxchannels .or. &
      lbound(mass_min, 1) /= -nexternal .or. &
      ubound(mass_min, 1) /= nexternal .or. &
      size(mass_min, 2) /= maxchannels .or. &
      size(cbw_fks, 1) /= fks_configs .or. &
      lbound(cbw_fks, 2) /= -nexternal .or. &
      ubound(cbw_fks, 2) /= -1 .or. &
      size(cbw_fks, 3) /= maxchannels .or. &
      lbound(cbw_fks_mass, 2) /= -1 .or. &
      ubound(cbw_fks_mass, 2) /= 1 .or. &
      lbound(s_mass_fks, 2) /= -nexternal .or. &
      ubound(s_mass_fks, 2) /= nexternal .or. &
      size(firsttime_chans) /= maxchannels .or. &
      lbound(saved_schan_order, 1) /= -nexternal .or. &
      ubound(saved_schan_order, 1) /= 0) then
    call fail_setcuts('module storage has inconsistent bounds')
  end if
end subroutine validate_setcuts_state


subroutine validate_setcuts_inputs(pmass, idup, etmin, etmax, mxxmin, &
     is_a_j_compat, is_a_lp_compat, is_a_lm_compat, is_a_ph_compat)
  implicit none
  double precision, intent(in) :: pmass(:)
  integer, intent(in) :: idup(:, :)
  double precision, intent(in) :: etmin(nincoming + 1:)
  double precision, intent(in) :: etmax(nincoming + 1:)
  double precision, intent(in) :: mxxmin(nincoming + 1:, nincoming + 1:)
  logical, intent(inout) :: is_a_j_compat(:), is_a_lp_compat(:)
  logical, intent(inout) :: is_a_lm_compat(:), is_a_ph_compat(:)

  if (size(pmass) /= nexternal .or. size(idup, 1) /= nexternal .or. &
      size(idup, 2) < 1) then
    call fail_setcuts('setcuts received inconsistent process data')
  end if
  if (size(etmin) /= nexternal - nincoming - 1 .or. &
      size(etmax) /= nexternal - nincoming - 1 .or. &
      size(mxxmin, 1) /= nexternal - nincoming - 1 .or. &
      size(mxxmin, 2) /= nexternal - nincoming - 1) then
    call fail_setcuts('setcuts received inconsistent cut storage')
  end if
  if (size(is_a_j_compat) /= nexternal .or. &
      size(is_a_lp_compat) /= nexternal .or. &
      size(is_a_lm_compat) /= nexternal .or. &
      size(is_a_ph_compat) /= nexternal) then
    call fail_setcuts('setcuts received inconsistent classification storage')
  end if
end subroutine validate_setcuts_inputs


subroutine validate_tau_inputs(pmass, pwidth, itree, iconf, &
     nfksprocess, idup, emass, etmin, etmax, mxxmin, is_a_j_compat, &
     is_a_lp_compat, is_a_lm_compat, is_a_ph_compat, cbw_mass, &
     cbw_width, cbw, cbw_level, s_mass)
  implicit none
  double precision, intent(in) :: pmass(-nexternal:, :)
  double precision, intent(in) :: pwidth(-nexternal:, :)
  integer, intent(in) :: itree(2, -max_branch:-1), iconf, nfksprocess
  integer, intent(in) :: idup(:, :)
  double precision, intent(in) :: emass(:)
  double precision, intent(in) :: etmin(nincoming + 1:)
  double precision, intent(in) :: etmax(nincoming + 1:)
  double precision, intent(in) :: mxxmin(nincoming + 1:, nincoming + 1:)
  logical, intent(in) :: is_a_j_compat(:), is_a_lp_compat(:)
  logical, intent(in) :: is_a_lm_compat(:), is_a_ph_compat(:)
  double precision, intent(inout) :: cbw_mass(-1:, -nexternal:)
  double precision, intent(inout) :: cbw_width(-1:, -nexternal:)
  integer, intent(out) :: cbw(-nexternal:), cbw_level(-nexternal:)
  double precision, intent(out) :: s_mass(-nexternal:)

  if (ichan < 1 .or. ichan > maxchannels) then
    call fail_setcuts('set_tau_min received an invalid MINT channel')
  end if
  if (nfksprocess < 1 .or. nfksprocess > fks_configs) then
    call fail_setcuts('set_tau_min received an invalid FKS configuration')
  end if
  if (iconf < 1 .or. iconf > lmaxconfigs) then
    call fail_setcuts('set_tau_min received an invalid phase-space config')
  end if
  if (size(pmass, 1) /= nexternal + 1 .or. &
      size(pmass, 2) /= lmaxconfigs .or. &
      any(shape(pwidth) /= shape(pmass))) then
    call fail_setcuts('set_tau_min received inconsistent propagator data')
  end if
  if (size(itree, 2) /= max_branch .or. &
      size(idup, 1) /= nexternal .or. size(idup, 2) /= maxproc .or. &
      size(emass) /= nexternal) then
    call fail_setcuts('set_tau_min received inconsistent process data')
  end if
  if (size(etmin) /= nexternal - nincoming - 1 .or. &
      size(etmax) /= nexternal - nincoming - 1 .or. &
      size(mxxmin, 1) /= nexternal - nincoming - 1 .or. &
      size(mxxmin, 2) /= nexternal - nincoming - 1) then
    call fail_setcuts('set_tau_min received inconsistent cut data')
  end if
  if (size(is_a_j_compat) /= nexternal .or. &
      size(is_a_lp_compat) /= nexternal .or. &
      size(is_a_lm_compat) /= nexternal .or. &
      size(is_a_ph_compat) /= nexternal) then
    call fail_setcuts('set_tau_min received inconsistent classifications')
  end if
  if (size(cbw_mass, 1) /= 3 .or. size(cbw_mass, 2) /= nexternal .or. &
      any(shape(cbw_width) /= shape(cbw_mass)) .or. &
      size(cbw) /= nexternal .or. size(cbw_level) /= nexternal .or. &
      size(s_mass) /= 2 * nexternal + 1) then
    call fail_setcuts('set_tau_min received inconsistent output storage')
  end if
end subroutine validate_tau_inputs


subroutine fail_setcuts(message)
  implicit none
  character(len=*), intent(in) :: message

  write(*, '(a)') 'ERROR in setcuts_module: ' // trim(message)
  stop 1
end subroutine fail_setcuts

end module setcuts_module
