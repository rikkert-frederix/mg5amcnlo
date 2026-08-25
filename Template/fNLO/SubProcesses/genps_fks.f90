module genps_fks
use boostwdir2_module, only: boostwdir2, boostwdir2_in_place
use process_dimensions, only: nexternal, nincoming, max_particles, &
  max_branch, lmaxconfigs, fks_configs, validate_process_dimensions
use run_state
use timing_state, only: tGenPS
use kin_functions_module, only: dot => dot_impl, rho => rho_impl
implicit none
private
public :: generate_momenta
public :: initialize_genps_fks_state

! Python writes the dimensions of these historical COMMON blocks for each
! subprocess.  The fixed-form bridge owns that generated storage; these
! pointers let the free-form module use the same live objects without
! embedding generated include files or automatic objects in COMMON.
double precision, pointer :: config_mass(:, :, :) => null()
double precision, pointer :: config_width(:, :, :) => null()
integer, pointer :: config_forest(:, :, :, :) => null()
integer, pointer :: config_tree(:, :) => null()
integer, pointer :: config_index => null()

double precision, pointer :: cnt_momenta(:, :, :) => null()
double precision, pointer :: cnt_weight(:) => null()
double precision, pointer :: cnt_psweight(:) => null()
double precision, pointer :: cnt_jacobian(:) => null()

integer, pointer :: born_tree(:, :) => null()
integer, pointer :: born_ns_channel => null()
integer, pointer :: born_nt_channel => null()
integer, pointer :: born_onebody_index => null()
integer, pointer :: born_nbranch => null()
logical, pointer :: born_one_body => null()

double precision, pointer :: born_momenta(:, :) => null()
double precision, pointer :: born_lab_momenta(:, :) => null()
double precision, pointer :: born_coll_momenta(:, :) => null()
double precision, pointer :: born_norad_momenta(:, :) => null()
double precision, pointer :: event_momenta(:, :) => null()

double precision, pointer :: cbw_mass_state(:, :) => null()
double precision, pointer :: cbw_width_state(:, :) => null()
integer, pointer :: cbw_level_max_state => null()
integer, pointer :: cbw_state(:) => null()
integer, pointer :: cbw_level_state(:) => null()

double precision, pointer :: particle_masses(:) => null()
double precision, pointer :: schannel_masses(:) => null()
double precision, allocatable :: saved_particle_masses(:)
integer :: saved_configuration = 0
logical :: first_configuration = .true.
double precision :: saved_stot = 0d0
double precision :: saved_initial_mass = 0d0
double precision :: saved_final_mass = 0d0
double precision :: saved_fks_mass = 0d0
double precision :: massive_xjac_cache = 1d0
logical :: genps_state_initialized = .false.

contains

subroutine initialize_genps_fks_state(config_mass_in, config_width_in, &
  config_forest_in, config_tree_in, config_index_in, cnt_momenta_in, &
  cnt_weight_in, cnt_psweight_in, cnt_jacobian_in, born_tree_in, &
  born_ns_in, born_nt_in, born_onebody_in, born_nbranch_in, &
  born_one_body_in, born_momenta_in, born_lab_momenta_in, &
  born_coll_momenta_in, born_norad_momenta_in, event_momenta_in, &
  cbw_mass_in, cbw_width_in, cbw_level_max_in, cbw_in, cbw_level_in, &
  particle_masses_in, schannel_masses_in)
implicit none
double precision, target, intent(inout) :: &
  config_mass_in(-nexternal:, 1:, 0:)
double precision, target, intent(inout) :: &
  config_width_in(-nexternal:, 1:, 0:)
integer, target, intent(inout) :: &
  config_forest_in(1:, -max_branch:, 1:, 0:)
integer, target, intent(inout) :: config_tree_in(1:, -max_branch:)
integer, target, intent(inout) :: config_index_in
double precision, target, intent(inout) :: &
  cnt_momenta_in(0:, 1:, -2:)
double precision, target, intent(inout) :: cnt_weight_in(-2:)
double precision, target, intent(inout) :: cnt_psweight_in(-2:)
double precision, target, intent(inout) :: cnt_jacobian_in(-2:)
integer, target, intent(inout) :: born_tree_in(1:, -max_branch:)
integer, target, intent(inout) :: born_ns_in, born_nt_in
integer, target, intent(inout) :: born_onebody_in, born_nbranch_in
logical, target, intent(inout) :: born_one_body_in
double precision, target, intent(inout) :: born_momenta_in(0:, 1:)
double precision, target, intent(inout) :: born_lab_momenta_in(0:, 1:)
double precision, target, intent(inout) :: born_coll_momenta_in(0:, 1:)
double precision, target, intent(inout) :: born_norad_momenta_in(0:, 1:)
double precision, target, intent(inout) :: event_momenta_in(0:, 1:)
double precision, target, intent(inout) :: &
  cbw_mass_in(-1:, -nexternal:), cbw_width_in(-1:, -nexternal:)
integer, target, intent(inout) :: cbw_level_max_in
integer, target, intent(inout) :: cbw_in(-nexternal:)
integer, target, intent(inout) :: cbw_level_in(-nexternal:)
double precision, target, intent(inout) :: particle_masses_in(1:)
double precision, target, intent(inout) :: &
  schannel_masses_in(-nexternal:)

call validate_process_dimensions()
config_mass => config_mass_in
config_width => config_width_in
config_forest => config_forest_in
config_tree => config_tree_in
config_index => config_index_in
cnt_momenta => cnt_momenta_in
cnt_weight => cnt_weight_in
cnt_psweight => cnt_psweight_in
cnt_jacobian => cnt_jacobian_in
born_tree => born_tree_in
born_ns_channel => born_ns_in
born_nt_channel => born_nt_in
born_onebody_index => born_onebody_in
born_nbranch => born_nbranch_in
born_one_body => born_one_body_in
born_momenta => born_momenta_in
born_lab_momenta => born_lab_momenta_in
born_coll_momenta => born_coll_momenta_in
born_norad_momenta => born_norad_momenta_in
event_momenta => event_momenta_in
cbw_mass_state => cbw_mass_in
cbw_width_state => cbw_width_in
cbw_level_max_state => cbw_level_max_in
cbw_state => cbw_in
cbw_level_state => cbw_level_in
particle_masses => particle_masses_in
schannel_masses => schannel_masses_in

call validate_bound_genps_state()

if (.not. allocated(saved_particle_masses)) then
  allocate(saved_particle_masses(-max_branch:max_particles))
else if (lbound(saved_particle_masses, 1) /= -max_branch .or. &
    ubound(saved_particle_masses, 1) /= max_particles) then
  call fail_genps_state('saved mass storage has inconsistent bounds')
end if
genps_state_initialized = .true.
end subroutine initialize_genps_fks_state

subroutine generate_momenta(ndim,iconfig,wgt,x,p)
implicit none
integer ndim,iconfig
double precision wgt,x(99),p(0:3,nexternal)
real :: tBefore, tAfter
double precision qmass(-nexternal:0),qwidth(-nexternal:0),jac
integer i,j
double precision zero
parameter (zero=0d0)
integer iconfig0
common/ciconfig0/iconfig0
integer            this_config
common/to_mconfigs/this_config
!
call require_genps_state()
call cpu_time(tBefore)
this_config=iconfig
config_index=iconfig
iconfig0=iconfig
do i=-max_branch,-1
do j=1,2
config_tree(j,i)=config_forest(j,i,iconfig,0)
enddo
enddo

do i=-nexternal,0
qmass(i)=config_mass(i,iconfig,0)
qwidth(i)=config_width(i,iconfig,0)
enddo
!
call generate_momenta_conf_wrapper(ndim,jac,x,config_tree,qmass,qwidth,p)
! If the input weight 'wgt' to this subroutine was not equal to one,
! make sure we update all the (counter-event) jacobians and return also
! the updated wgt (i.e. the jacobian for the event)
do i=-2,2
cnt_jacobian(i)=cnt_jacobian(i)*wgt
enddo
wgt=wgt*jac
!
call cpu_time(tAfter)
tGenPS=tGenPS+(tAfter-tBefore)
return
end subroutine generate_momenta

subroutine generate_momenta_conf_wrapper(nndim,jac,x,itree,qmass &
  & ,qwidth,p)
!     Standard FKS generation: establish the tau bound and generate the
!     ordinary event and counterevent contexts in one call.
implicit none
integer nndim
double precision jac,x(99),p(0:3,nexternal)
integer itree(2,-max_branch:-1)
double precision qmass(-nexternal:0),qwidth(-nexternal:0)

call set_tau_min()
call generate_momenta_conf(nndim,jac,x,itree,qmass,qwidth,p)
return
end subroutine generate_momenta_conf_wrapper


subroutine generate_tau_y_wrapper( &
  & qmass,qwidth,totmass,stot,rndx,tau_born,ycm_born,ycmhat,xjac)
! generates tau and y, calling the functions that correpsond to the
! case at hand
implicit none

double precision qmass(-nexternal:0),qwidth(-nexternal:0)
double precision totmass, stot
double precision rndx(2)
double precision tau_born, ycm_born, ycmhat, xjac
!
integer ndim_dummy
double precision fksmass

integer i_fks,j_fks
common/fks_indices/i_fks,j_fks

logical softtest,colltest
common/sctests/softtest,colltest

ndim_dummy=-1 ! this is actually not used anymore

if (abs(lpp(1)).gt.1 .or. abs(lpp(2)).gt.1) then
write(*,*) 'The fNLO template supports only lpp=0,+1,-1', lpp
stop 1
endif

if (abs(lpp(1)).ne.abs(lpp(2))) then
write(*,*) 'Different beams not implemented', lpp
stop 1
endif

if (abs(lpp(1)).ge.1 .and. abs(lpp(2)).ge.1 .and. &
  & .not.(softtest.or.colltest)) then
! x(ndim-1) -> tau_cnt(0); x(ndim) -> ycm_cnt(0)
!  rndx(1) -> tau; rndx(2) -> ycm
if (born_one_body) then
! tau is fixed by the mass of the final state particle
call compute_tau_one_body(totmass,stot,tau_born,xjac)
else
if(born_nt_channel.eq.0 .and. &
  & qwidth(-born_ns_channel-1).ne.0.d0 .and. &
  & cbw_state(-born_ns_channel-1).ne.2)then
! Generate tau according to a Breit-Wiger function
call generate_tau_BW(stot,ndim_dummy,rndx(1),qmass( &
  & -born_ns_channel-1),qwidth(-born_ns_channel-1), &
  & cbw_state(-born_ns_channel-1), &
  & cbw_mass_state(-1:1, -born_ns_channel-1), &
  & cbw_width_state(-1:1, -born_ns_channel-1),tau_born,xjac)
else
!     not a Breit Wigner
call generate_tau(stot,ndim_dummy,rndx(1),tau_born,xjac)
endif
endif

! Generate the rapditity of the Born system
call generate_y(tau_born,rndx(2),ycm_born,ycmhat,xjac)
elseif (abs(lpp(1)).ge.1 .and. &
  & .not.(softtest.or.colltest)) then
write(*,*)'Option x1 not implemented in one_tree'
stop
elseif (abs(lpp(2)).ge.1 .and. &
  & .not.(softtest.or.colltest)) then
write(*,*)'Option x2 not implemented in one_tree'
stop
else
! No PDFs (also use fixed energy when performing tests)
call compute_tau_y_epem(j_fks,born_one_body,totmass,stot, &
  & tau_born,ycm_born,ycmhat)
if (j_fks.le.nincoming .and. .not.(softtest.or.colltest)) then
write (*,*) 'Process has incoming j_fks, but fixed shat: '// &
  & 'not allowed for processes generated at NLO.'
stop 1
endif
endif

return
end subroutine generate_tau_y_wrapper


subroutine generate_momenta_born(x,shat_born,sqrtshat_born,totmass, &
  & m,s, &
  & qmass,qwidth,m_born,xpswgt0,xjac0)
! generate the momenta for the reduced born system
implicit none

double precision x(99), shat_born, sqrtshat_born, totmass
double precision S(-max_branch:max_particles),M(-max_branch:max_particles)
double precision xpswgt0, xjac0
double precision qmass(-nexternal:0),qwidth(-nexternal:0), &
  & m_born(nexternal-1)

logical pass
double precision pb(0:3,-max_branch:nexternal-1),p_born_CHECK(0:3,nexternal-1)
integer i,j

pass = .true.

! Generate the momenta for the initial state of the Born system
if(nincoming.eq.2) then
call mom2cx(sqrtshat_born,m(1),m(2),1d0,0d0,pb(0,1),pb(0,2))
else
pb(0,1)=sqrtshat_born
do i=1,2
pb(i,1)=0d0
enddo
endif
s(-born_nbranch)  = shat_born
m(-born_nbranch)  = sqrtshat_born
pb(0,-born_nbranch)= m(-born_nbranch)
pb(1,-born_nbranch)= 0d0
pb(2,-born_nbranch)= 0d0
pb(3,-born_nbranch)= 0d0
!
! Generate Born-level momenta
!
! Start by generating all the invariant masses of the s-channels
call generate_inv_mass_sch(born_ns_channel,born_tree,m,sqrtshat_born &
  & ,totmass,qwidth,qmass,cbw_state,cbw_mass_state, &
  & cbw_width_state,s,x,xjac0 &
  & ,pass)
if (.not.pass) then
xjac0=-139
return
endif
! If only s-channels, also set the p1+p2 s-channel
if (born_nt_channel .eq. 0 .and. nincoming .eq. 2) then
s(-born_nbranch+1)=s(-born_nbranch)
m(-born_nbranch+1)=m(-born_nbranch)       !Basic s-channel has s_hat
pb(0,-born_nbranch+1) = m(-born_nbranch+1)!and 0 momentum
pb(1,-born_nbranch+1) = 0d0
pb(2,-born_nbranch+1) = 0d0
pb(3,-born_nbranch+1) = 0d0
endif
!
!     Next do the T-channel branchings
!
if (born_nt_channel.ne.0) then
call generate_t_channel_branchings(born_ns_channel,born_nbranch,born_tree &
  & ,m,s,x,pb,xjac0,xpswgt0,pass)
if (.not.pass) then
xjac0=-140
return
endif
endif
!
!     Now generate momentum for all intermediate and final states
!     being careful to calculate from more massive to less massive states
!     so the last states done are the final particle states.
!
call fill_born_momenta(born_nbranch,born_nt_channel,born_one_body, &
  & born_onebody_index,x,born_tree,m,s,pb,xjac0,xpswgt0,pass)
if (.not.pass) then
xjac0=-141
return
endif
!
!  Now I have the Born momenta
!
do i=1,nexternal-1
do j=0,3
born_lab_momenta(j,i)=pb(j,i)
p_born_CHECK(j,i)=pb(j,i)
enddo
m_born(i)=m(i)
enddo
call phspncheck_born(sqrtshat_born,m_born,p_born_CHECK,pass)
if (.not.pass) then
xjac0=-142
return
endif

do i=1,nexternal-1
do j=0,3
born_momenta(j,i)=born_lab_momenta(j,i)
enddo
enddo

return
end subroutine generate_momenta_born



subroutine generate_momenta_conf(ndim,jac,x,itree,qmass,qwidth,p)
!
! x(1)...x(ndim-5) --> invariant mass & angles for the Born
! x(ndim-4) --> tau_born
! x(ndim-3) --> y_born
! x(ndim-2) --> xi_i_fks
! x(ndim-1) --> y_ij_fks
! x(ndim) --> phi_i
!
implicit none
! arguments
integer ndim
double precision jac,x(99),p(0:3,nexternal)
integer itree(2,-max_branch:-1)
double precision qmass(-nexternal:0),qwidth(-nexternal:0)
! common
!     Arguments have the following meanings:
!     -2 soft-collinear, incoming leg, - direction as in FKS paper
!     -1 collinear, incoming leg, - direction as in FKS paper
!     0 soft
!     1 collinear
!     2 soft-collinear
integer i_fks,j_fks
common/fks_indices/i_fks,j_fks
logical nocntevents
common/cnocntevents/nocntevents
integer iconfig0
common/ciconfig0/iconfig0
! local
integer i,j
double precision xjac0,S(-max_branch:max_particles) &
  & ,tau_born,ycm_born,ycmhat,xbjrk_born(2),shat_born &
  & ,sqrtshat_born,xpswgt0,m_born(nexternal-1)
logical pass
logical use_evpr
common /to_use_evpr/use_evpr
pass=.true.
do i=1,nexternal-1
if (i.lt.i_fks) then
saved_particle_masses(i)=particle_masses(i)
else
saved_particle_masses(i)=particle_masses(i+1)
endif
enddo
if( first_configuration .or. iconfig0.ne.saved_configuration ) then
if (nincoming.eq.2) then
saved_stot = 4d0*ebeam(1)*ebeam(2)
else
saved_stot=particle_masses(1)**2
endif
! Make sure have enough mass for external particles
saved_initial_mass=0d0
do i=1,nincoming
saved_initial_mass=saved_initial_mass+saved_particle_masses(i)
enddo
saved_final_mass=0d0
do i=nincoming+1,nexternal-1
saved_final_mass=saved_final_mass+saved_particle_masses(i)
enddo
saved_fks_mass=saved_final_mass
if (saved_stot .lt. max(saved_final_mass,saved_initial_mass)**2) then
write (*,*) 'Fatal error #0 in one_tree:'// &
  & 'insufficient collider energy'
stop
endif

call fill_genmom_born_commons(itree,saved_particle_masses)

first_configuration=.false.
saved_configuration=iconfig0
endif                     ! first configuration
!
xjac0=1d0
xpswgt0=1d0

! generate tau and y
call generate_tau_y_wrapper( &
  & qmass,qwidth,saved_final_mass,saved_stot,x(ndim-4:ndim-3), &
  & tau_born,ycm_born,ycmhat,xjac0)
! filter unphysical configurations
if (xjac0.lt.0d0) goto 222

! Compute Bjorken x's from tau and y
xbjrk_born(1)=sqrt(tau_born)*exp(ycm_born)
xbjrk_born(2)=sqrt(tau_born)*exp(-ycm_born)
! Compute shat and sqrt(shat)
if(.not.born_one_body)then
shat_born=tau_born*saved_stot
sqrtshat_born=sqrt(shat_born)
else
! Trivial, but prevents loss of accuracy
shat_born=saved_final_mass**2
sqrtshat_born=saved_final_mass
endif

! fNLO supports only the standard event-projection mapping.
use_evpr = .true.
call generate_momenta_born(x,shat_born,sqrtshat_born, &
  & saved_final_mass, &
  & saved_particle_masses,s, &
  & qmass,qwidth,m_born,xpswgt0,xjac0)

call generate_FKS_kinematics(x,ndim,xjac0,xpswgt0, &
  & saved_stot,shat_born,sqrtshat_born,tau_born,ycm_born,ycmhat, &
  & xbjrk_born,saved_particle_masses,m_born,jac,p,pass)

!MZ check that adding .or.xjac0<0 does not screw things up
if(.not.pass.or.xjac0.lt.0d0)goto 222
return

222 continue
!
! Born momenta have not been generated. Neither events nor counterevents exist.
! Set all to negative values and exit
jac=-222
cnt_jacobian(0)=-222
cnt_jacobian(1)=-222
cnt_jacobian(2)=-222
p(0,1)=-99
do i=-2,2
cnt_momenta(0,1,i)=-99
enddo
born_momenta(0,1)=-99
nocntevents=.true.

return
end subroutine generate_momenta_conf



subroutine generate_FKS_kinematics(x,ndim,xjac0,xpswgt0, &
  & stot,shat_born,sqrtshat_born,tau_born,ycm_born,ycmhat, &
  & xbjrk_born,m,m_born,jac,p,pass)
implicit none


double precision xjac0,xpswgt0,x(99),p(0:3,nexternal), &
  & stot,shat_born,sqrtshat_born,tau_born,ycm_born,ycmhat,jac
double precision xbjrk_born(2)
double precision M(-max_branch:max_particles),m_born(nexternal-1)
integer ndim
logical pass

integer icountevts
integer ixEi,ixyij,ixpi,imother
double precision xmrec2,m_j_fks,phi_i_fks,rat_xi,tau, &
  & xi_i_fks,y_ij_fks,xi_i_hat,xiimax,xinorm,xjac,xpswgt, &
  & ycm,xp(0:3,nexternal),xbjrk(2),p_i_fks(0:3)
integer i,j

double precision pi
parameter (pi=3.1415926535897932d0)

logical nocntevents
common/cnocntevents/nocntevents

logical nbody
common/cnbody/nbody

double precision xi_i_hat_ev,xi_i_hat_cnt(-2:2)
common /cxi_i_hat/xi_i_hat_ev,xi_i_hat_cnt

complex(kind=kind(0d0)) xij_aor
common/cxij_aor/xij_aor

integer i_fks,j_fks
common/fks_indices/i_fks,j_fks

double precision ybst_til_tolab,ybst_til_tocm,sqrtshat,shat
common/parton_cms_stuff/ybst_til_tolab,ybst_til_tocm, &
  & sqrtshat,shat

double precision xi_i_fks_ev,y_ij_fks_ev
double precision p_i_fks_ev(0:3),p_i_fks_cnt(0:3,-2:2)
common/fksvariables/xi_i_fks_ev,y_ij_fks_ev,p_i_fks_ev,p_i_fks_cnt

integer isolsign

double precision xiimax_ev
common /cxiimaxev/xiimax_ev
double precision xiimax_cnt(-2:2)
common /cxiimaxcnt/xiimax_cnt

logical fks_as_is
parameter (fks_as_is=.false.)

! check that the starting PS point is meaningful
if (xjac0.lt.0d0) then
pass = .false.
return
endif
!
! Here we start with the FKS Stuff
!
! icountevts=-100 is the event, -2 to 2 the counterevents
icountevts = -100
! Event/counterevent values stay negative until their ordinary context is
! generated.  This also forces a failure if a skipped context is consumed.
p_i_fks_ev(0)=-1.d0
xiimax_ev=-1.d0
do i=-2,2
p_i_fks_cnt(0,i)=-1.d0
xiimax_cnt(i)=-1.d0
cnt_jacobian(i)=-1.d0
enddo
! set cm stuff to values to make the program crash if not set elsewhere
ybst_til_tolab=1.d14
ybst_til_tocm=1.d14
sqrtshat=0.d0
shat=0.d0
! If the collinear counterevent is not generated, this stays zero.
xij_aor=(0.d0,0.d0)
!
! These will correspond to the vegas x's for the FKS variables xi_i,
! y_ij and phi_i (changing this also requires changing folding parameters)
ixEi=ndim-2
ixyij=ndim-1
ixpi=ndim
!
imother=min(j_fks,i_fks)
m_j_fks=particle_masses(j_fks)
!
! For final state j_fks, compute the recoil invariant mass
if (j_fks.gt.nincoming) then
call get_recoil(born_lab_momenta,imother,shat_born,xmrec2,pass)
if (.not.pass) then
xjac0=-44
return
endif
endif

! Here is the beginning of the loop over the momenta for the event and
! counter-events. This will fill the xp momenta with the event and
! counter-event momenta.
111 continue
xjac   = xjac0
xpswgt = xpswgt0
!
! Put the Born momenta in the xp momenta, making sure that the mapping
! is correct; put i_fks momenta equal to zero.
do i=1,nexternal
if(i.lt.i_fks) then
do j=0,3
xp(j,i)=born_lab_momenta(j,i)
enddo
m(i)=m_born(i)
elseif(i.eq.i_fks) then
do j=0,3
xp(j,i)=0d0
enddo
m(i)=0d0
elseif(i.ge.i_fks) then
do j=0,3
xp(j,i)=born_lab_momenta(j,i-1)
enddo
m(i)=m_born(i-1)
endif
enddo
!
! set-up phi_i_fks
!
phi_i_fks=2d0*pi*x(ixpi)
xjac=xjac*2d0*pi
! To keep track of the special phase-space region with massive j_fks
isolsign=0
!
! consider the three cases:
! case 1: j_fks is massless final state
! case 2: j_fks is massive final state
! case 3: j_fks is initial state
if (j_fks.gt.nincoming) then
shat=shat_born
sqrtshat=sqrtshat_born
tau=tau_born
ycm=ycm_born
xbjrk(1)=xbjrk_born(1)
xbjrk(2)=xbjrk_born(2)
if (m_j_fks.eq.0d0) then
isolsign=1
call generate_momenta_massless_final(icountevts,i_fks,j_fks &
  & ,born_lab_momenta(0:3,imother),shat,sqrtshat,x(ixEi),xmrec2,xp &
  & ,phi_i_fks,xiimax,xinorm,xi_i_fks,y_ij_fks,xi_i_hat &
  & ,p_i_fks,xjac,xpswgt,pass)
if (.not.pass) goto 112
elseif(m_j_fks.gt.0d0) then
call generate_momenta_massive_final(icountevts,isolsign &
  & ,rat_xi,i_fks,j_fks,born_lab_momenta(0:3,imother) &
  & ,shat,sqrtshat,m_j_fks,x(ixEi),xmrec2,xp,phi_i_fks &
  & ,xiimax,xinorm,xi_i_fks,y_ij_fks,xi_i_hat,p_i_fks,xjac &
  & ,xpswgt,pass)
if (.not.pass) goto 112
endif
elseif(j_fks.le.nincoming) then
isolsign=1
call generate_momenta_initial(icountevts,i_fks,j_fks,xbjrk_born &
  & ,tau_born,ycm_born,ycmhat,shat_born,phi_i_fks,xp,x(ixEi) &
  & ,shat,stot,sqrtshat,tau,ycm,xbjrk,p_i_fks,xiimax,xinorm &
  & ,xi_i_fks,y_ij_fks,xi_i_hat,xpswgt,xjac ,pass)
if (.not.pass) goto 112
else
write (*,*) 'Error #2 in genps_fks.f',j_fks
stop
endif
! At this point, the phase space lacks a factor xi_i_fks, which need be
! excluded in an NLO computation according to FKS, being taken into
! account elsewhere
!$$$      xpswgt=xpswgt*xi_i_fks
!
! All done, so check four-momentum conservation
if(xjac.gt.0.d0)then
call phspncheck_nocms(nexternal,sqrtshat,m,xp,pass)
if (.not.pass) then
xjac=-199
goto 112
endif
endif

call compute_flux(shat,sqrtshat,m(1),m(2),xpswgt,xjac)
!
112 continue

call fill_FKS_commons(icountevts,tau,ycm,ycm_born,shat,sqrtshat,xbjrk, &
  & xiimax,xinorm,xi_i_fks,xi_i_hat,p_i_fks,y_ij_fks,xp,p,xjac,jac)
!
if(icountevts.eq.-100)then
if( (j_fks.eq.1.or.j_fks.eq.2).and.fks_as_is )then
icountevts=-2
else
icountevts=0
endif
! skips counterevents when integrating over second fold for massive
! j_fks
if( isolsign.eq.-1 )icountevts=5
else
icountevts=icountevts+1
endif
if( (icountevts.le.2.and.m_j_fks.eq.0.d0.and.(.not.nbody)).or. &
  & (icountevts.eq.0.and.m_j_fks.eq.0.d0.and.nbody) .or. &
  & (icountevts.eq.0.and.m_j_fks.ne.0.d0) )then
goto 111 ! back to the top of the loop

elseif(icountevts.eq.5) then
! icountevts=5 only when integrating over the second fold with j_fks
! massive. The counterevents have been skipped, so make sure their
! momenta are unphysical. Born are physical if event was generated, and
! must stay so for the computation of enhancement factors.
do i=0,2
cnt_jacobian(i)=-299
cnt_momenta(0,1,i)=-99
enddo
endif

nocntevents=(cnt_jacobian(0).le.0.d0) .and. &
  & (cnt_jacobian(1).le.0.d0) .and. &
  & (cnt_jacobian(2).le.0.d0)
call xmom_compare(i_fks,j_fks,jac,cnt_jacobian,p,cnt_momenta,pass)
!
return
end subroutine generate_FKS_kinematics







subroutine compute_flux(shat,sqrtshat,m1,m2,xpswgt,xjac)
implicit none
double precision shat,sqrtshat,m1,m2,xpswgt,xjac

double precision pwgt, flux

double precision pi
parameter (pi=3.1415926535897932d0)

if(nincoming.eq.2)then
flux  = 1d0 /(2.D0*SQRT(LAMBDA(shat,m1**2,m2**2)))
else                      ! Decays
flux = 1d0/(2d0*sqrtshat)
endif
! The pi-dependent factor inserted below is due to the fact that the
! weight computed above is relevant to R_n, as defined in Kajantie's
! book, eq.(III.3.1), while we need the full n-body phase space
flux  = flux / (2d0*pi)**(3 * (nexternal-nincoming) - 4)
! This extra pi-dependent factor is due to the fact that the phase-space
! part relevant to i_fks and j_fks does contain all the pi's needed for
! the correct normalization of the phase space
flux  = flux * (2d0*pi)**3
pwgt=max(xjac*xpswgt,1d-99)
xjac = pwgt*flux
!
return
end subroutine compute_flux



subroutine fill_genmom_born_commons(itree,m)
implicit none
! arguments
integer itree(2,-max_branch:-1)
double precision M(-max_branch:max_particles)


born_tree(:,:) = itree(:,:)

born_nbranch = nexternal-3 ! nexternal is for n+1-body, while itree uses n-body

! Determine number of s- and t-channel branches, at this point it
! includes the s-channel p1+p2
born_ns_channel=1
do while(itree(1,-born_ns_channel).ne.1 .and. &
  & itree(1,-born_ns_channel).ne.2 .and. &
  & born_ns_channel.lt.born_nbranch)
m(-born_ns_channel)=0d0
born_ns_channel=born_ns_channel+1
enddo
born_ns_channel=born_ns_channel - 1
born_nt_channel=born_nbranch-born_ns_channel-1
! If no t-channles, ns_channels is one less, because we want to exclude
! the s-channel p1+p2
if (born_nt_channel .eq. 0 .and. nincoming .eq. 2) then
born_ns_channel=born_ns_channel-1
endif
! Set one_body to true if it's a 2->1 process at the Born (i.e. 2->2 for the n+1-body)
if((nexternal-nincoming).eq.2)then
born_one_body=.true.
born_onebody_index=nexternal-1
born_ns_channel=0
born_nt_channel=0
elseif((nexternal-nincoming).gt.2)then
born_one_body=.false.
else
write(*,*)'Error #1 in genps_fks.f',nexternal,nincoming
stop
endif

return
end subroutine fill_genmom_born_commons


subroutine fill_FKS_commons(icountevts,tau,ycm,ycm_born,shat,sqrtshat,xbjrk, &
  & xiimax,xinorm,xi_i_fks,xi_i_hat,p_i_fks,y_ij_fks,xp,p,xjac,jac)

implicit none
integer icountevts
double precision tau,ycm,ycm_born,shat,sqrtshat,xbjrk(2),xiimax,xinorm, &
  & xi_i_fks,xi_i_hat,p_i_fks(0:3),y_ij_fks,xp(0:3,nexternal),p(0:3,nexternal), &
  & xjac,jac

integer i,j

double precision xi_i_fks_ev,y_ij_fks_ev
double precision p_i_fks_ev(0:3),p_i_fks_cnt(0:3,-2:2)
common/fksvariables/xi_i_fks_ev,y_ij_fks_ev,p_i_fks_ev,p_i_fks_cnt

double precision xi_i_fks_cnt(-2:2)
common /cxiifkscnt/xi_i_fks_cnt

double precision xi_i_hat_ev,xi_i_hat_cnt(-2:2)
common /cxi_i_hat/xi_i_hat_ev,xi_i_hat_cnt

double precision xbjrk_ev(2),xbjrk_cnt(2,-2:2)
common/cbjorkenx/xbjrk_ev,xbjrk_cnt

double precision sqrtshat_ev,shat_ev
common/parton_cms_ev/sqrtshat_ev,shat_ev
double precision sqrtshat_cnt(-2:2),shat_cnt(-2:2)
common/parton_cms_cnt/sqrtshat_cnt,shat_cnt

double precision tau_ev,ycm_ev
common/cbjrk12_ev/tau_ev,ycm_ev
double precision tau_cnt(-2:2),ycm_cnt(-2:2)
common/cbjrk12_cnt/tau_cnt,ycm_cnt

double precision xiimax_ev
common /cxiimaxev/xiimax_ev
double precision xiimax_cnt(-2:2)
common /cxiimaxcnt/xiimax_cnt

double precision xinorm_ev
common /cxinormev/xinorm_ev
double precision xinorm_cnt(-2:2)
common /cxinormcnt/xinorm_cnt

! Catch the points for which there is no viable phase-space generation
! (still fill the common blocks with some information that is needed
! (e.g. ycm_cnt)).
if (xjac .le. 0d0 ) then
xp(0,1)=-99d0
endif
!
! Fill common blocks
if (icountevts.eq.-100) then
tau_ev=tau
ycm_ev=ycm
shat_ev=shat
sqrtshat_ev=sqrtshat
xbjrk_ev(1)=xbjrk(1)
xbjrk_ev(2)=xbjrk(2)
xiimax_ev=xiimax
xinorm_ev=xinorm
xi_i_fks_ev=xi_i_fks
xi_i_hat_ev=xi_i_hat
do i=0,3
p_i_fks_ev(i)=p_i_fks(i)
enddo
y_ij_fks_ev=y_ij_fks
do i=1,nexternal
do j=0,3
p(j,i)=xp(j,i)
event_momenta(j,i)=xp(j,i)
enddo
enddo
jac=xjac
else
tau_cnt(icountevts)=tau
! Special fix in the case the soft counter-events are not generated but
! the Born and real are. (This can happen if ptj>0 in the
! run_card). This fix is needed for set_cms_stuff to work properly.
if (icountevts.eq.0) then
ycm=ycm_born
endif
ycm_cnt(icountevts)=ycm
shat_cnt(icountevts)=shat
sqrtshat_cnt(icountevts)=sqrtshat
xbjrk_cnt(1,icountevts)=xbjrk(1)
xbjrk_cnt(2,icountevts)=xbjrk(2)
xiimax_cnt(icountevts)=xiimax
xinorm_cnt(icountevts)=xinorm
xi_i_fks_cnt(icountevts)=xi_i_fks
xi_i_hat_cnt(icountevts)=xi_i_hat
do i=0,3
p_i_fks_cnt(i,icountevts)=p_i_fks(i)
enddo
do i=1,nexternal
do j=0,3
cnt_momenta(j,i,icountevts)=xp(j,i)
enddo
enddo
cnt_jacobian(icountevts)=xjac
! the following two are obsolete, but still part of some common block:
! so give some non-physical values
cnt_weight(icountevts)=-1d99
cnt_psweight(icountevts)=-1d99
endif

return
end subroutine fill_FKS_commons




subroutine generate_momenta_massless_final(icountevts,i_fks,j_fks &
  & ,p_born_imother,shat,sqrtshat,x,xmrec2,xp,phi_i_fks,xiimax &
  & ,xinorm,xi_i_fks,y_ij_fks,xi_i_hat,p_i_fks,xjac,xpswgt &
  & ,pass)
implicit none
! arguments
integer icountevts,i_fks,j_fks
double precision shat,sqrtshat,x(2),xmrec2,xp(0:3,nexternal) &
  & ,y_ij_fks,p_born_imother(0:3),phi_i_fks,xi_i_hat
double precision xiimax,xinorm,xi_i_fks,p_i_fks(0:3),xjac,xpswgt
logical pass
! common blocks
double precision  veckn_ev,veckbarn_ev,xp0jfks
common/cgenps_fks/veckn_ev,veckbarn_ev,xp0jfks
complex(kind=kind(0d0)) xij_aor
common/cxij_aor/xij_aor
logical softtest,colltest
common/sctests/softtest,colltest
double precision xi_i_fks_fix,y_ij_fks_fix
common/cxiyfix/xi_i_fks_fix,y_ij_fks_fix
! local
integer i,j
double precision E_i_fks,x3len_i_fks,x3len_j_fks,x3len_fks_mother &
  & ,costh_i_fks,sinth_i_fks,xpifksred(0:3),th_mother_fks &
  & ,costh_mother_fks,sinth_mother_fks, phi_mother_fks &
  & ,cosphi_mother_fks,sinphi_mother_fks,recoil(0:3),sumrec &
  & ,sumrec2,betabst,gammabst,shybst,chybst,chybstmo,xdir(3) &
  & ,veckn,veckbarn,xp_mother(0:3),cosphi_i_fks &
  & ,sinphi_i_fks
complex(kind=kind(0d0)) resAoR0
! external
! parameters
double precision pi
parameter (pi=3.1415926535897932d0)
double precision xi_i_fks_matrix(-2:2)
data xi_i_fks_matrix/0.d0,-1.d8,0.d0,-1.d8,0.d0/
double precision y_ij_fks_matrix(-2:2)
data y_ij_fks_matrix/-1.d0,-1.d0,-1.d8,1.d0,1.d0/
double precision stiny,sstiny,qtiny,ctiny,cctiny
complex(kind=kind(0d0)) ximag
parameter (stiny=1d-6)
parameter (qtiny=1d-7)
parameter (ctiny=5d-7)
parameter (ximag=(0d0,1d0))
!
pass=.true.
if(softtest)then
sstiny=0.d0
else
sstiny=stiny
endif
if(colltest)then
cctiny=0.d0
else
cctiny=ctiny
endif
!
! set-up y_ij_fks
!
if( (icountevts.eq.-100.or.icountevts.eq.0) .and. &
  & ((.not.softtest) .or. &
  & (softtest.and.y_ij_fks_fix.eq.-2.d0)) .and. &
  & (.not.colltest)  )then
! importance sampling towards collinear singularity
! insert here further importance sampling towards y_ij_fks->1
y_ij_fks = -2d0*(cctiny+(1-cctiny)*x(2)**2)+1d0
elseif( (icountevts.eq.-100.or.icountevts.eq.0) .and. &
  & ((softtest.and.y_ij_fks_fix.ne.-2.d0) .or. &
  & colltest)  )then
y_ij_fks=y_ij_fks_fix
elseif(abs(icountevts).eq.2.or.abs(icountevts).eq.1)then
y_ij_fks=y_ij_fks_matrix(icountevts)
else
write(*,*)'Error #3 in genps_fks.f',icountevts
stop
endif
! importance sampling towards collinear singularity
xjac=xjac*2d0*x(2)*2d0

call getangles(p_born_imother, &
  & th_mother_fks,costh_mother_fks,sinth_mother_fks, &
  & phi_mother_fks,cosphi_mother_fks,sinphi_mother_fks)
!
! Compute maximum allowed xi_i_fks
xiimax=1-xmrec2/shat
xinorm=xiimax
!
! Define xi_i_fks
!
if( (icountevts.eq.-100.or.abs(icountevts).eq.1) .and. &
  & ((.not.colltest) .or. &
  & (colltest.and.xi_i_fks_fix.eq.-2.d0)) .and. &
  & (.not.softtest)  )then
if(icountevts.eq.-100)then
! importance sampling towards soft singularity
! insert here further importance sampling towards xi_i_hat->0
xi_i_hat=sstiny+(1-sstiny)*x(1)**2
endif
! in the case of counter events, xi_i_hat is an input to this function
xi_i_fks=xi_i_hat*xiimax
elseif( (icountevts.eq.-100.or.abs(icountevts).eq.1) .and. &
  & (colltest.and.xi_i_fks_fix.ne.-2.d0) .and. &
  & (.not.softtest)  )then
! This is to keep xi_i_hat, rather than xi_i, fixed in the tests.
if(xi_i_fks_fix.lt.xiimax)then
xi_i_fks=xi_i_fks_fix*xiimax
else
xi_i_fks=xi_i_fks_fix*xiimax
endif
elseif( (icountevts.eq.-100.or.abs(icountevts).eq.1) .and. &
  & softtest )then
if(xi_i_fks_fix.lt.1d0)then
xi_i_fks=xi_i_fks_fix*xiimax
else
xjac=-102
pass=.false.
return
endif
elseif(abs(icountevts).eq.2.or.icountevts.eq.0)then
xi_i_fks=xi_i_fks_matrix(icountevts)
else
write(*,*)'Error #4 in genps_fks.f',icountevts
stop
endif
! remove the following if no importance sampling towards soft
! singularity is performed when integrating over xi_i_hat
xjac=xjac*2d0*x(1)

! Check that xii is in the allowed range
if( icountevts.eq.-100 .or. abs(icountevts).eq.1 )then
if(xi_i_fks.gt.(1-xmrec2/shat))then
xjac=-101
pass=.false.
return
endif
elseif(icountevts.eq.0 .or. abs(icountevts).eq.2)then
! May insert here a check on whether xii<xicut, rather than doing it
! in the cross sections
continue
endif
!
! Compute costh_i_fks from xi_i_fks et al.
!
E_i_fks=xi_i_fks*sqrtshat/2d0
x3len_i_fks=E_i_fks
x3len_j_fks=(shat-xmrec2-2*sqrtshat*x3len_i_fks)/ &
  & (2*(sqrtshat-x3len_i_fks*(1-y_ij_fks)))
x3len_fks_mother=sqrt( x3len_i_fks**2+x3len_j_fks**2+ &
  & 2*x3len_i_fks*x3len_j_fks*y_ij_fks )
if(xi_i_fks.lt.qtiny)then
costh_i_fks=y_ij_fks+shat*(1-y_ij_fks**2)*xi_i_fks/ &
  & (shat-xmrec2)
if(abs(costh_i_fks).gt.1.d0)costh_i_fks=y_ij_fks
elseif(1-y_ij_fks.lt.qtiny)then
costh_i_fks=1-(shat*(1-xi_i_fks)-xmrec2)**2*(1-y_ij_fks)/ &
  & (shat-xmrec2)**2
if(abs(costh_i_fks).gt.1.d0)costh_i_fks=1.d0
else
costh_i_fks=(x3len_fks_mother**2-x3len_j_fks**2+x3len_i_fks**2) &
  & /(2*x3len_fks_mother*x3len_i_fks)
if(abs(costh_i_fks).gt.1.d0)then
if(abs(costh_i_fks).le.(1.d0+1.d-5))then
costh_i_fks=sign(1.d0,costh_i_fks)
else
write(*,*)'Fatal error #5 in one_tree', &
  & costh_i_fks,xi_i_fks,y_ij_fks,xmrec2
stop
endif
endif
endif
sinth_i_fks=sqrt(1-costh_i_fks**2)
cosphi_i_fks=cos(phi_i_fks)
sinphi_i_fks=sin(phi_i_fks)
xpifksred(1)=sinth_i_fks*cosphi_i_fks
xpifksred(2)=sinth_i_fks*sinphi_i_fks
xpifksred(3)=costh_i_fks
!
! The momentum if i_fks and j_fks
!
xp(0,i_fks)=E_i_fks
xp(0,j_fks)=sqrt(x3len_j_fks**2)
p_i_fks(0)=sqrtshat/2d0
do j=1,3
p_i_fks(j)=sqrtshat/2d0*xpifksred(j)
xp(j,i_fks)=E_i_fks*xpifksred(j)
if(j.ne.3)then
xp(j,j_fks)=-xp(j,i_fks)
else
xp(j,j_fks)=x3len_fks_mother-xp(j,i_fks)
endif
enddo
!
call rotate_invar(xp(0,i_fks),xp(0,i_fks), &
  & costh_mother_fks,sinth_mother_fks, &
  & cosphi_mother_fks,sinphi_mother_fks)
call rotate_invar(xp(0,j_fks),xp(0,j_fks), &
  & costh_mother_fks,sinth_mother_fks, &
  & cosphi_mother_fks,sinphi_mother_fks)
call rotate_invar(p_i_fks,p_i_fks, &
  & costh_mother_fks,sinth_mother_fks, &
  & cosphi_mother_fks,sinphi_mother_fks)
!
! Now the xp four vectors of all partons except i_fks and j_fks will be
! boosted along the direction of the mother; start by redefining the
! mother four momenta
do i=0,3
xp_mother(i)=xp(i,i_fks)+xp(i,j_fks)
if (nincoming.eq.2) then
recoil(i)=xp(i,1)+xp(i,2)-xp_mother(i)
else
recoil(i)=xp(i,1)-xp_mother(i)
endif
enddo
sumrec=recoil(0)+rho(recoil)
sumrec2=sumrec**2
betabst=-(shat-sumrec2)/(shat+sumrec2)
gammabst=1/sqrt(1-betabst**2)
shybst=-(shat-sumrec2)/(2*sumrec*sqrtshat)
chybst=(shat+sumrec2)/(2*sumrec*sqrtshat)
! cosh(y) is very often close to one, so define cosh(y)-1 as well
chybstmo=(sqrtshat-sumrec)**2/(2*sumrec*sqrtshat)
do j=1,3
xdir(j)=xp_mother(j)/x3len_fks_mother
enddo
! Perform the boost here
do i=nincoming+1,nexternal
if(i.ne.i_fks.and.i.ne.j_fks.and.shybst.ne.0.d0) &
  & call boostwdir2_in_place(chybst,shybst,chybstmo,xdir, &
  & xp(0,i))
enddo
!
! Collinear limit of <ij>/[ij]. See innerp3.m.
if( ( icountevts.eq.-100 .or. &
  & (icountevts.eq.1.and.xij_aor.eq.0) ) )then
resAoR0=-exp( 2*ximag*(phi_mother_fks+phi_i_fks) )
! The term O(srt(1-y)) is formally correct but may be numerically large
! Set it to zero
!$$$          resAoR5=-ximag*sqrt(2.d0)*
!$$$       &          sinphi_i_fks*tan(th_mother_fks/2.d0)*
!$$$       &          exp( 2*ximag*(phi_mother_fks+phi_i_fks) )
!$$$          xij_aor=resAoR0+resAoR5*sqrt(1-y_ij_fks)
xij_aor=resAoR0
endif
!
! Phase-space factor for (xii,yij,phii)
veckn=rho(xp(0,j_fks))
veckbarn=rho(p_born_imother)
!
! Store event-kinematics quantities.
if(icountevts.eq.-100)then
veckn_ev=veckn
veckbarn_ev=veckbarn
xp0jfks=xp(0,j_fks)
endif
!
xpswgt=xpswgt*2*shat/(4*pi)**3*veckn/veckbarn/ &
  & ( 2-xi_i_fks*(1-xp(0,j_fks)/veckn*y_ij_fks) )
xpswgt=abs(xpswgt)
return
end subroutine generate_momenta_massless_final

subroutine generate_momenta_massive_final(icountevts,isolsign &
  & ,rat_xi,i_fks,j_fks,p_born_imother,shat &
  & ,sqrtshat,m_j_fks,x,xmrec2,xp,phi_i_fks,xiimax,xinorm &
  & ,xi_i_fks,y_ij_fks,xi_i_hat,p_i_fks,xjac,xpswgt,pass)
implicit none
! arguments
integer icountevts,i_fks,j_fks,isolsign
double precision shat,sqrtshat,x(2),xmrec2,xp(0:3,nexternal) &
  & ,y_ij_fks,p_born_imother(0:3),m_j_fks,phi_i_fks,xi_i_hat
double precision xiimax,xinorm,xi_i_fks,p_i_fks(0:3),xjac,xpswgt
logical pass
! common blocks
double precision  veckn_ev,veckbarn_ev,xp0jfks
common/cgenps_fks/veckn_ev,veckbarn_ev,xp0jfks
logical softtest,colltest
common/sctests/softtest,colltest
double precision xi_i_fks_fix,y_ij_fks_fix
common/cxiyfix/xi_i_fks_fix,y_ij_fks_fix
! local
integer i,j
double precision xmj,xmj2,xmjhat,xmhat,xim,cffA2,cffB2,cffC2 &
  & ,cffDEL2,xiBm,ximax,xirplus,xirminus,rat_xi,xitmp1 &
  & ,E_i_fks,x3len_i_fks,b2m4ac,x3len_j_fks_num,x3len_j_fks_den &
  & ,x3len_j_fks,x3len_fks_mother,costh_i_fks,sinth_i_fks &
  & ,xpifksred(0:3),recoil(0:3),xp_mother(0:3),sumrec,expybst &
  & ,shybst,chybst,chybstmo,xdir(3),veckn,veckbarn ,cosphi_i_fks &
  & ,sinphi_i_fks,cosphi_mother_fks,costh_mother_fks &
  & ,phi_mother_fks,sinphi_mother_fks,th_mother_fks,xitmp2 &
  & ,sinth_mother_fks
! external
! parameters
double precision pi
parameter (pi=3.1415926535897932d0)
double precision xi_i_fks_matrix(-2:2)
data xi_i_fks_matrix/0.d0,-1.d8,0.d0,-1.d8,0.d0/
double precision y_ij_fks_matrix(-2:2)
data y_ij_fks_matrix/-1.d0,-1.d0,-1.d8,1.d0,1.d0/
double precision stiny,sstiny,qtiny,ctiny,cctiny
parameter (stiny=1d-6)
parameter (qtiny=1d-7)
parameter (ctiny=5d-7)
!
if(colltest .or. &
  & abs(icountevts).eq.1.or.abs(icountevts).eq.2)then
write(*,*)'Error #5 in genps_fks.f:'
write(*,*) &
  & 'This parametrization cannot be used in FS coll limits'
stop
endif
!
pass=.true.
if(softtest)then
sstiny=0.d0
else
sstiny=stiny
endif
if(colltest)then
cctiny=0.d0
else
cctiny=ctiny
endif
!
! set-up y_ij_fks
!
if( (icountevts.eq.-100.or.icountevts.eq.0) .and. &
  & ((.not.softtest) .or. &
  & (softtest.and.y_ij_fks_fix.eq.-2.d0)) .and. &
  & (.not.colltest)  )then
! importance sampling towards collinear singularity
! insert here further importance sampling towards y_ij_fks->1
y_ij_fks = -2d0*(cctiny+(1-cctiny)*x(2)**2)+1d0
elseif( (icountevts.eq.-100.or.icountevts.eq.0) .and. &
  & ((softtest.and.y_ij_fks_fix.ne.-2.d0) .or. &
  & colltest)  )then
y_ij_fks=y_ij_fks_fix
elseif(abs(icountevts).eq.2.or.abs(icountevts).eq.1)then
y_ij_fks=y_ij_fks_matrix(icountevts)
else
write(*,*)'Error #6 in genps_fks.f',icountevts
stop
endif
! importance sampling towards collinear singularity
xjac=xjac*2d0*x(2)*2d0

call getangles(p_born_imother, &
  & th_mother_fks,costh_mother_fks,sinth_mother_fks, &
  & phi_mother_fks,cosphi_mother_fks,sinphi_mother_fks)
!
! Compute the maximum allowed xi_i_fks
!
xmj=m_j_fks
xmj2=xmj**2
xmjhat=xmj/sqrtshat
xmhat=sqrt(xmrec2)/sqrtshat
xim=(1-xmhat**2-2*xmjhat+xmjhat**2)/(1-xmjhat)
cffA2=1-xmjhat**2*(1-y_ij_fks**2)
cffB2=-2*(1-xmhat**2-xmjhat**2)
cffC2=(1-(xmhat-xmjhat)**2)*(1-(xmhat+xmjhat)**2)
cffDEL2=cffB2**2-4*cffA2*cffC2
xiBm=(-cffB2-sqrt(cffDEL2))/(2*cffA2)
ximax=1-(xmhat+xmjhat)**2
if(xiBm.lt.(xim-1.d-8).or.xim.lt.0.d0.or.xiBm.lt.0.d0.or. &
  & xiBm.gt.(ximax+1.d-8).or.ximax.gt.1.or.ximax.lt.0.d0)then
write(*,*)'WARNING #4 in one_tree',xim,xiBm,ximax
xjac=-104d0
pass=.false.
return
endif
if(y_ij_fks.ge.0.d0)then
xirplus=xim
xirminus=0.d0
else
xirplus=xiBm
xirminus=xiBm-xim
endif
xiimax=xirplus
xinorm=xirplus+xirminus
rat_xi=xiimax/xinorm
!
! Generate xi_i_fks
!
if( icountevts.eq.-100 .and. &
  & ((.not.colltest) .or. &
  & (colltest.and.xi_i_fks_fix.eq.-2.d0)) .and. &
  & (.not.softtest)  )then
massive_xjac_cache=1.d0
xitmp1=x(1)
! Map regions (0,A) and (A,1) in xitmp1 onto regions (0,rat_xi) and (rat_xi,1)
! in xi_i_hat respectively. The parameter A is free, but it appears to be
! convenient to choose A=rat_xi
if(xitmp1.le.rat_xi)then
xitmp1=xitmp1/rat_xi
massive_xjac_cache=massive_xjac_cache/rat_xi
! importance sampling towards soft singularity
! insert here further importance samplings
xitmp2=sstiny+(1-sstiny)*xitmp1**2
massive_xjac_cache=massive_xjac_cache*2*xitmp1
xi_i_hat=xitmp2*rat_xi
massive_xjac_cache=massive_xjac_cache*rat_xi
xi_i_fks=xinorm*xi_i_hat
isolsign=1
else
! insert here further importance samplings
xi_i_hat=xitmp1
xi_i_fks=-xinorm*xi_i_hat+2*xiimax
isolsign=-1
endif
elseif( icountevts.eq.-100 .and. &
  & (colltest.and.xi_i_fks_fix.ne.-2.d0) .and. &
  & (.not.softtest)  )then
massive_xjac_cache=1.d0
if(xi_i_fks_fix.lt.xiimax)then
xi_i_fks=xi_i_fks_fix
else
xi_i_fks=xi_i_fks_fix*xiimax
endif
isolsign=1
elseif( (icountevts.eq.-100) .and. &
  & softtest )then
massive_xjac_cache=1.d0
if(xi_i_fks_fix.lt.xiimax)then
xi_i_fks=xi_i_fks_fix
else
xjac=-102
pass=.false.
return
endif
isolsign=1
elseif(icountevts.eq.0)then
! Keep the event Jacobian cache here for the matching counterevent.
! used for the (real-emission) event
xi_i_fks=xi_i_fks_matrix(icountevts)
isolsign=1
else
write(*,*)'Error #7 in genps_fks.f',icountevts
stop
endif
xjac=xjac*massive_xjac_cache
!
if(isolsign.eq.0)then
write(*,*)'Fatal error #11 in one_tree',isolsign
stop
endif
!
! Compute costh_i_fks
!
E_i_fks=xi_i_fks*sqrtshat/2d0
x3len_i_fks=E_i_fks
b2m4ac=xi_i_fks**2*cffA2 + xi_i_fks*cffB2 + cffC2
if(b2m4ac.le.0.d0)then
if(abs(b2m4ac).lt.1.d-3)then
b2m4ac=0.d0
else
write(*,*)'Fatal error #6 in one_tree'
write(*,*)b2m4ac,xi_i_fks,cffA2,cffB2,cffC2
write(*,*)y_ij_fks,xim,xiBm
stop
endif
endif
x3len_j_fks_num=-xi_i_fks*y_ij_fks* &
  & (1-xmhat**2+xmjhat**2-xi_i_fks) + &
  & (2-xi_i_fks)*sqrt(b2m4ac)*isolsign
x3len_j_fks_den=(2-xi_i_fks*(1-y_ij_fks))* &
  & (2-xi_i_fks*(1+y_ij_fks))
x3len_j_fks=sqrtshat*x3len_j_fks_num/x3len_j_fks_den
if(x3len_j_fks.lt.0.d0)then
write(*,*)'WARNING #7 in one_tree', &
  & x3len_j_fks_num,x3len_j_fks_den,xi_i_fks,y_ij_fks
xjac=-107d0
pass=.false.
return
endif
x3len_fks_mother=sqrt( x3len_i_fks**2+x3len_j_fks**2+ &
  & 2*x3len_i_fks*x3len_j_fks*y_ij_fks )
if(xi_i_fks.lt.qtiny)then
costh_i_fks=y_ij_fks+(1-y_ij_fks**2)*xi_i_fks/sqrt(cffC2)
if(abs(costh_i_fks).gt.1.d0)costh_i_fks=y_ij_fks
else
costh_i_fks=(x3len_fks_mother**2-x3len_j_fks**2+x3len_i_fks**2) &
  & /(2*x3len_fks_mother*x3len_i_fks)
if(abs(costh_i_fks).gt.1.d0+qtiny)then
write(*,*)'Fatal error #8 in one_tree', &
  & costh_i_fks,xi_i_fks,y_ij_fks,xmrec2
stop
elseif(abs(costh_i_fks).gt.1.d0)then
costh_i_fks = sign(1d0,costh_i_fks)
endif
endif
sinth_i_fks=sqrt(1-costh_i_fks**2)
cosphi_i_fks=cos(phi_i_fks)
sinphi_i_fks=sin(phi_i_fks)
xpifksred(1)=sinth_i_fks*cosphi_i_fks
xpifksred(2)=sinth_i_fks*sinphi_i_fks
xpifksred(3)=costh_i_fks
!
! Generate momenta for j_fks and i_fks
!
xp(0,i_fks)=E_i_fks
xp(0,j_fks)=sqrt(x3len_j_fks**2+m_j_fks**2)
p_i_fks(0)=sqrtshat/2d0
do j=1,3
p_i_fks(j)=sqrtshat/2d0*xpifksred(j)
xp(j,i_fks)=E_i_fks*xpifksred(j)
if(j.ne.3)then
xp(j,j_fks)=-xp(j,i_fks)
else
xp(j,j_fks)=x3len_fks_mother-xp(j,i_fks)
endif
enddo
!
call rotate_invar(xp(0,i_fks),xp(0,i_fks), &
  & costh_mother_fks,sinth_mother_fks, &
  & cosphi_mother_fks,sinphi_mother_fks)
call rotate_invar(xp(0,j_fks),xp(0,j_fks), &
  & costh_mother_fks,sinth_mother_fks, &
  & cosphi_mother_fks,sinphi_mother_fks)
call rotate_invar(p_i_fks,p_i_fks, &
  & costh_mother_fks,sinth_mother_fks, &
  & cosphi_mother_fks,sinphi_mother_fks)
!
! Now the xp four vectors of all partons except i_fks and j_fks will be
! boosted along the direction of the mother; start by redefining the
! mother four momenta
do i=0,3
xp_mother(i)=xp(i,i_fks)+xp(i,j_fks)
if (nincoming.eq.2) then
recoil(i)=xp(i,1)+xp(i,2)-xp_mother(i)
else
recoil(i)=xp(i,1)-xp_mother(i)
endif
enddo
!
sumrec=recoil(0)+rho(recoil)
if(xmrec2.lt.1.d-16*shat)then
expybst=sqrtshat*sumrec/(shat-xmj2)* &
  & (1+xmj2*xmrec2/(shat-xmj2)**2)
else
expybst=sumrec/(2*sqrtshat*xmrec2)* &
  & (shat+xmrec2-xmj2-shat*sqrt(cffC2))
endif
if(expybst.le.0.d0)then
write(*,*)'Fatal error #10 in one_tree',expybst
stop
endif
shybst=(expybst-1/expybst)/2.d0
chybst=(expybst+1/expybst)/2.d0
chybstmo=chybst-1.d0
!
do j=1,3
xdir(j)=xp_mother(j)/x3len_fks_mother
enddo
! Boost the momenta
do i=nincoming+1,nexternal
if(i.ne.i_fks.and.i.ne.j_fks.and.shybst.ne.0.d0) &
  & call boostwdir2_in_place(chybst,shybst,chybstmo,xdir, &
  & xp(0,i))
enddo
!
! Phase-space factor for (xii,yij,phii)
veckn=rho(xp(0,j_fks))
veckbarn=rho(p_born_imother)
!
! Store event-kinematics quantities.
if(icountevts.eq.-100)then
veckn_ev=veckn
veckbarn_ev=veckbarn
xp0jfks=xp(0,j_fks)
endif
!
xpswgt=xpswgt*2*shat/(4*pi)**3*veckn/veckbarn/ &
  & ( 2-xi_i_fks*(1-xp(0,j_fks)/veckn*y_ij_fks) )
xpswgt=abs(xpswgt)
return
end subroutine generate_momenta_massive_final


subroutine generate_momenta_initial(icountevts,i_fks,j_fks, &
  & xbjrk_born,tau_born,ycm_born,ycmhat,shat_born,phi_i_fks ,xp,x &
  & , shat,stot,sqrtshat,tau,ycm,xbjrk ,p_i_fks,xiimax,xinorm &
  & ,xi_i_fks,y_ij_fks,xi_i_hat,xpswgt ,xjac ,pass)
implicit none
! arguments
integer icountevts,i_fks,j_fks
double precision xbjrk_born(2),tau_born,ycm_born,ycmhat,shat_born &
  & ,phi_i_fks,xpswgt,xjac,xiimax,xinorm,xp(0:3,nexternal),stot &
  & ,x(2),y_ij_fks,xi_i_hat
double precision shat,sqrtshat,tau,ycm,xbjrk(2),p_i_fks(0:3)
logical pass
! common blocks
double precision tau_Born_lower_bound,tau_lower_bound_resonance &
  & ,tau_lower_bound
common/ctau_lower_bound/tau_Born_lower_bound &
  & ,tau_lower_bound_resonance,tau_lower_bound
double precision  veckn_ev,veckbarn_ev,xp0jfks
common/cgenps_fks/veckn_ev,veckbarn_ev,xp0jfks
complex(kind=kind(0d0)) xij_aor
common/cxij_aor/xij_aor
logical softtest,colltest
common/sctests/softtest,colltest
double precision xi_i_fks_fix,y_ij_fks_fix
common/cxiyfix/xi_i_fks_fix,y_ij_fks_fix
! local
integer i,j,idir
double precision yijdir,costh_i_fks,x1bar2,x2bar2,yij_sol,xi1,xi2 &
  & ,ximaxtmp,omega,bstfact,shy_tbst,chy_tbst,chy_tbstmo &
  & ,xdir_t(3),cosphi_i_fks,sinphi_i_fks,shy_lbst,chy_lbst &
  & ,encmso2,E_i_fks,sinth_i_fks,xpifksred(0:3),xi_i_fks &
  & ,xiimin,yij_upp,yij_low,y_ij_fks_upp,y_ij_fks_low
complex(kind=kind(0d0)) resAoR0

double precision omx1bar2, omx2bar2
double precision ltau_born,e2ycm_born,em2ycm_born
! external
!
! parameters
double precision pi
parameter (pi=3.1415926535897932d0)
double precision xi_i_fks_matrix(-2:2)
data xi_i_fks_matrix/0.d0,-1.d8,0.d0,-1.d8,0.d0/
double precision y_ij_fks_matrix(-2:2)
data y_ij_fks_matrix/-1.d0,-1.d0,-1.d8,1.d0,1.d0/
logical fks_as_is
parameter (fks_as_is=.false.)
complex(kind=kind(0d0)) ximag
parameter (ximag=(0d0,1d0))
double precision stiny,sstiny,qtiny,zero,ctiny,cctiny
parameter (stiny=1d-6)
parameter (qtiny=1d-7)
parameter (zero=0d0)
parameter (ctiny=5d-7)
!
pass=.true.
if(softtest)then
sstiny=0.d0
else
sstiny=stiny
endif
!
! FKS for left or right incoming parton
!
idir=0
if(.not.fks_as_is)then
if(j_fks.eq.1)then
idir=1
elseif(j_fks.eq.2)then
idir=-1
endif
else
idir=1
write(*,*)'One_tree: option not checked'
stop
endif

if (1d0-tau_born.gt.stiny) then
ltau_born = log(tau_born)
else
ltau_born = tau_born-1d0
endif
if (abs(ycm_born).gt.stiny) then
e2ycm_born = exp(2*ycm_born)
em2ycm_born = exp(-2*ycm_born)
else
e2ycm_born = 1d0 + 2*ycm_born + 2*ycm_born**2
em2ycm_born = 1d0 - 2*ycm_born + 2*ycm_born**2
endif

!
! set-up lower and upper bounds on y_ij_fks
!
if( tau_born.le.tau_lower_bound .and.ycm_born.gt. &
  & (0.5d0*ltau_born-log(tau_lower_bound)) )then
yij_upp= (tau_lower_bound+tau_born)* &
  & ( 1-e2ycm_born*tau_lower_bound ) / &
  & ( (tau_lower_bound-tau_born)* &
  & (1+e2ycm_born*tau_lower_bound) )
else
yij_upp=1.d0
endif
if( tau_born.le.tau_lower_bound .and. ycm_born.lt. &
  & (-0.5d0*ltau_born+log(tau_lower_bound)) )then
yij_low=-(tau_lower_bound+tau_born)* &
  & ( 1-em2ycm_born*tau_lower_bound ) / &
  & ( (tau_lower_bound-tau_born)* &
  & (1+em2ycm_born*tau_lower_bound) )
else
yij_low=-1.d0
endif
!
if(idir.eq.1)then
y_ij_fks_upp=yij_upp
y_ij_fks_low=yij_low
elseif(idir.eq.-1)then
y_ij_fks_upp=-yij_low
y_ij_fks_low=-yij_upp
endif

!
! set-up y_ij_fks
!
if(colltest)then
cctiny=0.d0
else
cctiny=ctiny
endif
if( (icountevts.eq.-100.or.icountevts.eq.0) .and. &
  & ((.not.softtest) .or. &
  & (softtest.and.y_ij_fks_fix.eq.-2.d0)) .and. &
  & (.not.colltest)  )then
! importance sampling towards collinear singularity
! insert here further importance sampling towards y_ij_fks->1
y_ij_fks = y_ij_fks_upp - &
  & (y_ij_fks_upp-y_ij_fks_low)*(cctiny+(1-cctiny)*x(2)**2)
elseif( (icountevts.eq.-100.or.icountevts.eq.0) .and. &
  & ((softtest.and.y_ij_fks_fix.ne.-2.d0) .or. &
  & colltest)  )then
y_ij_fks = y_ij_fks_fix
if ( y_ij_fks.gt.y_ij_fks_upp+1d-12 .or. &
  & y_ij_fks.lt.y_ij_fks_low-1d-12) then
xjac=-33d0
pass=.false.
return
endif
elseif(abs(icountevts).eq.2.or.abs(icountevts).eq.1)then
y_ij_fks=y_ij_fks_matrix(icountevts)
! Check that y_ij_fks is in the allowed range. If not, counter events
! cannot be generated
if ( y_ij_fks.gt.y_ij_fks_upp+1d-12 .or. &
  & y_ij_fks.lt.y_ij_fks_low-1d-12) then
xjac=-33d0
pass=.false.
return
endif
else
write(*,*)'Error #8 in genps_fks.f',icountevts
stop
endif
! importance sampling towards collinear singularity
xjac=xjac*(y_ij_fks_upp-y_ij_fks_low)*x(2)*2d0
!
! Compute costh_i_fks
!
yijdir=idir*y_ij_fks
costh_i_fks=yijdir
!
! Compute maximal xi_i_fks
!
x1bar2 = xbjrk_born(1)**2
omx1bar2 = 1d0-x1bar2
x2bar2 = xbjrk_born(2)**2
omx2bar2 = 1d0-x2bar2

if(1-tau_born.gt.1.d-5)then
yij_sol=-sinh(ycm_born)*(1+tau_born)/ &
  & ( cosh(ycm_born)*(1-tau_born) )
else
yij_sol=-ycmhat
endif
if(abs(yij_sol).gt.1.d0)then
if (abs(yij_sol).lt.1d0+qtiny) then
yij_sol = sign(1d0, yij_sol)
else
write(*,*)'Error #9 in genps_fks.f',yij_sol,icountevts
write(*,*)xbjrk_born(1),xbjrk_born(2),yijdir
endif
endif

if(yijdir.eq.yij_sol)then
ximaxtmp=1-xbjrk_born(1)*xbjrk_born(2)
elseif(yijdir.ge.yij_sol)then
!this is an expansion when both yij->-1 and x1->1
! in this case there may be precision loosses
! from the argument in the sqrt
if (abs(yijdir+1d0).lt.ctiny.and.omx1bar2.lt.ctiny) then
xi1=(4*x1bar2 + yijdir + 11*x1bar2*yijdir - 5*x1bar2**2*yijdir + &
  & x1bar2**3*yijdir+4*yijdir**2)/(2*(1 + yijdir)**2)
ximaxtmp=1-xi1
else if (omx1bar2.lt.ctiny) then
! compute directly ximaxtmp
ximaxtmp = omx1bar2 / (1+yijdir)
else
xi1=2*(1+yijdir)*x1bar2/( &
  & sqrt( ((1+x1bar2)*(1-yijdir))**2+16*yijdir*x1bar2 ) + &
  & (1-yijdir)*(omx1bar2) )
ximaxtmp=1-xi1
endif
elseif(yijdir.lt.yij_sol)then
!this is an expansion when both yij->+1 and x1->1
! in this case there may be precision loosses
! from the argument in the sqrt
if (abs(yijdir-1d0).lt.ctiny.and.omx2bar2.lt.ctiny) then
xi2=(4*x2bar2 - yijdir - 11*x2bar2*yijdir + 5*x2bar2**2*yijdir - &
  & x2bar2**3*yijdir +4*yijdir**2)/(4*(-1 + yijdir)**2)
ximaxtmp=1-xi2
else if (omx2bar2.lt.ctiny) then
! compute directly ximaxtmp
ximaxtmp = omx2bar2 / (1-yijdir)
else
xi2=2*(1-yijdir)*x2bar2/( &
  & sqrt( ((1+x2bar2)*(1+yijdir))**2-16*yijdir*x2bar2 ) + &
  & (1+yijdir)*(omx2bar2) )
ximaxtmp=1-xi2
endif
else
write(*,*)'Fatal error #14 in one_tree: unknown option'
write(*,*)y_ij_fks,yij_sol,idir
stop
endif
xiimax=ximaxtmp
!
! Lower bound on xi_i_fks
!
if (tau_born.lt.tau_lower_bound) then
xiimin=1d0-tau_born/tau_lower_bound
else
xiimin=0d0
endif
if (xiimax.lt.xiimin) then
write (*,*) 'WARNING #10 in genps_fks.f',icountevts,xiimax &
  & ,xiimin
xjac=-342d0
pass=.false.
return
endif

xinorm=xiimax-xiimin
if( icountevts.ge.1 .and. &
  & ( (idir.eq.1.and. &
  & abs(ximaxtmp-(1-xbjrk_born(1))).gt.1.d-5) .or. &
  & (idir.eq.-1.and. &
  & abs(ximaxtmp-(1-xbjrk_born(2))).gt.1.d-5) ) )then
write(*,*)'Fatal error #15 in one_tree'
write(*,*)ximaxtmp,xbjrk_born(1),xbjrk_born(2),idir
stop
endif
!
! Define xi_i_fks
!
if( (icountevts.eq.-100.or.abs(icountevts).eq.1) .and. &
  & ((.not.colltest) .or. &
  & (colltest.and.xi_i_fks_fix.eq.-2.d0)) .and. &
  & (.not.softtest)  )then
if(icountevts.eq.-100)then
! importance sampling towards soft singularity
! insert here further importance sampling towards xi_i_hat->0
xi_i_hat=sstiny+(1-sstiny)*x(1)**2
endif
xi_i_fks=xiimin+(xiimax-xiimin)*xi_i_hat
elseif( (icountevts.eq.-100.or.abs(icountevts).eq.1) .and. &
  & (colltest.and.xi_i_fks_fix.ne.-2.d0) .and. &
  & (.not.softtest)  )then
if(xi_i_fks_fix.lt.xiimax)then
xi_i_fks=xi_i_fks_fix
else
xi_i_fks=xi_i_fks_fix*xiimax
endif
elseif( (icountevts.eq.-100.or.abs(icountevts).eq.1) .and. &
  & softtest )then
if(xi_i_fks_fix.lt.xiimax)then
xi_i_fks=xi_i_fks_fix
else
xjac=-102
pass=.false.
return
endif
elseif(abs(icountevts).eq.2.or.icountevts.eq.0)then
xi_i_fks=xi_i_fks_matrix(icountevts)
! Check that xi_i_fks is in the allowed range. If not, counter events
! cannot be generated
if ( xi_i_fks.gt.xiimax+1d-12 .or. &
  & xi_i_fks.lt.xiimin-1d-12 ) then
xjac=-34d0
pass=.false.
return
endif
else
write(*,*)'Error #11 in genps_fks.f',icountevts
stop
endif
! remove the following if no importance sampling towards soft
! singularity is performed when integrating over xi_i_hat
xjac=xjac*2d0*x(1)
!
! Initial state variables are different for events and counterevents. Update them here.
!
omega=sqrt( (2-xi_i_fks*(1+yijdir))/ &
  & (2-xi_i_fks*(1-yijdir)) )
if (icountevts.ne.0) then
tau=tau_born/(1-xi_i_fks)
ycm=ycm_born-log(omega)
shat=tau*stot
sqrtshat=sqrt(shat)
xbjrk(1)=xbjrk_born(1)/(sqrt(1-xi_i_fks)*omega)
xbjrk(2)=xbjrk_born(2)*omega/sqrt(1-xi_i_fks)
else
tau=tau_born
ycm=ycm_born
shat=shat_born
sqrtshat=sqrt(shat)
xbjrk(1)=xbjrk_born(1)
xbjrk(2)=xbjrk_born(2)
endif
!
! Define the boost factor here
!
bstfact=sqrt( (2-xi_i_fks*(1-yijdir))*(2-xi_i_fks*(1+yijdir)) )
shy_tbst=-xi_i_fks*sqrt(1-yijdir**2)/(2*sqrt(1-xi_i_fks))
chy_tbst=bstfact/(2*sqrt(1-xi_i_fks))
chy_tbstmo=chy_tbst-1.d0
cosphi_i_fks=cos(phi_i_fks)
sinphi_i_fks=sin(phi_i_fks)
xdir_t(1)=-cosphi_i_fks
xdir_t(2)=-sinphi_i_fks
xdir_t(3)=zero
!
shy_lbst=-xi_i_fks*yijdir/bstfact
chy_lbst=(2-xi_i_fks)/bstfact
! Boost the momenta
do i=3,nexternal
if(i.ne.i_fks.and.shy_tbst.ne.0.d0) &
  & call boostwdir2_in_place(chy_tbst,shy_tbst,chy_tbstmo, &
  & xdir_t,xp(0,i))
enddo
!
encmso2=sqrtshat/2.d0
p_i_fks(0)=encmso2
E_i_fks=xi_i_fks*encmso2
sinth_i_fks=sqrt(1-costh_i_fks**2)
!
xp(0,1)=encmso2*(chy_lbst-shy_lbst)
xp(1,1)=0.d0
xp(2,1)=0.d0
xp(3,1)=xp(0,1)
!
xp(0,2)=encmso2*(chy_lbst+shy_lbst)
xp(1,2)=0.d0
xp(2,2)=0.d0
xp(3,2)=-xp(0,2)
!
xp(0,i_fks)=E_i_fks*(chy_lbst-shy_lbst*yijdir)
p_i_fks(0)=p_i_fks(0)*(chy_lbst-shy_lbst*yijdir)
xpifksred(1)=sinth_i_fks*cosphi_i_fks
xpifksred(2)=sinth_i_fks*sinphi_i_fks
xpifksred(3)=chy_lbst*yijdir-shy_lbst
!
do j=1,3
xp(j,i_fks)=E_i_fks*xpifksred(j)
p_i_fks(j)=encmso2*xpifksred(j)
enddo
!
! Collinear limit of <ij>/[ij]. See innerpin.m.
if( icountevts.eq.-100 .or. &
  & (icountevts.eq.1.and.xij_aor.eq.0) )then
resAoR0=-exp( 2*idir*ximag*phi_i_fks )
xij_aor=resAoR0
endif
!
! Phase-space factor for (xii,yij,phii) * (tau,ycm)
xpswgt=xpswgt*shat
xpswgt=xpswgt/(4*pi)**3/(1-xi_i_fks)
xpswgt=abs(xpswgt)
!
return
end subroutine generate_momenta_initial


subroutine getangles(pin,th,cth,sth,phi,cphi,sphi)
implicit none
double precision pin(0:3),th,cth,sth,phi,cphi,sphi,xlength
!
xlength=pin(1)**2+pin(2)**2+pin(3)**2
if(xlength.eq.0)then
th=0.d0
cth=1.d0
sth=0.d0
phi=0.d0
cphi=1.d0
sphi=0.d0
else
xlength=sqrt(xlength)
cth=pin(3)/xlength
th=acos(cth)
if(cth.ne.1.d0)then
sth=sqrt(1-cth**2)
phi=atan2(pin(2),pin(1))
cphi=cos(phi)
sphi=sin(phi)
else
sth=0.d0
phi=0.d0
cphi=1.d0
sphi=0.d0
endif
endif
return
end subroutine getangles

subroutine gentcms(pa,pb,t,phi,m1,m2,p1,pr,jac)
!*************************************************************************
!     Generates 4 momentum for particle 1, and remainder pr
!     given the values t, and phi
!     Assuming incoming particles with momenta pa, pb
!     And outgoing particles with mass m1,m2
!     s = (pa+pb)^2  t=(pa-p1)^2
!*************************************************************************
implicit none
!
!     Arguments
!
double precision t,phi,m1,m2               !inputs
double precision pa(0:3),pb(0:3),jac
double precision p1(0:3),pr(0:3)           !outputs
!
!     local
!
double precision ptot(0:3),E_acms,p_acms,pa_cms(0:3)
double precision esum,ed,pp,md2,ma2,pt,ptotm(0:3)
integer i
!
!     External
!
!-----
!  Begin Code
!-----
do i=0,3
ptot(i)  = pa(i)+pb(i)
if (i .gt. 0) then
ptotm(i) = -ptot(i)
else
ptotm(i) = ptot(i)
endif
enddo
ma2 = dot(pa,pa)
!
!     determine magnitude of p1 in cms frame (from dhelas routine mom2cx)
!
ESUM = sqrt(max(0d0,dot(ptot,ptot)))
if (esum .eq. 0d0) then
jac=-8d0             !Failed esum must be > 0
return
endif
MD2=(M1-M2)*(M1+M2)
ED=MD2/ESUM
IF (M1*M2.EQ.0.) THEN
PP=(ESUM-ABS(ED))*0.5d0
ELSE
PP=(MD2/ESUM)**2-2.0d0*(M1**2+M2**2)+ESUM**2
if (pp .gt. 0) then
PP=SQRT(pp)*0.5d0
else
write(*,*) 'Warning #12 in genps_fks.f',pp
jac=-1
return
endif
ENDIF
!
!     Energy of pa in pa+pb cms system
!
call boostx(pa,ptotm,pa_cms)
E_acms = pa_cms(0)
p_acms = dsqrt(pa_cms(1)**2+pa_cms(2)**2+pa_cms(3)**2)
!
p1(0) = MAX((ESUM+ED)*0.5d0,0.d0)
p1(3) = -(m1*m1+ma2-t-2d0*p1(0)*E_acms)/(2d0*p_acms)
pt = dsqrt(max(pp*pp-p1(3)*p1(3),0d0))
p1(1) = pt*cos(phi)
p1(2) = pt*sin(phi)
!
call rotxxx(p1,pa_cms,p1)          !Rotate back to pa_cms frame
call boostx(p1,ptot,p1)            !boost back to lab fram
do i=0,3
pr(i)=pa(i)-p1(i)               !Return remainder of momentum
enddo
end subroutine gentcms


DOUBLE PRECISION FUNCTION LAMBDA(S,MA2,MB2)
IMPLICIT NONE
!****************************************************************************
!     THIS IS THE LAMBDA FUNCTION FROM VERNONS BOOK COLLIDER PHYSICS P 662
!     MA2 AND MB2 ARE THE MASS SQUARED OF THE FINAL STATE PARTICLES
!     2-D PHASE SPACE = .5*PI*SQRT(1.,MA2/S^2,MB2/S^2)*(D(OMEGA)/4PI)
!****************************************************************************
DOUBLE PRECISION MA2,MB2,S,tiny,tmp,rat
parameter (tiny=1.d-8)
!
tmp=S**2+MA2**2+MB2**2-2d0*S*MA2-2d0*MA2*MB2-2d0*S*MB2
if(tmp.le.0.d0)then
if(ma2.lt.0.d0.or.mb2.lt.0.d0)then
write(6,*)'Error #1 in function Lambda:',s,ma2,mb2
stop
endif
rat=1-(sqrt(ma2)+sqrt(mb2))/s
if(rat.gt.-tiny)then
tmp=0.d0
else
write(6,*)'Error #2 in function Lambda:',s,ma2,mb2,rat
endif
endif
LAMBDA=tmp
RETURN
end function LAMBDA


SUBROUTINE YMINMAX(X,Y,Z,U,V,W,YMIN,YMAX)
!**************************************************************************
!     This is the G function from Particle Kinematics by
!     E. Byckling and K. Kajantie, Chapter 4 p. 91 eqs 5.28
!     It is used to determine physical limits for Y based on inputs
!**************************************************************************
implicit none
!
!     Constant
!
double precision tiny
parameter       (tiny=1d-199)
!
!     Arguments
!
Double precision x,y,z,u,v,w              !inputs  y is dummy
Double precision ymin,ymax                !output
!
!     Local
!
double precision y1,y2,yr,ysqr
!
!-----
!  Begin Code
!-----
ysqr = lambda(x,u,v)*lambda(x,w,z)
if (ysqr .ge. 0d0) then
yr = dsqrt(ysqr)
else
print*,'Error in yminymax sqrt(-x)',lambda(x,u,v),lambda(x,w,z)
yr=0d0
endif
y1 = u+w -.5d0* ((x+u-v)*(x+w-z) - yr)/(x+tiny)
y2 = u+w -.5d0* ((x+u-v)*(x+w-z) + yr)/(x+tiny)
ymin = min(y1,y2)
ymax = max(y1,y2)
end subroutine YMINMAX


subroutine compute_tau_one_body(totmass,stot,tau,jac)
implicit none
double precision totmass,stot,tau,jac,roH
roH=totmass**2/stot
tau=roH
! Jacobian due to delta() of tau_born
jac=jac*2*totmass/stot
return
end subroutine compute_tau_one_body


subroutine generate_tau_BW(stot,idim,x,mass,width,cBW,BWmass &
  & ,BWwidth,tau,jac)
implicit none
integer cBW,idim
double precision stot,x,tau,jac,mass,width,BWmass(-1:1),BWwidth( &
  & -1:1),s_mass,s
double precision smax,smin
double precision tau_Born_lower_bound,tau_lower_bound_resonance &
  & ,tau_lower_bound
common/ctau_lower_bound/tau_Born_lower_bound &
  & ,tau_lower_bound_resonance,tau_lower_bound
if (cBW.eq.1 .and. width.gt.0d0 .and. BWwidth(1).gt.0d0) then
smin=tau_Born_lower_bound*stot
smax=stot
s_mass=smin
call trans_x(5,idim,x,smin,smax,s_mass,mass,width,BWmass( &
  & -1),BWwidth(-1),jac,s)
tau=s/stot
jac=jac/stot
else
smin=tau_Born_lower_bound*stot
smax=stot
s_mass=smin
call trans_x(3,idim,x,smin,smax,s_mass,mass,width,BWmass( &
  & -1),BWwidth(-1),jac,s)
tau=s/stot
jac=jac/stot
endif
return
end subroutine generate_tau_BW


subroutine generate_tau(stot,idim,x,tau,jac)
implicit none
integer idim
double precision x,tau,jac,smin,smax,s_mass,s,tiny,dum,dum3(-1:1) &
  & ,stot
parameter (tiny=1d-8)
double precision tau_Born_lower_bound,tau_lower_bound_resonance &
  & ,tau_lower_bound
common/ctau_lower_bound/tau_Born_lower_bound &
  & ,tau_lower_bound_resonance,tau_lower_bound
smin=tau_born_lower_bound*stot
smax=stot
s_mass=tau_lower_bound_resonance*stot
if (s_mass.gt.smin*(1d0+tiny)) then
call trans_x(2,idim,x,smin,smax,s_mass,dum,dum &
  & ,dum3,dum3,jac,s)
elseif(abs(s_mass-smin).lt.tiny*smin) then
call trans_x(7,idim,x,smin,smax,s_mass,dum,dum &
  & ,dum3,dum3,jac,s)
else
write (*,*) 'ERROR #39 in genps_fks.f',s_mass,smin,smax
jac=-1d0
endif
tau=s/stot
jac=jac/stot
return
end subroutine generate_tau


subroutine generate_y(tau,x,ycm,ycmhat,jac)
implicit none
double precision tau,x,ycm,jac
double precision ylim,ycmhat
ylim=-0.5d0*log(tau)
ycmhat=2*x-1
ycm=ylim*ycmhat
jac=jac*ylim*2
return
end subroutine generate_y


subroutine compute_tau_y_epem(j_fks,one_body,fksmass, &
  & stot,tau,ycm,ycmhat)
implicit none
integer j_fks
logical one_body
double precision fksmass,stot,tau,ycm,ycmhat
if(j_fks.le.nincoming)then
! This should never happen in normal integration: when no PDFs, j_fks
! cannot be initial state (but needed for testing). If tau set to one,
! integration range in xi_i_fks will be zero, so lower it artificially
! when too large
if(one_body)then
tau=fksmass**2/stot
else
tau=max((0.85d0)**2,fksmass**2/stot)
endif
ycm=0.d0
else
! For e+e- collisions, set tau to one and y to zero
tau=1.d0
ycm=0.d0
endif
ycmhat=0.d0
return
end subroutine compute_tau_y_epem


subroutine generate_inv_mass_sch(ns_channel,itree,m,sqrtshat_born &
  & ,totmass,qwidth,qmass,cBW,cBW_mass,cBW_width,s,x,xjac0,pass)
implicit none
integer ns_channel
double precision qmass(-nexternal:0),qwidth(-nexternal:0)
double precision M(-max_branch:max_particles),x(99)
double precision s(-max_branch:max_particles)
double precision sqrtshat_born,totmass,xjac0
integer itree(2,-max_branch:-1)
integer i,j,ii,order(-nexternal:0)
double precision smin,smax,totalmass
logical pass
integer cBW(-nexternal:-1)
double precision cBW_mass(-1:1,-nexternal:-1),cBW_width(-1:1, &
  & -nexternal:-1)
pass=.true.
totalmass=totmass
do ii = -1,-ns_channel,-1
! Randomize the order with which to generate the s-channel masses:
call sChan_order(ns_channel,order)
i=order(ii)
!     Generate invariant masses for all s-channel branchings of the Born
smin = (m(itree(1,i))+m(itree(2,i)))**2
smax = (sqrtshat_born-totalmass+sqrt(smin))**2
if(smax.lt.smin.or.smax.lt.0.d0.or.smin.lt.0.d0)then
write(*,*)'Error #13 in genps_fks.f'
write(*,*)smin,smax,i
stop
endif
call generate_si(i,smin,smax,s,cBW,cBW_width,cBW_mass,qmass &
  & ,qwidth,x,xjac0,schannel_masses)
! If numerical inaccuracy, quit loop
if (xjac0 .lt. 0d0) then
if ((xjac0.gt.-400d0 .or. xjac0.le.-500d0) .and. &
  & xjac0.ne.0d0)then
write (*,*) 'WARNING #31 in genps_fks.f',i,s(i),smin,smax &
  & ,xjac0
endif
xjac0 = -6
pass=.false.
return
endif
if (s(i) .lt. smin) then
write (*,*) 'WARNING #32 in genps_fks.f',i,s(i),smin,smax,x( &
  & -i)
xjac0=-5
pass=.false.
return
endif
!
!     fill masses, update totalmass
!
m(i) = sqrt(s(i))
totalmass=totalmass+m(i)- &
  & m(itree(1,i))-m(itree(2,i))
if ( totalmass.gt.sqrtshat_born )then
write (*,*) 'WARNING #33 in genps_fks.f',i,totalmass &
  & ,sqrtshat_born,s(i)
xjac0 = -4
pass=.false.
return
endif
enddo
return
end subroutine generate_inv_mass_sch


subroutine generate_si(i,smin,smax,s,cBW,cBW_width,cBW_mass,qmass &
  & ,qwidth,x,xjac0,s_mass)
implicit none
integer i
double precision smin,smax,s(-max_branch:max_particles),qwidth( &
  & -nexternal:0),qmass(-nexternal:0),cBW_width(-1:1,-nexternal: &
  & -1),cBW_mass(-1:1,-nexternal:-1),xjac0,x(99),s_mass( &
  & -nexternal:nexternal)
integer cBW(-nexternal:-1)
! Choose the appropriate s given our constraints smin,smax
if(qwidth(i).ne.0.d0 .and. cBW(i).ne.2)then
! Breit Wigner
if (cBW(i).eq.1 .and. &
  & cBW_width(1,i).gt.0d0 .and. cBW_width(-1,i).gt.0d0) then
!     conflicting BW on both sides
call trans_x(6,-i,x(-i),smin,smax,s_mass(i),qmass(i) &
  & ,qwidth(i),cBW_mass(-1,i),cBW_width(-1,i),xjac0,s(i))
elseif (cBW(i).eq.1.and.cBW_width(1,i).gt.0d0) then
!     conflicting BW with alternative mass larger
call trans_x(5,-i,x(-i),smin,smax,s_mass(i),qmass(i) &
  & ,qwidth(i),cBW_mass(-1,i),cBW_width(-1,i),xjac0,s(i))
elseif (cBW(i).eq.1.and.cBW_width(-1,i).gt.0d0) then
!     conflicting BW with alternative mass smaller
call trans_x(4,-i,x(-i),smin,smax,s_mass(i),qmass(i) &
  & ,qwidth(i),cBW_mass(-1,i),cBW_width(-1,i),xjac0,s(i))
else
!     normal BW
call trans_x(3,-i,x(-i),smin,smax,s_mass(i),qmass(i) &
  & ,qwidth(i),cBW_mass(-1,i),cBW_width(-1,i),xjac0,s(i))
endif
else
! not a Breit Wigner
if (smin.eq.0d0 .and. s_mass(i).eq.0d0) then
!     no lower limit on invariant mass from cuts or final state masses:
!     use flat distribution
call trans_x(1,-i,x(-i),smin,smax,s_mass(i),qmass(i) &
  & ,qwidth(i),cBW_mass(-1,i),cBW_width(-1,i),xjac0,s(i))
elseif (smin.ge.s_mass(i) .and. smin.gt.0d0) then
!     A lower limit on smin, which is larger than lower limit from cuts
!     or masses. Use 1/x importance sampling
call trans_x(7,-i,x(-i),smin,smax,s_mass(i),qmass(i) &
  & ,qwidth(i),cBW_mass(-1,i),cBW_width(-1,i),xjac0,s(i))
elseif (smin.lt.s_mass(i) .and. s_mass(i).gt.0d0) then
!     Use flat grid between smin and s_mass(i), and 1/x^nsamp above
!     s_mass(i)
call trans_x(2,-i,x(-i),smin,smax,s_mass(i),qmass(i) &
  & ,qwidth(i),cBW_mass(-1,i),cBW_width(-1,i),xjac0,s(i))
else
write (*,*) "ERROR in genps_fks.f:"// &
  & " cannot set s-channel without BW",i,smin,s_mass(i)
stop 1
endif
endif
return
end subroutine generate_si

subroutine generate_t_channel_branchings(ns_channel,nbranch,itree &
  & ,m,s,x,pb,xjac0,xpswgt0,pass)
! First we need to determine the energy of the remaining particles this
! is essentially in place of the cos(theta) degree of freedom we have
! with the s channel decay sequence
implicit none
double precision pi
parameter (pi=3.1415926535897932d0)
double precision xjac0,xpswgt0
double precision M(-max_branch:max_particles),x(99)
double precision s(-max_branch:max_particles)
double precision pb(0:3,-max_branch:nexternal-1)
integer itree(2,-max_branch:-1)
integer ns_channel,nbranch
logical pass
!
double precision totalmass,smin,smax,s1,ma2,mbq,m12,mnq,tmin,tmax &
  & ,t,tmax_temp,phi,dum,dum3(-1:1),s_m,tm,tiny
parameter (tiny=1d-8)
integer i,ibranch,idim
!
pass=.true.
totalmass=0d0
s_m=0d0
do ibranch = -ns_channel-1,-nbranch,-1
totalmass=totalmass+m(itree(2,ibranch))
s_m=s_m+sqrt(schannel_masses(itree(2,ibranch)))
enddo
m(-ns_channel-1) = dsqrt(S(-nbranch))
!
! Choose invariant masses of the pseudoparticles obtained by taking together
! all final-state particles or pseudoparticles found from the current
! t-channel propagator down to the initial-state particle found at the end
! of the t-channel line.
do ibranch = -ns_channel-1,-nbranch+2,-1
totalmass=totalmass-m(itree(2,ibranch))
smin = totalmass**2
smax = (m(ibranch) - m(itree(2,ibranch)))**2
if (smin .gt. smax) then
xjac0=-3d0
pass=.false.
return
endif
idim=(nbranch-1+(-ibranch)*2)
s_m=s_m-sqrt(schannel_masses(itree(2,ibranch)))
if (abs(smin-s_m**2).lt.tiny) then
call trans_x(1,idim,x(idim),smin,smax,s_m**2,dum &
  & ,dum,dum3(-1),dum3(-1),xjac0,s1)
else
call trans_x(1,idim,x(idim),smin,smax,s_m**2,dum &
  & ,dum,dum3(-1),dum3(-1),xjac0,s1)
endif
if (xjac0.le.0d0) then
if ((xjac0.gt.-400d0 .or. xjac0.le.-500d0) .and. &
  & xjac0.ne.0d0)then
write (*,*) 'WARNING #31a in genps_fks.f',ibranch,s1 &
  & ,smin,smax,s_m**2,xjac0
endif
xjac0 = -6
pass=.false.
return
endif
m(ibranch-1)=sqrt(s1)
if (m(ibranch-1)**2.lt.smin.or.m(ibranch-1)**2.gt.smax &
  & .or.m(ibranch-1).ne.m(ibranch-1)) then
xjac0=-1d0
pass=.false.
return
endif
enddo
!
! Set m(-nbranch) equal to the mass of the particle or pseudoparticle P
! attached to the vertex (P,t,p2), with t being the last t-channel propagator
! in the t-channel line, and p2 the incoming particle opposite to that from
! which the t-channel line starts
m(-nbranch) = m(itree(2,-nbranch))
!
!     Now perform the t-channel decay sequence. Most of this comes from:
!     Particle Kinematics Chapter 6 section 3 page 166
!
!     From here, on we can just pretend this is a 2->2 scattering with
!     Pa                    + Pb     -> P1          + P2
!     p(0,itree(ibranch,1)) + p(0,2) -> p(0,ibranch)+ p(0,itree(ibranch,2))
!     M(ibranch) is the total mass available (Pa+Pb)^2
!     M(ibranch-1) is the mass of P2  (all the remaining particles)
!
do ibranch=-ns_channel-1,-nbranch+1,-1
s1  = m(ibranch)**2    !Total mass available
ma2 = m(2)**2
mbq = dot(pb(0,itree(1,ibranch)),pb(0,itree(1,ibranch)))
m12 = m(itree(2,ibranch))**2
mnq = m(ibranch-1)**2
call yminmax(s1,t,m12,ma2,mbq,mnq,tmin,tmax)
call trans_x(1,-ibranch,x(-ibranch),-tmax,-tmin, &
  & schannel_masses(ibranch) &
  & ,dum,dum,dum3(-1),dum3(-1),xjac0,tm)
if (xjac0.le.0d0) then
if ((xjac0.gt.-400d0 .or. xjac0.le.-500d0) .and. &
  & xjac0.ne.0d0)then
write (*,*) 'WARNING #31b in genps_fks.f',ibranch,tm &
  & ,-tmax,-tmin,xjac0
endif
xjac0 = -6
pass=.false.
return
endif
t=-tm
if (t .lt. tmin .or. t .gt. tmax) then
write (*,*) "WARNING #35 in genps_fks.f",t,tmin,tmax
xjac0=-3d0
pass=.false.
return
endif
phi = 2d0*pi*x(nbranch+(-ibranch-1)*2)
xjac0 = xjac0*2d0*pi
! Finally generate the momentum. The call is of the form
! pa+pb -> p1+ p2; t=(pa-p1)**2;   pr = pa-p1
! gentcms(pa,pb,t,phi,m1,m2,p1,pr)
call gentcms(pb(0,itree(1,ibranch)),pb(0,2),t,phi, &
  & m(itree(2,ibranch)),m(ibranch-1),pb(0,itree(2,ibranch)), &
  & pb(0,ibranch),xjac0)
!
if (xjac0 .lt. 0d0) then
write(*,*) 'Failed gentcms',ibranch,xjac0
pass=.false.
return
endif
xpswgt0 = xpswgt0/(4d0*dsqrt(lambda(s1,ma2,mbq)))
enddo
! We need to get the momentum of the last external particle.  This
! should just be the sum of p(0,2) and the remaining momentum from our
! last t channel 2->2
do i=0,3
pb(i,itree(2,-nbranch)) = pb(i,-nbranch+1)+pb(i,2)
enddo
return
end subroutine generate_t_channel_branchings


subroutine fill_born_momenta(nbranch,nt_channel,one_body,ionebody &
  & ,x,itree,m,s,pb,xjac0,xpswgt0,pass)
implicit none
double precision pi
parameter (pi=3.1415926535897932d0)
integer nbranch,nt_channel,ionebody
double precision M(-max_branch:max_particles),x(99)
double precision s(-max_branch:max_particles)
double precision pb(0:3,-max_branch:nexternal-1)
integer itree(2,-max_branch:-1)
double precision xjac0,xpswgt0
logical pass,one_body
!
double precision one
parameter (one=1d0)
double precision costh,phi,xa2,xb2
integer i,ix
double precision vtiny
parameter (vtiny=1d-12)
!
pass=.true.
do i = -nbranch+nt_channel+(nincoming-1),-1
ix = nbranch+(-i-1)*2+(2-nincoming)
if (nt_channel .eq. 0) ix=ix-1
costh= 2d0*x(ix)-1d0
phi  = 2d0*pi*x(ix+1)
xjac0 = xjac0 * 4d0*pi
xa2 = m(itree(1,i))*m(itree(1,i))/s(i)
xb2 = m(itree(2,i))*m(itree(2,i))/s(i)
if (m(itree(1,i))+m(itree(2,i)) .ge. m(i)) then
xjac0=-8
pass=.false.
return
endif
xpswgt0 = xpswgt0*.5D0*PI*SQRT(LAMBDA(ONE,XA2,XB2))/(4.D0*PI)
call mom2cx(m(i),m(itree(1,i)),m(itree(2,i)),costh,phi, &
  & pb(0,itree(1,i)),pb(0,itree(2,i)))
! If there is an extremely large boost needed here, skip the phase-space point
! because of numerical stabilities.
if (dsqrt(abs(dot(pb(0,i),pb(0,i))))/pb(0,i) &
  & .lt.vtiny) then
xjac0=-81
pass=.false.
return
else
call boostm(pb(0,itree(1,i)),pb(0,i),m(i),pb(0,itree(1,i)))
call boostm(pb(0,itree(2,i)),pb(0,i),m(i),pb(0,itree(2,i)))
endif
enddo
!
!
! Special phase-space fix for the one_body
if (one_body) then
! Factor due to the delta function in dphi_1
xpswgt0=pi/m(ionebody)
! Kajantie's normalization of phase space (compensated below in flux)
xpswgt0=xpswgt0/(2*pi)
do i=0,3
pb(i,3) = pb(i,1)+pb(i,2)
enddo
endif
return
end subroutine fill_born_momenta


subroutine get_recoil(p_born,imother,shat_born,xmrec2,pass)
implicit none
double precision p_born(0:3,nexternal-1),xmrec2,shat_born
logical pass
integer imother,i
double precision recoilbar(0:3)
pass=.true.
do i=0,3
if (nincoming.eq.2) then
recoilbar(i)=p_born(i,1)+p_born(i,2)-p_born(i,imother)
else
recoilbar(i)=p_born(i,1)-p_born(i,imother)
endif
enddo
xmrec2=dot(recoilbar,recoilbar)
if(xmrec2.lt.0.d0)then
if(abs(xmrec2).gt.(1.d-4*shat_born))then
write(*,*)'Fatal error #14 in genps_fks.f',xmrec2,imother
stop
else
write(*,*)'Error #15 in genps_fks.f',xmrec2,imother
pass=.false.
return
endif
endif
if (xmrec2.ne.xmrec2) then
write (*,*) 'Error #16 in setting up event in genps_fks.f,'// &
  & ' skipping event'
pass=.false.
return
endif
return
end subroutine get_recoil



subroutine trans_x(itype,idim,x,smin,smax,s_mass,qmass,qwidth &
  & ,cBW_mass,cBW_width,jac,s)
! Given the input random number 'x', returns the corresponding value of
! the invariant mass squared 's'.
!
!     itype=1: flat transformation
!     itype=2: flat between 0 and s_mass/stot, 1/x above
!     itype=3: Breit-Wigner
!     itype=4: Conflicting BW, with alternative mass smaller
!     itype=5: Conflicting BW, with alternative mass larger
!     itype=6: Conflicting BW on both sides
!
implicit none
integer itype,idim
double precision x,smin,smax,s_mass,qmass,qwidth,cBW_mass(-1:1) &
  & ,cBW_width(-1:1),jac,s
double precision fract,A,B,C,bs(-1:1),maxi,mini
integer j
!
if (itype.eq.1) then
!     flat transformation:
A=smax-smin
B=smin
s=A*x+B
jac=jac*A
elseif (itype.eq.2) then
fract=0.25d0
if (s_mass.eq.0d0) then
write (*,*) 's_mass is zero',itype,idim
endif
if (x.lt.fract) then
!     flat transformation:
if (s_mass.lt.smin) then
jac=-421d0
return
endif
maxi=min(s_mass,smax)
A=(maxi-smin)/fract
B=smin
s=A*x+B
jac=jac*A
else
!     S=A/(B-x) transformation:
if (s_mass.ge.smax) then
jac=-422d0
return
endif
mini=max(s_mass,smin)
A=mini*smax*(1d0-fract)/(smax-mini)
B=(smax-fract*mini)/(smax-mini)
s=A/(B-x)
jac=jac*s**2/A
endif
elseif(itype.eq.3) then
!     Normal Breit-Wigner, i.e.
!        \int_smin^smax ds g(s)/((s-qmass^2)^2-qmass^2*qwidth^2) =
!        \int_0^1 dx g(s(x))
A=atan((qmass-smin/qmass)/qwidth)
B=atan((qmass-smax/qmass)/qwidth)
s=qmass*(qmass-qwidth*tan(A-(A-B)*x))
jac=jac*qmass*qwidth*(A-B)/(cos(A-(A-B)*x))**2
elseif(itype.eq.4) then
!     Conflicting BW, with alternative mass smaller than current
!     mass. That is, we need to throw also many events at smaller masses
!     than the peak of the current BW. Split 'x' at 'bs(-1)', using a
!     flat distribution below the split, and a BW above the split.
fract=0.3d0
bs(-1)=(cBW_mass(-1)-qmass)/ &
  & (qwidth+cBW_width(-1)) ! bs(-1) is negative here
bs(-1)=qmass+bs(-1)*qwidth
bs(-1)=bs(-1)**2
if (x.lt.fract) then
if(smin.gt.bs(-1)) then
jac=-441d0
return
endif
maxi=min(bs(-1),smax)
A=(maxi-smin)/fract
B=smin
s=A*x+B
jac=jac*A
else
if(smax.lt.bs(-1)) then
jac=-442d0
return
endif
mini=max(bs(-1),smin)
A=atan((qmass-mini/qmass)/qwidth)
B=atan((qmass-smax/qmass)/qwidth)
C=((1d0-x)*A+(x-fract)*B)/(1d0-fract)
s=qmass*(qmass-qwidth*tan(C))
jac=jac*qmass*qwidth*(A-B)/((cos(C))**2*(1d0-fract))
endif
elseif(itype.eq.5) then
!     Conflicting BW, with alternative mass larger than current
!     mass. That is, we need to throw also many events at larger masses
!     than the peak of the current BW. Split 'x' at 'bs(1)' and the
!     alternative mass. Use a BW below bs(1), a flat distribution
!     between bs(1) and the alternative mass, and a 1/x above the
!     alternative mass.
fract=0.35d0
bs(1)=(cBW_mass(1)-qmass)/ &
  & (qwidth+cBW_width(1))
bs(1)=qmass+bs(1)*qwidth
bs(1)=bs(1)**2
if (x.lt.fract) then
if(smin.gt.bs(1)) then
jac=-451d0
return
endif
maxi=min(bs(1),smax)
A=atan((qmass-smin/qmass)/qwidth)
B=atan((qmass-maxi/qmass)/qwidth)
C=((B-A)*x+fract*A)/fract
s=qmass*(qmass-qwidth*tan(C))
jac=jac*qmass*qwidth*(A-B)/((cos(C))**2*fract)
elseif (x.lt.1d0-fract) then
if(smin.gt.cBW_mass(1)**2 .or. smax.lt.bs(1)) then
jac=-452d0
return
endif
maxi=min(cBW_mass(1)**2,smax)
mini=max(bs(1),smin)
A=(maxi-mini)/(1d0-2d0*fract)
B=((1d0-fract)*mini-fract*maxi)/(1d0-2d0*fract)
s=A*x+B
jac=jac*A
else
if(smax.le.cBW_mass(1)**2) then
jac=-453d0
return
endif
mini=max(cBW_mass(1)**2,smin)
A=mini*smax*fract/(smax-mini)
B=(smax-(1d0-fract)*mini)/(smax-mini)
s=A/(B-x)
jac=jac*s**2/A
endif
elseif(itype.eq.6) then
fract=0.3d0
!     Conflicting BW on both sides. Use flat below bs(-1); BW between
!     bs(-1) and bs(1); flat between bs(1) and alternative mass; and 1/x
!     above alternative mass.
do j=-1,1,2
bs(j)=(cBW_mass(j)-qmass)/ &
  & (qwidth+cBW_width(j))
bs(j)=qmass+bs(j)*qwidth
bs(j)=bs(j)**2
enddo
if (x.lt.fract) then
if(smin.gt.bs(-1)) then
jac=-461d0
return
endif
maxi=min(bs(-1),smax)
A=(maxi-smin)/fract
B=smin
s=A*x+B
jac=jac*A
elseif(x.lt.1d0-fract) then
if(smin.gt.bs(1) .or. smax.lt.bs(-1)) then
jac=-462d0
return
endif
maxi=min(bs(1),smax)
mini=max(bs(-1),smin)
A=atan((qmass-mini/qmass)/qwidth)
B=atan((qmass-maxi/qmass)/qwidth)
C=((1d0-fract-x)*A+(x-fract)*B)/(1d0-2d0*fract)
s=qmass*(qmass-qwidth*tan(C))
jac=-jac*qmass*qwidth*(B-A)/((cos(C))**2*(1d0-2d0*fract))
elseif(x.lt.1d0-fract/2d0) then
if(smin.gt.cBW_mass(1)**2 .or. smax.lt.bs(1)) then
jac=-463d0
return
endif
maxi=min(cBW_mass(1)**2,smax)
mini=max(bs(1),smin)
A=2d0*(maxi-mini)/fract
B=2d0*maxi-mini-2d0*(maxi-mini)/fract
s=A*x+B
jac=jac*A
else
if(smax.le.cBW_mass(1)**2) then
jac=-464d0
return
endif
mini=max(cBW_mass(1)**2,smin)
A=mini*smax*fract/(2d0*(smax-mini))
B=(smax-(1d0-fract/2d0)*mini)/(smax-mini)
s=A/(B-x)
jac=jac*s**2/A
endif
elseif (itype.eq.7) then
!     S=A/(B-x) transformation:
if (smin.le.0d0) then
jac=-471d0
return
endif
A=smin*smax/(smax-smin)
B=smax/(smax-smin)
s=A/(B-x)
jac=jac*s**2/A
endif
return
end subroutine trans_x



subroutine validate_bound_genps_state()
implicit none

if (size(config_mass, 1) /= nexternal + 1 .or. &
    size(config_mass, 2) /= lmaxconfigs .or. &
    size(config_mass, 3) /= fks_configs + 1) then
  call fail_genps_state('configuration masses have inconsistent bounds')
end if
if (any(shape(config_width) /= shape(config_mass))) then
  call fail_genps_state('configuration widths have inconsistent bounds')
end if
if (size(config_forest, 1) /= 2 .or. &
    size(config_forest, 2) /= max_branch .or. &
    size(config_forest, 3) /= lmaxconfigs .or. &
    size(config_forest, 4) /= fks_configs + 1) then
  call fail_genps_state('configuration forest has inconsistent bounds')
end if
if (size(config_tree, 1) /= 2 .or. &
    size(config_tree, 2) /= max_branch) then
  call fail_genps_state('configuration tree has inconsistent bounds')
end if
if (size(cnt_momenta, 1) /= 4 .or. &
    size(cnt_momenta, 2) /= nexternal .or. &
    size(cnt_momenta, 3) /= 5) then
  call fail_genps_state('counterevent momenta have inconsistent bounds')
end if
if (size(cnt_weight) /= 5 .or. size(cnt_psweight) /= 5 .or. &
    size(cnt_jacobian) /= 5) then
  call fail_genps_state('counterevent weights have inconsistent bounds')
end if
if (size(born_tree, 1) /= 2 .or. &
    size(born_tree, 2) /= max_branch) then
  call fail_genps_state('Born tree has inconsistent bounds')
end if
if (size(born_momenta, 1) /= 4 .or. &
    size(born_momenta, 2) /= nexternal - 1 .or. &
    any(shape(born_lab_momenta) /= shape(born_momenta)) .or. &
    any(shape(born_coll_momenta) /= shape(born_momenta)) .or. &
    any(shape(born_norad_momenta) /= shape(born_momenta))) then
  call fail_genps_state('Born momenta have inconsistent bounds')
end if
if (size(event_momenta, 1) /= 4 .or. &
    size(event_momenta, 2) /= nexternal) then
  call fail_genps_state('event momenta have inconsistent bounds')
end if
if (size(cbw_mass_state, 1) /= 3 .or. &
    size(cbw_mass_state, 2) /= nexternal .or. &
    any(shape(cbw_width_state) /= shape(cbw_mass_state)) .or. &
    size(cbw_state) /= nexternal .or. &
    size(cbw_level_state) /= nexternal) then
  call fail_genps_state('conflicting-BW state has inconsistent bounds')
end if
if (size(particle_masses) /= nexternal) then
  call fail_genps_state('particle masses have inconsistent bounds')
end if
if (size(schannel_masses) /= 2 * nexternal + 1) then
  call fail_genps_state('s-channel masses have inconsistent bounds')
end if
end subroutine validate_bound_genps_state


subroutine require_genps_state()
implicit none

if (.not. genps_state_initialized) then
  call fail_genps_state('module state has not been initialized')
end if
if (.not. associated(config_mass) .or. &
    .not. associated(config_width) .or. &
    .not. associated(config_forest) .or. &
    .not. associated(config_tree) .or. &
    .not. associated(config_index) .or. &
    .not. associated(cnt_momenta) .or. &
    .not. associated(cnt_weight) .or. &
    .not. associated(cnt_psweight) .or. &
    .not. associated(cnt_jacobian) .or. &
    .not. associated(born_tree) .or. &
    .not. associated(born_ns_channel) .or. &
    .not. associated(born_nt_channel) .or. &
    .not. associated(born_onebody_index) .or. &
    .not. associated(born_nbranch) .or. &
    .not. associated(born_one_body) .or. &
    .not. associated(born_momenta) .or. &
    .not. associated(born_lab_momenta) .or. &
    .not. associated(born_coll_momenta) .or. &
    .not. associated(born_norad_momenta) .or. &
    .not. associated(event_momenta) .or. &
    .not. associated(cbw_mass_state) .or. &
    .not. associated(cbw_width_state) .or. &
    .not. associated(cbw_level_max_state) .or. &
    .not. associated(cbw_state) .or. &
    .not. associated(cbw_level_state) .or. &
    .not. associated(particle_masses) .or. &
    .not. associated(schannel_masses) .or. &
    .not. allocated(saved_particle_masses)) then
  call fail_genps_state('module state is incomplete')
end if
end subroutine require_genps_state


subroutine fail_genps_state(message)
implicit none
character(len=*), intent(in) :: message

write (*,*) 'genps_fks: ', trim(message)
stop 1
end subroutine fail_genps_state

end module genps_fks
