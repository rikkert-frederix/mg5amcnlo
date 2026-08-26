module genps_born
  use process_dimensions, only: nexternal, nincoming, max_particles, &
                                max_branch, lmaxconfigs, fks_configs, &
                                validate_process_dimensions
  use run_state
  use kin_functions_module, only: dot => dot_impl
  use fks_singular_module, only: phspncheck_born
  implicit none
  private

  public :: born_phase_space
  public :: generate_born_phase_space
  public :: initialize_genps_born_state
  public :: invalidate_born_phase_space
  public :: phase_space_lambda

  type :: born_phase_space
    logical :: valid = .false.
    double precision :: xjac = -222d0
    double precision :: xpswgt = 1d0
    double precision :: stot = 0d0
    double precision :: shat = 0d0
    double precision :: sqrtshat = 0d0
    double precision :: tau = 0d0
    double precision :: ycm = 0d0
    double precision :: ycmhat = 0d0
    double precision :: xbjrk(2) = 0d0
    double precision, pointer :: masses(:) => null()
    double precision, pointer :: external_masses(:) => null()
  end type born_phase_space

! Python writes the dimensions of these historical COMMON blocks for each
! subprocess. The fixed-form bridge owns that generated storage; these
! pointers let this module use the same live objects.
  double precision, pointer :: config_mass(:, :, :) => null()
  double precision, pointer :: config_width(:, :, :) => null()
  integer, pointer :: config_forest(:, :, :, :) => null()
  integer, pointer :: config_tree(:, :) => null()
  integer, pointer :: config_index => null()

  integer, pointer :: born_tree(:, :) => null()
  integer, pointer :: born_ns_channel => null()
  integer, pointer :: born_nt_channel => null()
  integer, pointer :: born_onebody_index => null()
  integer, pointer :: born_nbranch => null()
  logical, pointer :: born_one_body => null()

  double precision, pointer :: born_momenta(:, :) => null()
  double precision, pointer :: born_lab_momenta(:, :) => null()

  double precision, pointer :: cbw_mass_state(:, :) => null()
  double precision, pointer :: cbw_width_state(:, :) => null()
  integer, pointer :: cbw_level_max_state => null()
  integer, pointer :: cbw_state(:) => null()
  integer, pointer :: cbw_level_state(:) => null()

  double precision, pointer :: particle_masses(:) => null()
  double precision, pointer :: schannel_masses(:) => null()

  double precision, allocatable, target :: saved_particle_masses(:)
  double precision, allocatable, target :: saved_external_masses(:)
  integer :: saved_configuration = 0
  logical :: first_configuration = .true.
  double precision :: saved_stot = 0d0
  double precision :: saved_initial_mass = 0d0
  double precision :: saved_final_mass = 0d0
  logical :: born_state_initialized = .false.

contains

  subroutine initialize_genps_born_state(config_mass_in, config_width_in, &
                                         config_forest_in, config_tree_in, config_index_in, &
                                         born_tree_in, born_ns_in, born_nt_in, &
                                         born_onebody_in, born_nbranch_in, born_one_body_in, &
                                         born_momenta_in, born_lab_momenta_in, &
                                         cbw_mass_in, cbw_width_in, cbw_level_max_in, &
                                         cbw_in, cbw_level_in, particle_masses_in, &
                                         schannel_masses_in)
    implicit none
    double precision, target, intent(inout) :: config_mass_in(-nexternal:, 1:, 0:)
    double precision, target, intent(inout) :: config_width_in(-nexternal:, 1:, 0:)
    integer, target, intent(inout) :: config_forest_in(1:, -max_branch:, 1:, 0:)
    integer, target, intent(inout) :: config_tree_in(1:, -max_branch:)
    integer, target, intent(inout) :: config_index_in
    integer, target, intent(inout) :: born_tree_in(1:, -max_branch:)
    integer, target, intent(inout) :: born_ns_in, born_nt_in
    integer, target, intent(inout) :: born_onebody_in, born_nbranch_in
    logical, target, intent(inout) :: born_one_body_in
    double precision, target, intent(inout) :: born_momenta_in(0:, 1:)
    double precision, target, intent(inout) :: born_lab_momenta_in(0:, 1:)
    double precision, target, intent(inout) :: cbw_mass_in(-1:, -nexternal:)
    double precision, target, intent(inout) :: cbw_width_in(-1:, -nexternal:)
    integer, target, intent(inout) :: cbw_level_max_in
    integer, target, intent(inout) :: cbw_in(-nexternal:)
    integer, target, intent(inout) :: cbw_level_in(-nexternal:)
    double precision, target, intent(inout) :: particle_masses_in(1:)
    double precision, target, intent(inout) :: schannel_masses_in(-nexternal:)

    call validate_process_dimensions()
    config_mass => config_mass_in
    config_width => config_width_in
    config_forest => config_forest_in
    config_tree => config_tree_in
    config_index => config_index_in
    born_tree => born_tree_in
    born_ns_channel => born_ns_in
    born_nt_channel => born_nt_in
    born_onebody_index => born_onebody_in
    born_nbranch => born_nbranch_in
    born_one_body => born_one_body_in
    born_momenta => born_momenta_in
    born_lab_momenta => born_lab_momenta_in
    cbw_mass_state => cbw_mass_in
    cbw_width_state => cbw_width_in
    cbw_level_max_state => cbw_level_max_in
    cbw_state => cbw_in
    cbw_level_state => cbw_level_in
    particle_masses => particle_masses_in
    schannel_masses => schannel_masses_in

    call validate_bound_born_state()

    if (.not. allocated(saved_particle_masses)) then
      allocate (saved_particle_masses(-max_branch:max_particles))
    else if (lbound(saved_particle_masses, 1) /= -max_branch .or. &
             ubound(saved_particle_masses, 1) /= max_particles) then
      call fail_born_state('saved mass storage has inconsistent bounds')
    end if
    if (.not. allocated(saved_external_masses)) then
      allocate (saved_external_masses(nexternal - 1))
    else if (size(saved_external_masses) /= nexternal - 1) then
      call fail_born_state('saved external mass storage has inconsistent bounds')
    end if
    born_state_initialized = .true.
  end subroutine initialize_genps_born_state

  subroutine generate_born_phase_space(ndim, iconfig, x, point)
    implicit none
    integer, intent(in) :: ndim, iconfig
    double precision, intent(inout) :: x(99)
    type(born_phase_space), intent(out) :: point

    double precision :: qmass(-nexternal:0), qwidth(-nexternal:0)
    double precision :: s(-max_branch:max_particles)
    integer :: i, j
    integer :: i_fks, j_fks
    common/fks_indices/i_fks, j_fks
    integer :: iconfig0
    common/ciconfig0/iconfig0
    integer :: this_config
    common/to_mconfigs/this_config

    call require_born_state()

    point%valid = .false.
    point%xjac = 1d0
    point%xpswgt = 1d0
    point%stot = 0d0
    point%shat = 0d0
    point%sqrtshat = 0d0
    point%tau = 0d0
    point%ycm = 0d0
    point%ycmhat = 0d0
    point%xbjrk = 0d0
    point%masses => saved_particle_masses
    point%external_masses => saved_external_masses
    saved_external_masses = 0d0

    this_config = iconfig
    config_index = iconfig
    iconfig0 = iconfig
    do i = -max_branch, -1
      do j = 1, 2
        config_tree(j, i) = config_forest(j, i, iconfig, 0)
      end do
    end do
    do i = -nexternal, 0
      qmass(i) = config_mass(i, iconfig, 0)
      qwidth(i) = config_width(i, iconfig, 0)
    end do

    call set_tau_min()

    do i = 1, nexternal - 1
      if (i < i_fks) then
        saved_particle_masses(i) = particle_masses(i)
      else
        saved_particle_masses(i) = particle_masses(i + 1)
      end if
    end do

    if (first_configuration .or. iconfig0 /= saved_configuration) then
      if (nincoming == 2) then
        saved_stot = 4d0*ebeam(1)*ebeam(2)
      else
        saved_stot = particle_masses(1)**2
      end if

      saved_initial_mass = sum(saved_particle_masses(1:nincoming))
      saved_final_mass = sum(saved_particle_masses(nincoming + 1:nexternal - 1))
      if (saved_stot < max(saved_final_mass, saved_initial_mass)**2) then
        write (*, *) 'Fatal error #0 in genps_born: insufficient collider energy'
        stop 1
      end if

      call fill_genmom_born_commons(config_tree, saved_particle_masses)
      first_configuration = .false.
      saved_configuration = iconfig0
    end if

    point%stot = saved_stot
    call generate_tau_y_wrapper(qmass, qwidth, saved_final_mass, saved_stot, &
                                x(ndim - 4:ndim - 3), point%tau, point%ycm, &
                                point%ycmhat, point%xjac)
    if (point%xjac < 0d0) then
      call invalidate_born_phase_space()
      return
    end if

    point%xbjrk(1) = sqrt(point%tau)*exp(point%ycm)
    point%xbjrk(2) = sqrt(point%tau)*exp(-point%ycm)
    if (.not. born_one_body) then
      point%shat = point%tau*saved_stot
      point%sqrtshat = sqrt(point%shat)
    else
      point%shat = saved_final_mass**2
      point%sqrtshat = saved_final_mass
    end if

    call generate_momenta_born(x, point%shat, point%sqrtshat, &
                               saved_final_mass, saved_particle_masses, s, &
                               qmass, qwidth, saved_external_masses, &
                               point%xpswgt, point%xjac)
    if (point%xjac < 0d0) then
      call invalidate_born_phase_space()
      return
    end if
    point%valid = .true.
  end subroutine generate_born_phase_space

  subroutine invalidate_born_phase_space()
    implicit none

    if (.not. associated(born_momenta)) then
      call fail_born_state('Born momenta are not bound')
    end if
    born_momenta(0, 1) = -99d0
  end subroutine invalidate_born_phase_space

  subroutine generate_tau_y_wrapper( &
    & qmass, qwidth, totmass, stot, rndx, tau_born, ycm_born, ycmhat, xjac)
! generates tau and y, calling the functions that correpsond to the
! case at hand
    implicit none

    double precision qmass(-nexternal:0), qwidth(-nexternal:0)
    double precision totmass, stot
    double precision rndx(2)
    double precision tau_born, ycm_born, ycmhat, xjac
!
    integer ndim_dummy

    integer i_fks, j_fks
    common/fks_indices/i_fks, j_fks

    logical softtest, colltest
    common/sctests/softtest, colltest

    ndim_dummy = -1 ! this is actually not used anymore

    if (abs(lpp(1)) .gt. 1 .or. abs(lpp(2)) .gt. 1) then
      write (*, *) 'The fNLO template supports only lpp=0,+1,-1', lpp
      stop 1
    end if

    if (abs(lpp(1)) .ne. abs(lpp(2))) then
      write (*, *) 'Different beams not implemented', lpp
      stop 1
    end if

    if (abs(lpp(1)) .ge. 1 .and. abs(lpp(2)) .ge. 1 .and. &
      & .not. (softtest .or. colltest)) then
! x(ndim-1) -> tau_cnt(0); x(ndim) -> ycm_cnt(0)
!  rndx(1) -> tau; rndx(2) -> ycm
      if (born_one_body) then
! tau is fixed by the mass of the final state particle
        call compute_tau_one_body(totmass, stot, tau_born, xjac)
      else
        if (born_nt_channel .eq. 0 .and. &
          & qwidth(-born_ns_channel - 1) .ne. 0.d0 .and. &
          & cbw_state(-born_ns_channel - 1) .ne. 2) then
! Generate tau according to a Breit-Wiger function
          call generate_tau_BW(stot, ndim_dummy, rndx(1), qmass( &
            & -born_ns_channel - 1), qwidth(-born_ns_channel - 1), &
            & cbw_state(-born_ns_channel - 1), &
            & cbw_mass_state(-1:1, -born_ns_channel - 1), &
            & cbw_width_state(-1:1, -born_ns_channel - 1), tau_born, xjac)
        else
!     not a Breit Wigner
          call generate_tau(stot, ndim_dummy, rndx(1), tau_born, xjac)
        end if
      end if

! Generate the rapditity of the Born system
      call generate_y(tau_born, rndx(2), ycm_born, ycmhat, xjac)
    elseif (abs(lpp(1)) .ge. 1 .and. &
      & .not. (softtest .or. colltest)) then
      write (*, *) 'Option x1 not implemented in one_tree'
      stop
    elseif (abs(lpp(2)) .ge. 1 .and. &
      & .not. (softtest .or. colltest)) then
      write (*, *) 'Option x2 not implemented in one_tree'
      stop
    else
! No PDFs (also use fixed energy when performing tests)
      call compute_tau_y_epem(j_fks, born_one_body, totmass, stot, &
        & tau_born, ycm_born, ycmhat)
      if (j_fks .le. nincoming .and. .not. (softtest .or. colltest)) then
        write (*, *) 'Process has incoming j_fks, but fixed shat: '// &
          & 'not allowed for processes generated at NLO.'
        stop 1
      end if
    end if

    return
  end subroutine generate_tau_y_wrapper

  subroutine generate_momenta_born(x, shat_born, sqrtshat_born, totmass, &
    & m, s, &
    & qmass, qwidth, m_born, xpswgt0, xjac0)
! generate the momenta for the reduced born system
    implicit none

    double precision x(99), shat_born, sqrtshat_born, totmass
    double precision S(-max_branch:max_particles), M(-max_branch:max_particles)
    double precision xpswgt0, xjac0
    double precision qmass(-nexternal:0), qwidth(-nexternal:0), &
      & m_born(nexternal - 1)

    logical pass
    double precision pb(0:3, -max_branch:nexternal - 1), p_born_CHECK(0:3, nexternal - 1)
    integer i, j

    pass = .true.

! Generate the momenta for the initial state of the Born system
    if (nincoming .eq. 2) then
      call mom2cx(sqrtshat_born, m(1), m(2), 1d0, 0d0, pb(0, 1), pb(0, 2))
    else
      pb(0, 1) = sqrtshat_born
      do i = 1, 2
        pb(i, 1) = 0d0
      end do
    end if
    s(-born_nbranch) = shat_born
    m(-born_nbranch) = sqrtshat_born
    pb(0, -born_nbranch) = m(-born_nbranch)
    pb(1, -born_nbranch) = 0d0
    pb(2, -born_nbranch) = 0d0
    pb(3, -born_nbranch) = 0d0
!
! Generate Born-level momenta
!
! Start by generating all the invariant masses of the s-channels
    call generate_inv_mass_sch(born_ns_channel, born_tree, m, sqrtshat_born &
      & , totmass, qwidth, qmass, cbw_state, cbw_mass_state, &
      & cbw_width_state, s, x, xjac0 &
      & , pass)
    if (.not. pass) then
      xjac0 = -139
      return
    end if
! If only s-channels, also set the p1+p2 s-channel
    if (born_nt_channel .eq. 0 .and. nincoming .eq. 2) then
      s(-born_nbranch + 1) = s(-born_nbranch)
      m(-born_nbranch + 1) = m(-born_nbranch)       !Basic s-channel has s_hat
      pb(0, -born_nbranch + 1) = m(-born_nbranch + 1)!and 0 momentum
      pb(1, -born_nbranch + 1) = 0d0
      pb(2, -born_nbranch + 1) = 0d0
      pb(3, -born_nbranch + 1) = 0d0
    end if
!
!     Next do the T-channel branchings
!
    if (born_nt_channel .ne. 0) then
      call generate_t_channel_branchings(born_ns_channel, born_nbranch, born_tree &
        & , m, s, x, pb, xjac0, xpswgt0, pass)
      if (.not. pass) then
        xjac0 = -140
        return
      end if
    end if
!
!     Now generate momentum for all intermediate and final states
!     being careful to calculate from more massive to less massive states
!     so the last states done are the final particle states.
!
    call fill_born_momenta(born_nbranch, born_nt_channel, born_one_body, &
      & born_onebody_index, x, born_tree, m, s, pb, xjac0, xpswgt0, pass)
    if (.not. pass) then
      xjac0 = -141
      return
    end if
!
!  Now I have the Born momenta
!
    do i = 1, nexternal - 1
      do j = 0, 3
        born_lab_momenta(j, i) = pb(j, i)
        p_born_CHECK(j, i) = pb(j, i)
      end do
      m_born(i) = m(i)
    end do
    call phspncheck_born(sqrtshat_born, m_born, p_born_CHECK, pass)
    if (.not. pass) then
      xjac0 = -142
      return
    end if

    do i = 1, nexternal - 1
      do j = 0, 3
        born_momenta(j, i) = born_lab_momenta(j, i)
      end do
    end do

    return
  end subroutine generate_momenta_born

  subroutine fill_genmom_born_commons(itree, m)
    implicit none
! arguments
    integer itree(2, -max_branch:-1)
    double precision M(-max_branch:max_particles)

    born_tree(:, :) = itree(:, :)

    born_nbranch = nexternal - 3 ! nexternal is for n+1-body, while itree uses n-body

! Determine number of s- and t-channel branches, at this point it
! includes the s-channel p1+p2
    born_ns_channel = 1
    do while (itree(1, -born_ns_channel) .ne. 1 .and. &
      & itree(1, -born_ns_channel) .ne. 2 .and. &
      & born_ns_channel .lt. born_nbranch)
      m(-born_ns_channel) = 0d0
      born_ns_channel = born_ns_channel + 1
    end do
    born_ns_channel = born_ns_channel - 1
    born_nt_channel = born_nbranch - born_ns_channel - 1
! If no t-channles, ns_channels is one less, because we want to exclude
! the s-channel p1+p2
    if (born_nt_channel .eq. 0 .and. nincoming .eq. 2) then
      born_ns_channel = born_ns_channel - 1
    end if
! Set one_body to true if it's a 2->1 process at the Born (i.e. 2->2 for the n+1-body)
    if ((nexternal - nincoming) .eq. 2) then
      born_one_body = .true.
      born_onebody_index = nexternal - 1
      born_ns_channel = 0
      born_nt_channel = 0
    elseif ((nexternal - nincoming) .gt. 2) then
      born_one_body = .false.
    else
      write (*, *) 'Error #1 in genps_born.f', nexternal, nincoming
      stop
    end if

    return
  end subroutine fill_genmom_born_commons

  subroutine gentcms(pa, pb, t, phi, m1, m2, p1, pr, jac)
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
    double precision t, phi, m1, m2               !inputs
    double precision pa(0:3), pb(0:3), jac
    double precision p1(0:3), pr(0:3)           !outputs
!
!     local
!
    double precision ptot(0:3), E_acms, p_acms, pa_cms(0:3)
    double precision esum, ed, pp, md2, ma2, pt, ptotm(0:3)
    integer i
!
!     External
!
!-----
!  Begin Code
!-----
    do i = 0, 3
      ptot(i) = pa(i) + pb(i)
      if (i .gt. 0) then
        ptotm(i) = -ptot(i)
      else
        ptotm(i) = ptot(i)
      end if
    end do
    ma2 = dot(pa, pa)
!
!     determine magnitude of p1 in cms frame (from dhelas routine mom2cx)
!
    ESUM = sqrt(max(0d0, dot(ptot, ptot)))
    if (esum .eq. 0d0) then
      jac = -8d0             !Failed esum must be > 0
      return
    end if
    MD2 = (M1 - M2)*(M1 + M2)
    ED = MD2/ESUM
    if (M1*M2 .eq. 0.) then
      PP = (ESUM - abs(ED))*0.5d0
    else
      PP = (MD2/ESUM)**2 - 2.0d0*(M1**2 + M2**2) + ESUM**2
      if (pp .gt. 0) then
        PP = sqrt(pp)*0.5d0
      else
        write (*, *) 'Warning #12 in genps_born.f', pp
        jac = -1
        return
      end if
    end if
!
!     Energy of pa in pa+pb cms system
!
    call boostx(pa, ptotm, pa_cms)
    E_acms = pa_cms(0)
    p_acms = dsqrt(pa_cms(1)**2 + pa_cms(2)**2 + pa_cms(3)**2)
!
    p1(0) = max((ESUM + ED)*0.5d0, 0.d0)
    p1(3) = -(m1*m1 + ma2 - t - 2d0*p1(0)*E_acms)/(2d0*p_acms)
    pt = dsqrt(max(pp*pp - p1(3)*p1(3), 0d0))
    p1(1) = pt*cos(phi)
    p1(2) = pt*sin(phi)
!
    call rotxxx(p1, pa_cms, p1)          !Rotate back to pa_cms frame
    call boostx(p1, ptot, p1)            !boost back to lab fram
    do i = 0, 3
      pr(i) = pa(i) - p1(i)               !Return remainder of momentum
    end do
  end subroutine gentcms

  double precision function phase_space_lambda(S, MA2, MB2)
    implicit none
!****************************************************************************
!     THIS IS THE phase_space_lambda FUNCTION FROM VERNONS BOOK COLLIDER PHYSICS P 662
!     MA2 AND MB2 ARE THE MASS SQUARED OF THE FINAL STATE PARTICLES
!     2-D PHASE SPACE = .5*PI*SQRT(1.,MA2/S^2,MB2/S^2)*(D(OMEGA)/4PI)
!****************************************************************************
    double precision MA2, MB2, S, tiny, tmp, rat
    parameter(tiny=1.d-8)
!
    tmp = S**2 + MA2**2 + MB2**2 - 2d0*S*MA2 - 2d0*MA2*MB2 - 2d0*S*MB2
    if (tmp .le. 0.d0) then
      if (ma2 .lt. 0.d0 .or. mb2 .lt. 0.d0) then
        write (6, *) 'Error #1 in function phase_space_lambda:', s, ma2, mb2
        stop
      end if
      rat = 1 - (sqrt(ma2) + sqrt(mb2))/s
      if (rat .gt. -tiny) then
        tmp = 0.d0
      else
        write (6, *) 'Error #2 in function phase_space_lambda:', s, ma2, mb2, rat
      end if
    end if
    phase_space_lambda = tmp
    return
  end function phase_space_lambda

  subroutine YMINMAX(X, Z, U, V, W, YMIN, YMAX)
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
    parameter(tiny=1d-199)
!
!     Arguments
!
    double precision x, z, u, v, w                ! inputs
    double precision ymin, ymax                !output
!
!     Local
!
    double precision y1, y2, yr, ysqr
!
!-----
!  Begin Code
!-----
    ysqr = phase_space_lambda(x, u, v)*phase_space_lambda(x, w, z)
    if (ysqr .ge. 0d0) then
      yr = dsqrt(ysqr)
    else
      print *, 'Error in yminymax sqrt(-x)', phase_space_lambda(x, u, v), &
        phase_space_lambda(x, w, z)
      yr = 0d0
    end if
    y1 = u + w - .5d0*((x + u - v)*(x + w - z) - yr)/(x + tiny)
    y2 = u + w - .5d0*((x + u - v)*(x + w - z) + yr)/(x + tiny)
    ymin = min(y1, y2)
    ymax = max(y1, y2)
  end subroutine YMINMAX

  subroutine compute_tau_one_body(totmass, stot, tau, jac)
    implicit none
    double precision totmass, stot, tau, jac, roH
    roH = totmass**2/stot
    tau = roH
! Jacobian due to delta() of tau_born
    jac = jac*2*totmass/stot
    return
  end subroutine compute_tau_one_body

  subroutine generate_tau_BW(stot, idim, x, mass, width, cBW, BWmass &
    & , BWwidth, tau, jac)
    implicit none
    integer cBW, idim
    double precision stot, x, tau, jac, mass, width, BWmass(-1:1), BWwidth( &
      & -1:1), s_mass, s
    double precision smax, smin
    double precision tau_Born_lower_bound, tau_lower_bound_resonance &
      & , tau_lower_bound
    common/ctau_lower_bound/tau_Born_lower_bound &
      & , tau_lower_bound_resonance, tau_lower_bound
    if (cBW .eq. 1 .and. width .gt. 0d0 .and. BWwidth(1) .gt. 0d0) then
      smin = tau_Born_lower_bound*stot
      smax = stot
      s_mass = smin
      call trans_x(5, idim, x, smin, smax, s_mass, mass, width, BWmass( &
        & -1), BWwidth(-1), jac, s)
      tau = s/stot
      jac = jac/stot
    else
      smin = tau_Born_lower_bound*stot
      smax = stot
      s_mass = smin
      call trans_x(3, idim, x, smin, smax, s_mass, mass, width, BWmass( &
        & -1), BWwidth(-1), jac, s)
      tau = s/stot
      jac = jac/stot
    end if
    return
  end subroutine generate_tau_BW

  subroutine generate_tau(stot, idim, x, tau, jac)
    implicit none
    integer idim
    double precision x, tau, jac, smin, smax, s_mass, s, tiny, dum, dum3(-1:1) &
      & , stot
    parameter(tiny=1d-8)
    double precision tau_Born_lower_bound, tau_lower_bound_resonance &
      & , tau_lower_bound
    common/ctau_lower_bound/tau_Born_lower_bound &
      & , tau_lower_bound_resonance, tau_lower_bound
    smin = tau_born_lower_bound*stot
    smax = stot
    s_mass = tau_lower_bound_resonance*stot
    if (s_mass .gt. smin*(1d0 + tiny)) then
      call trans_x(2, idim, x, smin, smax, s_mass, dum, dum &
        & , dum3, dum3, jac, s)
    elseif (abs(s_mass - smin) .lt. tiny*smin) then
      call trans_x(7, idim, x, smin, smax, s_mass, dum, dum &
        & , dum3, dum3, jac, s)
    else
      write (*, *) 'ERROR #39 in genps_born.f', s_mass, smin, smax
      jac = -1d0
    end if
    tau = s/stot
    jac = jac/stot
    return
  end subroutine generate_tau

  subroutine generate_y(tau, x, ycm, ycmhat, jac)
    implicit none
    double precision tau, x, ycm, jac
    double precision ylim, ycmhat
    ylim = -0.5d0*log(tau)
    ycmhat = 2*x - 1
    ycm = ylim*ycmhat
    jac = jac*ylim*2
    return
  end subroutine generate_y

  subroutine compute_tau_y_epem(j_fks, one_body, fksmass, &
    & stot, tau, ycm, ycmhat)
    implicit none
    integer j_fks
    logical one_body
    double precision fksmass, stot, tau, ycm, ycmhat
    if (j_fks .le. nincoming) then
! This should never happen in normal integration: when no PDFs, j_fks
! cannot be initial state (but needed for testing). If tau set to one,
! integration range in xi_i_fks will be zero, so lower it artificially
! when too large
      if (one_body) then
        tau = fksmass**2/stot
      else
        tau = max((0.85d0)**2, fksmass**2/stot)
      end if
      ycm = 0.d0
    else
! For e+e- collisions, set tau to one and y to zero
      tau = 1.d0
      ycm = 0.d0
    end if
    ycmhat = 0.d0
    return
  end subroutine compute_tau_y_epem

  subroutine generate_inv_mass_sch(ns_channel, itree, m, sqrtshat_born &
    & , totmass, qwidth, qmass, cBW, cBW_mass, cBW_width, s, x, xjac0, pass)
    implicit none
    integer ns_channel
    double precision qmass(-nexternal:0), qwidth(-nexternal:0)
    double precision M(-max_branch:max_particles), x(99)
    double precision s(-max_branch:max_particles)
    double precision sqrtshat_born, totmass, xjac0
    integer itree(2, -max_branch:-1)
    integer i, ii, order(-nexternal:0)
    double precision smin, smax, totalmass
    logical pass
    integer cBW(-nexternal:-1)
    double precision cBW_mass(-1:1, -nexternal:-1), cBW_width(-1:1, &
      & -nexternal:-1)
    pass = .true.
    totalmass = totmass
    do ii = -1, -ns_channel, -1
! Randomize the order with which to generate the s-channel masses:
      call sChan_order(ns_channel, order)
      i = order(ii)
!     Generate invariant masses for all s-channel branchings of the Born
      smin = (m(itree(1, i)) + m(itree(2, i)))**2
      smax = (sqrtshat_born - totalmass + sqrt(smin))**2
      if (smax .lt. smin .or. smax .lt. 0.d0 .or. smin .lt. 0.d0) then
        write (*, *) 'Error #13 in genps_born.f'
        write (*, *) smin, smax, i
        stop
      end if
      call generate_si(i, smin, smax, s, cBW, cBW_width, cBW_mass, qmass &
        & , qwidth, x, xjac0, schannel_masses)
! If numerical inaccuracy, quit loop
      if (xjac0 .lt. 0d0) then
        if ((xjac0 .gt. -400d0 .or. xjac0 .le. -500d0) .and. &
          & xjac0 .ne. 0d0) then
          write (*, *) 'WARNING #31 in genps_born.f', i, s(i), smin, smax &
            & , xjac0
        end if
        xjac0 = -6
        pass = .false.
        return
      end if
      if (s(i) .lt. smin) then
        write (*, *) 'WARNING #32 in genps_born.f', i, s(i), smin, smax, x( &
          & -i)
        xjac0 = -5
        pass = .false.
        return
      end if
!
!     fill masses, update totalmass
!
      m(i) = sqrt(s(i))
      totalmass = totalmass + m(i) - &
        & m(itree(1, i)) - m(itree(2, i))
      if (totalmass .gt. sqrtshat_born) then
        write (*, *) 'WARNING #33 in genps_born.f', i, totalmass &
          & , sqrtshat_born, s(i)
        xjac0 = -4
        pass = .false.
        return
      end if
    end do
    return
  end subroutine generate_inv_mass_sch

  subroutine generate_si(i, smin, smax, s, cBW, cBW_width, cBW_mass, qmass &
    & , qwidth, x, xjac0, s_mass)
    implicit none
    integer i
    double precision smin, smax, s(-max_branch:max_particles), qwidth( &
      & -nexternal:0), qmass(-nexternal:0), cBW_width(-1:1, -nexternal: &
      & -1), cBW_mass(-1:1, -nexternal:-1), xjac0, x(99), s_mass( &
      & -nexternal:nexternal)
    integer cBW(-nexternal:-1)
! Choose the appropriate s given our constraints smin,smax
    if (qwidth(i) .ne. 0.d0 .and. cBW(i) .ne. 2) then
! Breit Wigner
      if (cBW(i) .eq. 1 .and. &
        & cBW_width(1, i) .gt. 0d0 .and. cBW_width(-1, i) .gt. 0d0) then
!     conflicting BW on both sides
        call trans_x(6, -i, x(-i), smin, smax, s_mass(i), qmass(i) &
          & , qwidth(i), cBW_mass(-1, i), cBW_width(-1, i), xjac0, s(i))
      elseif (cBW(i) .eq. 1 .and. cBW_width(1, i) .gt. 0d0) then
!     conflicting BW with alternative mass larger
        call trans_x(5, -i, x(-i), smin, smax, s_mass(i), qmass(i) &
          & , qwidth(i), cBW_mass(-1, i), cBW_width(-1, i), xjac0, s(i))
      elseif (cBW(i) .eq. 1 .and. cBW_width(-1, i) .gt. 0d0) then
!     conflicting BW with alternative mass smaller
        call trans_x(4, -i, x(-i), smin, smax, s_mass(i), qmass(i) &
          & , qwidth(i), cBW_mass(-1, i), cBW_width(-1, i), xjac0, s(i))
      else
!     normal BW
        call trans_x(3, -i, x(-i), smin, smax, s_mass(i), qmass(i) &
          & , qwidth(i), cBW_mass(-1, i), cBW_width(-1, i), xjac0, s(i))
      end if
    else
! not a Breit Wigner
      if (smin .eq. 0d0 .and. s_mass(i) .eq. 0d0) then
!     no lower limit on invariant mass from cuts or final state masses:
!     use flat distribution
        call trans_x(1, -i, x(-i), smin, smax, s_mass(i), qmass(i) &
          & , qwidth(i), cBW_mass(-1, i), cBW_width(-1, i), xjac0, s(i))
      elseif (smin .ge. s_mass(i) .and. smin .gt. 0d0) then
!     A lower limit on smin, which is larger than lower limit from cuts
!     or masses. Use 1/x importance sampling
        call trans_x(7, -i, x(-i), smin, smax, s_mass(i), qmass(i) &
          & , qwidth(i), cBW_mass(-1, i), cBW_width(-1, i), xjac0, s(i))
      elseif (smin .lt. s_mass(i) .and. s_mass(i) .gt. 0d0) then
!     Use flat grid between smin and s_mass(i), and 1/x^nsamp above
!     s_mass(i)
        call trans_x(2, -i, x(-i), smin, smax, s_mass(i), qmass(i) &
          & , qwidth(i), cBW_mass(-1, i), cBW_width(-1, i), xjac0, s(i))
      else
        write (*, *) "ERROR in genps_born.f:"// &
          & " cannot set s-channel without BW", i, smin, s_mass(i)
        stop 1
      end if
    end if
    return
  end subroutine generate_si

  subroutine generate_t_channel_branchings(ns_channel, nbranch, itree &
    & , m, s, x, pb, xjac0, xpswgt0, pass)
! First we need to determine the energy of the remaining particles this
! is essentially in place of the cos(theta) degree of freedom we have
! with the s channel decay sequence
    implicit none
    double precision pi
    parameter(pi=3.1415926535897932d0)
    double precision xjac0, xpswgt0
    double precision M(-max_branch:max_particles), x(99)
    double precision s(-max_branch:max_particles)
    double precision pb(0:3, -max_branch:nexternal - 1)
    integer itree(2, -max_branch:-1)
    integer ns_channel, nbranch
    logical pass
!
    double precision totalmass, smin, smax, s1, ma2, mbq, m12, mnq, tmin, tmax &
      & , t, phi, dum, dum3(-1:1), s_m, tm, tiny
    parameter(tiny=1d-8)
    integer i, ibranch, idim
!
    pass = .true.
    totalmass = 0d0
    s_m = 0d0
    do ibranch = -ns_channel - 1, -nbranch, -1
      totalmass = totalmass + m(itree(2, ibranch))
      s_m = s_m + sqrt(schannel_masses(itree(2, ibranch)))
    end do
    m(-ns_channel - 1) = dsqrt(S(-nbranch))
!
! Choose invariant masses of the pseudoparticles obtained by taking together
! all final-state particles or pseudoparticles found from the current
! t-channel propagator down to the initial-state particle found at the end
! of the t-channel line.
    do ibranch = -ns_channel - 1, -nbranch + 2, -1
      totalmass = totalmass - m(itree(2, ibranch))
      smin = totalmass**2
      smax = (m(ibranch) - m(itree(2, ibranch)))**2
      if (smin .gt. smax) then
        xjac0 = -3d0
        pass = .false.
        return
      end if
      idim = (nbranch - 1 + (-ibranch)*2)
      s_m = s_m - sqrt(schannel_masses(itree(2, ibranch)))
      if (abs(smin - s_m**2) .lt. tiny) then
        call trans_x(1, idim, x(idim), smin, smax, s_m**2, dum &
          & , dum, dum3(-1), dum3(-1), xjac0, s1)
      else
        call trans_x(1, idim, x(idim), smin, smax, s_m**2, dum &
          & , dum, dum3(-1), dum3(-1), xjac0, s1)
      end if
      if (xjac0 .le. 0d0) then
        if ((xjac0 .gt. -400d0 .or. xjac0 .le. -500d0) .and. &
          & xjac0 .ne. 0d0) then
          write (*, *) 'WARNING #31a in genps_born.f', ibranch, s1 &
            & , smin, smax, s_m**2, xjac0
        end if
        xjac0 = -6
        pass = .false.
        return
      end if
      m(ibranch - 1) = sqrt(s1)
      if (m(ibranch - 1)**2 .lt. smin .or. m(ibranch - 1)**2 .gt. smax &
        & .or. m(ibranch - 1) .ne. m(ibranch - 1)) then
        xjac0 = -1d0
        pass = .false.
        return
      end if
    end do
!
! Set m(-nbranch) equal to the mass of the particle or pseudoparticle P
! attached to the vertex (P,t,p2), with t being the last t-channel propagator
! in the t-channel line, and p2 the incoming particle opposite to that from
! which the t-channel line starts
    m(-nbranch) = m(itree(2, -nbranch))
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
    do ibranch = -ns_channel - 1, -nbranch + 1, -1
      s1 = m(ibranch)**2    !Total mass available
      ma2 = m(2)**2
      mbq = dot(pb(0, itree(1, ibranch)), pb(0, itree(1, ibranch)))
      m12 = m(itree(2, ibranch))**2
      mnq = m(ibranch - 1)**2
      call yminmax(s1, m12, ma2, mbq, mnq, tmin, tmax)
      call trans_x(1, -ibranch, x(-ibranch), -tmax, -tmin, &
        & schannel_masses(ibranch) &
        & , dum, dum, dum3(-1), dum3(-1), xjac0, tm)
      if (xjac0 .le. 0d0) then
        if ((xjac0 .gt. -400d0 .or. xjac0 .le. -500d0) .and. &
          & xjac0 .ne. 0d0) then
          write (*, *) 'WARNING #31b in genps_born.f', ibranch, tm &
            & , -tmax, -tmin, xjac0
        end if
        xjac0 = -6
        pass = .false.
        return
      end if
      t = -tm
      if (t .lt. tmin .or. t .gt. tmax) then
        write (*, *) "WARNING #35 in genps_born.f", t, tmin, tmax
        xjac0 = -3d0
        pass = .false.
        return
      end if
      phi = 2d0*pi*x(nbranch + (-ibranch - 1)*2)
      xjac0 = xjac0*2d0*pi
! Finally generate the momentum. The call is of the form
! pa+pb -> p1+ p2; t=(pa-p1)**2;   pr = pa-p1
! gentcms(pa,pb,t,phi,m1,m2,p1,pr)
      call gentcms(pb(0, itree(1, ibranch)), pb(0, 2), t, phi, &
        & m(itree(2, ibranch)), m(ibranch - 1), pb(0, itree(2, ibranch)), &
        & pb(0, ibranch), xjac0)
!
      if (xjac0 .lt. 0d0) then
        write (*, *) 'Failed gentcms', ibranch, xjac0
        pass = .false.
        return
      end if
      xpswgt0 = xpswgt0/(4d0*dsqrt(phase_space_lambda(s1, ma2, mbq)))
    end do
! We need to get the momentum of the last external particle.  This
! should just be the sum of p(0,2) and the remaining momentum from our
! last t channel 2->2
    do i = 0, 3
      pb(i, itree(2, -nbranch)) = pb(i, -nbranch + 1) + pb(i, 2)
    end do
    return
  end subroutine generate_t_channel_branchings

  subroutine fill_born_momenta(nbranch, nt_channel, one_body, ionebody &
    & , x, itree, m, s, pb, xjac0, xpswgt0, pass)
    implicit none
    double precision pi
    parameter(pi=3.1415926535897932d0)
    integer nbranch, nt_channel, ionebody
    double precision M(-max_branch:max_particles), x(99)
    double precision s(-max_branch:max_particles)
    double precision pb(0:3, -max_branch:nexternal - 1)
    integer itree(2, -max_branch:-1)
    double precision xjac0, xpswgt0
    logical pass, one_body
!
    double precision one
    parameter(one=1d0)
    double precision costh, phi, xa2, xb2
    integer i, ix
    double precision vtiny
    parameter(vtiny=1d-12)
!
    pass = .true.
    do i = -nbranch + nt_channel + (nincoming - 1), -1
      ix = nbranch + (-i - 1)*2 + (2 - nincoming)
      if (nt_channel .eq. 0) ix = ix - 1
      costh = 2d0*x(ix) - 1d0
      phi = 2d0*pi*x(ix + 1)
      xjac0 = xjac0*4d0*pi
      xa2 = m(itree(1, i))*m(itree(1, i))/s(i)
      xb2 = m(itree(2, i))*m(itree(2, i))/s(i)
      if (m(itree(1, i)) + m(itree(2, i)) .ge. m(i)) then
        xjac0 = -8
        pass = .false.
        return
      end if
      xpswgt0 = xpswgt0*.5d0*PI*sqrt(phase_space_lambda(ONE, XA2, XB2))/(4.d0*PI)
      call mom2cx(m(i), m(itree(1, i)), m(itree(2, i)), costh, phi, &
        & pb(0, itree(1, i)), pb(0, itree(2, i)))
! If there is an extremely large boost needed here, skip the phase-space point
! because of numerical stabilities.
      if (dsqrt(abs(dot(pb(0, i), pb(0, i))))/pb(0, i) &
        & .lt. vtiny) then
        xjac0 = -81
        pass = .false.
        return
      else
        call boostm(pb(0, itree(1, i)), pb(0, i), m(i), pb(0, itree(1, i)))
        call boostm(pb(0, itree(2, i)), pb(0, i), m(i), pb(0, itree(2, i)))
      end if
    end do
!
!
! Special phase-space fix for the one_body
    if (one_body) then
! Factor due to the delta function in dphi_1
      xpswgt0 = pi/m(ionebody)
! Kajantie's normalization of phase space (compensated below in flux)
      xpswgt0 = xpswgt0/(2*pi)
      do i = 0, 3
        pb(i, 3) = pb(i, 1) + pb(i, 2)
      end do
    end if
    return
  end subroutine fill_born_momenta

  subroutine trans_x(itype, idim, x, smin, smax, s_mass, qmass, qwidth &
    & , cBW_mass, cBW_width, jac, s)
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
    integer itype, idim
    double precision x, smin, smax, s_mass, qmass, qwidth, cBW_mass(-1:1) &
      & , cBW_width(-1:1), jac, s
    double precision fract, A, B, C, bs(-1:1), maxi, mini
    integer j
!
    if (itype .eq. 1) then
!     flat transformation:
      A = smax - smin
      B = smin
      s = A*x + B
      jac = jac*A
    elseif (itype .eq. 2) then
      fract = 0.25d0
      if (s_mass .eq. 0d0) then
        write (*, *) 's_mass is zero', itype, idim
      end if
      if (x .lt. fract) then
!     flat transformation:
        if (s_mass .lt. smin) then
          jac = -421d0
          return
        end if
        maxi = min(s_mass, smax)
        A = (maxi - smin)/fract
        B = smin
        s = A*x + B
        jac = jac*A
      else
!     S=A/(B-x) transformation:
        if (s_mass .ge. smax) then
          jac = -422d0
          return
        end if
        mini = max(s_mass, smin)
        A = mini*smax*(1d0 - fract)/(smax - mini)
        B = (smax - fract*mini)/(smax - mini)
        s = A/(B - x)
        jac = jac*s**2/A
      end if
    elseif (itype .eq. 3) then
!     Normal Breit-Wigner, i.e.
!        \int_smin^smax ds g(s)/((s-qmass^2)^2-qmass^2*qwidth^2) =
!        \int_0^1 dx g(s(x))
      A = atan((qmass - smin/qmass)/qwidth)
      B = atan((qmass - smax/qmass)/qwidth)
      s = qmass*(qmass - qwidth*tan(A - (A - B)*x))
      jac = jac*qmass*qwidth*(A - B)/(cos(A - (A - B)*x))**2
    elseif (itype .eq. 4) then
!     Conflicting BW, with alternative mass smaller than current
!     mass. That is, we need to throw also many events at smaller masses
!     than the peak of the current BW. Split 'x' at 'bs(-1)', using a
!     flat distribution below the split, and a BW above the split.
      fract = 0.3d0
      bs(-1) = (cBW_mass(-1) - qmass)/ &
        & (qwidth + cBW_width(-1)) ! bs(-1) is negative here
      bs(-1) = qmass + bs(-1)*qwidth
      bs(-1) = bs(-1)**2
      if (x .lt. fract) then
        if (smin .gt. bs(-1)) then
          jac = -441d0
          return
        end if
        maxi = min(bs(-1), smax)
        A = (maxi - smin)/fract
        B = smin
        s = A*x + B
        jac = jac*A
      else
        if (smax .lt. bs(-1)) then
          jac = -442d0
          return
        end if
        mini = max(bs(-1), smin)
        A = atan((qmass - mini/qmass)/qwidth)
        B = atan((qmass - smax/qmass)/qwidth)
        C = ((1d0 - x)*A + (x - fract)*B)/(1d0 - fract)
        s = qmass*(qmass - qwidth*tan(C))
        jac = jac*qmass*qwidth*(A - B)/((cos(C))**2*(1d0 - fract))
      end if
    elseif (itype .eq. 5) then
!     Conflicting BW, with alternative mass larger than current
!     mass. That is, we need to throw also many events at larger masses
!     than the peak of the current BW. Split 'x' at 'bs(1)' and the
!     alternative mass. Use a BW below bs(1), a flat distribution
!     between bs(1) and the alternative mass, and a 1/x above the
!     alternative mass.
      fract = 0.35d0
      bs(1) = (cBW_mass(1) - qmass)/ &
        & (qwidth + cBW_width(1))
      bs(1) = qmass + bs(1)*qwidth
      bs(1) = bs(1)**2
      if (x .lt. fract) then
        if (smin .gt. bs(1)) then
          jac = -451d0
          return
        end if
        maxi = min(bs(1), smax)
        A = atan((qmass - smin/qmass)/qwidth)
        B = atan((qmass - maxi/qmass)/qwidth)
        C = ((B - A)*x + fract*A)/fract
        s = qmass*(qmass - qwidth*tan(C))
        jac = jac*qmass*qwidth*(A - B)/((cos(C))**2*fract)
      elseif (x .lt. 1d0 - fract) then
        if (smin .gt. cBW_mass(1)**2 .or. smax .lt. bs(1)) then
          jac = -452d0
          return
        end if
        maxi = min(cBW_mass(1)**2, smax)
        mini = max(bs(1), smin)
        A = (maxi - mini)/(1d0 - 2d0*fract)
        B = ((1d0 - fract)*mini - fract*maxi)/(1d0 - 2d0*fract)
        s = A*x + B
        jac = jac*A
      else
        if (smax .le. cBW_mass(1)**2) then
          jac = -453d0
          return
        end if
        mini = max(cBW_mass(1)**2, smin)
        A = mini*smax*fract/(smax - mini)
        B = (smax - (1d0 - fract)*mini)/(smax - mini)
        s = A/(B - x)
        jac = jac*s**2/A
      end if
    elseif (itype .eq. 6) then
      fract = 0.3d0
!     Conflicting BW on both sides. Use flat below bs(-1); BW between
!     bs(-1) and bs(1); flat between bs(1) and alternative mass; and 1/x
!     above alternative mass.
      do j = -1, 1, 2
        bs(j) = (cBW_mass(j) - qmass)/ &
          & (qwidth + cBW_width(j))
        bs(j) = qmass + bs(j)*qwidth
        bs(j) = bs(j)**2
      end do
      if (x .lt. fract) then
        if (smin .gt. bs(-1)) then
          jac = -461d0
          return
        end if
        maxi = min(bs(-1), smax)
        A = (maxi - smin)/fract
        B = smin
        s = A*x + B
        jac = jac*A
      elseif (x .lt. 1d0 - fract) then
        if (smin .gt. bs(1) .or. smax .lt. bs(-1)) then
          jac = -462d0
          return
        end if
        maxi = min(bs(1), smax)
        mini = max(bs(-1), smin)
        A = atan((qmass - mini/qmass)/qwidth)
        B = atan((qmass - maxi/qmass)/qwidth)
        C = ((1d0 - fract - x)*A + (x - fract)*B)/(1d0 - 2d0*fract)
        s = qmass*(qmass - qwidth*tan(C))
        jac = -jac*qmass*qwidth*(B - A)/((cos(C))**2*(1d0 - 2d0*fract))
      elseif (x .lt. 1d0 - fract/2d0) then
        if (smin .gt. cBW_mass(1)**2 .or. smax .lt. bs(1)) then
          jac = -463d0
          return
        end if
        maxi = min(cBW_mass(1)**2, smax)
        mini = max(bs(1), smin)
        A = 2d0*(maxi - mini)/fract
        B = 2d0*maxi - mini - 2d0*(maxi - mini)/fract
        s = A*x + B
        jac = jac*A
      else
        if (smax .le. cBW_mass(1)**2) then
          jac = -464d0
          return
        end if
        mini = max(cBW_mass(1)**2, smin)
        A = mini*smax*fract/(2d0*(smax - mini))
        B = (smax - (1d0 - fract/2d0)*mini)/(smax - mini)
        s = A/(B - x)
        jac = jac*s**2/A
      end if
    elseif (itype .eq. 7) then
!     S=A/(B-x) transformation:
      if (smin .le. 0d0) then
        jac = -471d0
        return
      end if
      A = smin*smax/(smax - smin)
      B = smax/(smax - smin)
      s = A/(B - x)
      jac = jac*s**2/A
    end if
    return
  end subroutine trans_x

  subroutine validate_bound_born_state()
    implicit none

    if (size(config_mass, 1) /= nexternal + 1 .or. &
        size(config_mass, 2) /= lmaxconfigs .or. &
        size(config_mass, 3) /= fks_configs + 1) then
      call fail_born_state('configuration masses have inconsistent bounds')
    end if
    if (any(shape(config_width) /= shape(config_mass))) then
      call fail_born_state('configuration widths have inconsistent bounds')
    end if
    if (size(config_forest, 1) /= 2 .or. &
        size(config_forest, 2) /= max_branch .or. &
        size(config_forest, 3) /= lmaxconfigs .or. &
        size(config_forest, 4) /= fks_configs + 1) then
      call fail_born_state('configuration forest has inconsistent bounds')
    end if
    if (size(config_tree, 1) /= 2 .or. size(config_tree, 2) /= max_branch) then
      call fail_born_state('configuration tree has inconsistent bounds')
    end if
    if (size(born_tree, 1) /= 2 .or. size(born_tree, 2) /= max_branch) then
      call fail_born_state('Born tree has inconsistent bounds')
    end if
    if (size(born_momenta, 1) /= 4 .or. &
        size(born_momenta, 2) /= nexternal - 1 .or. &
        any(shape(born_lab_momenta) /= shape(born_momenta))) then
      call fail_born_state('Born momenta have inconsistent bounds')
    end if
    if (size(cbw_mass_state, 1) /= 3 .or. &
        size(cbw_mass_state, 2) /= nexternal .or. &
        any(shape(cbw_width_state) /= shape(cbw_mass_state)) .or. &
        size(cbw_state) /= nexternal .or. &
        size(cbw_level_state) /= nexternal) then
      call fail_born_state('conflicting-BW state has inconsistent bounds')
    end if
    if (size(particle_masses) /= nexternal) then
      call fail_born_state('particle masses have inconsistent bounds')
    end if
    if (size(schannel_masses) /= 2*nexternal + 1) then
      call fail_born_state('s-channel masses have inconsistent bounds')
    end if
  end subroutine validate_bound_born_state

  subroutine require_born_state()
    implicit none

    if (.not. born_state_initialized) then
      call fail_born_state('module state has not been initialized')
    end if
    if (.not. associated(config_mass) .or. &
        .not. associated(config_width) .or. &
        .not. associated(config_forest) .or. &
        .not. associated(config_tree) .or. &
        .not. associated(config_index) .or. &
        .not. associated(born_tree) .or. &
        .not. associated(born_ns_channel) .or. &
        .not. associated(born_nt_channel) .or. &
        .not. associated(born_onebody_index) .or. &
        .not. associated(born_nbranch) .or. &
        .not. associated(born_one_body) .or. &
        .not. associated(born_momenta) .or. &
        .not. associated(born_lab_momenta) .or. &
        .not. associated(cbw_mass_state) .or. &
        .not. associated(cbw_width_state) .or. &
        .not. associated(cbw_level_max_state) .or. &
        .not. associated(cbw_state) .or. &
        .not. associated(cbw_level_state) .or. &
        .not. associated(particle_masses) .or. &
        .not. associated(schannel_masses) .or. &
        .not. allocated(saved_particle_masses) .or. &
        .not. allocated(saved_external_masses)) then
      call fail_born_state('module state is incomplete')
    end if
  end subroutine require_born_state

  subroutine fail_born_state(message)
    implicit none
    character(len=*), intent(in) :: message

    write (*, *) 'genps_born: ', trim(message)
    stop 1
  end subroutine fail_born_state

end module genps_born
