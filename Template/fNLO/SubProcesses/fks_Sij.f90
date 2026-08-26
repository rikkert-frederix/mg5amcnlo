module fks_sij_module
  use process_dimensions, only: process_dimensions_initialized, &
                                process_nexternal => nexternal, process_nincoming => nincoming, &
                                validate_process_dimensions
  use kin_functions_module, only: dot => dot_impl
  implicit none
  private

  integer :: particle_count = 0
  integer :: incoming_count = 0
  double precision :: fks_energy_power = 0d0
  double precision :: fks_angle_power = 0d0
  double precision :: h_damp_power = 0d0
  double precision :: h_damp_edge = 0d0
  logical :: module_initialized = .false.
  logical :: partition_state_ready = .false.
  logical :: kinematic_state_ready = .false.

  integer, allocatable :: fks_j_from_i_state(:, :)
  integer, allocatable :: particle_type_state(:)
  logical, allocatable :: is_aorg_state(:)
  double precision, allocatable :: particle_mass_state(:)
  double precision, allocatable :: p_i_fks_cnt_state(:, :)
  integer, allocatable :: ijskip(:, :)

  integer :: i_fks_state = 0
  integer :: j_fks_state = 0
  double precision :: ybst_til_tocm_state = 0d0
  double precision :: sqrtshat_state = 0d0
  double precision :: shat_state = 0d0

  public :: initialize_fks_sij_module
  public :: set_fks_sij_partition_state
  public :: fks_sij_impl

contains

  subroutine initialize_fks_sij_module(nexternal_in, nincoming_in, &
                                       fks_a_in, fks_b_in, a_h_damp_in, &
                                       one_h_damp_in)
    implicit none
    integer, intent(in) :: nexternal_in, nincoming_in
    double precision, intent(in) :: fks_a_in, fks_b_in
    double precision, intent(in) :: a_h_damp_in, one_h_damp_in

    if (nexternal_in < 1) then
      call fail_validation('NEXTERNAL must be positive')
    end if
    if (nincoming_in < 1 .or. nincoming_in > 2 .or. &
        nincoming_in > nexternal_in) then
      call fail_validation('NINCOMING is outside the supported range')
    end if
    if (fks_a_in <= 0d0 .or. fks_b_in <= 0d0) then
      call fail_validation('the FKS partition powers must be positive')
    end if
    if (a_h_damp_in <= 0d0) then
      call fail_validation('the h-damping power must be positive')
    end if
    if (one_h_damp_in < 0d0 .or. one_h_damp_in >= 0.5d0) then
      call fail_validation('the h-damping edge is outside [0,1/2)')
    end if

    if (process_dimensions_initialized) then
      call validate_process_dimensions()
      if (process_nexternal /= nexternal_in .or. &
          process_nincoming /= nincoming_in) then
        call fail_validation('generated and process dimensions disagree')
      end if
    end if

    if (module_initialized) then
      if (particle_count /= nexternal_in .or. &
          incoming_count /= nincoming_in .or. &
          fks_energy_power /= fks_a_in .or. &
          fks_angle_power /= fks_b_in .or. &
          h_damp_power /= a_h_damp_in .or. &
          h_damp_edge /= one_h_damp_in) then
        call fail_validation('module was reinitialized with new values')
      end if
      return
    end if

    particle_count = nexternal_in
    incoming_count = nincoming_in
    fks_energy_power = fks_a_in
    fks_angle_power = fks_b_in
    h_damp_power = a_h_damp_in
    h_damp_edge = one_h_damp_in

    allocate (fks_j_from_i_state(particle_count, 0:particle_count))
    allocate (particle_type_state(particle_count))
    allocate (is_aorg_state(particle_count))
    allocate (particle_mass_state(particle_count))
    allocate (p_i_fks_cnt_state(0:3, 0:2))
    allocate (ijskip(particle_count, particle_count))

    fks_j_from_i_state = 0
    particle_type_state = 0
    is_aorg_state = .false.
    particle_mass_state = 0d0
    p_i_fks_cnt_state = -1d0
    ijskip = 0
    module_initialized = .true.
  end subroutine initialize_fks_sij_module

  subroutine set_fks_sij_partition_state(fks_j_from_i_in, &
                                         particle_type_in, is_aorg_in, i_fks_in, j_fks_in, &
                                         ybst_til_tocm_in, sqrtshat_in, shat_in, &
                                         p_i_fks_cnt_in, particle_mass_in)
    implicit none
    integer, intent(in) :: fks_j_from_i_in(:, 0:)
    integer, intent(in) :: particle_type_in(:)
    logical, intent(in) :: is_aorg_in(:)
    integer, intent(in) :: i_fks_in, j_fks_in
    double precision, intent(in) :: ybst_til_tocm_in
    double precision, intent(in) :: sqrtshat_in, shat_in
    double precision, intent(in) :: p_i_fks_cnt_in(0:, 0:)
    double precision, intent(in) :: particle_mass_in(:)

    call require_initialized()
    if (size(fks_j_from_i_in, 1) /= particle_count .or. &
        size(fks_j_from_i_in, 2) /= particle_count + 1) then
      call fail_validation('FKS_J_FROM_I has the wrong shape')
    end if
    if (size(particle_type_in) /= particle_count .or. &
        size(is_aorg_in) /= particle_count) then
      call fail_validation('the chooser arrays have the wrong shape')
    end if
    if (i_fks_in < 1 .or. i_fks_in > particle_count .or. &
        j_fks_in < 1 .or. j_fks_in > particle_count) then
      call fail_validation('the FKS indices are outside the process')
    end if
    if (size(p_i_fks_cnt_in, 1) /= 4 .or. &
        size(p_i_fks_cnt_in, 2) /= 3) then
      call fail_validation('counterevent momentum storage has wrong shape')
    end if

    fks_j_from_i_state = fks_j_from_i_in
    particle_type_state = particle_type_in
    is_aorg_state = is_aorg_in
    i_fks_state = i_fks_in
    j_fks_state = j_fks_in
    p_i_fks_cnt_state = p_i_fks_cnt_in
    call set_fks_sij_kinematic_state(ybst_til_tocm_in, sqrtshat_in, &
                                     shat_in, particle_mass_in)
    partition_state_ready = .true.
  end subroutine set_fks_sij_partition_state

  subroutine set_fks_sij_kinematic_state(ybst_til_tocm_in, sqrtshat_in, &
                                         shat_in, particle_mass_in)
    implicit none
    double precision, intent(in) :: ybst_til_tocm_in
    double precision, intent(in) :: sqrtshat_in, shat_in
    double precision, intent(in) :: particle_mass_in(:)

    call require_initialized()
    if (size(particle_mass_in) /= particle_count) then
      call fail_validation('the particle-mass array has the wrong shape')
    end if
    ybst_til_tocm_state = ybst_til_tocm_in
    sqrtshat_state = sqrtshat_in
    shat_state = shat_in
    particle_mass_state = particle_mass_in
    kinematic_state_ready = .true.
  end subroutine set_fks_sij_kinematic_state

  double precision function fks_sij_impl(p, ii_fks, jj_fks, &
                                         xi_i_fks, y_ij_fks)
    implicit none
    double precision, intent(in) :: p(0:, :)
    integer, intent(in) :: ii_fks, jj_fks
    double precision, intent(in) :: xi_i_fks, y_ij_fks

    double precision, parameter :: tiny = 1d-6
    double precision, parameter :: zero = 0d0
    double precision, parameter :: one = 1d0
    integer, parameter :: ione = 1
    integer, parameter :: itwo = 2
    double precision :: shattmp
    double precision :: hfact, z
    double precision :: phat_ii(0:3), e_ii_resc, xnum, xden
    double precision :: inverse_sij
    integer :: i, j, k, l, kk, ll, ihdamp, isorsc
    logical :: setsijzero

    call require_partition_state()
    call check_momenta_shape(p)
    if (ii_fks < 1 .or. ii_fks > particle_count .or. &
        jj_fks < 1 .or. jj_fks > particle_count) then
      call fail_validation('requested FKS indices are outside the process')
    end if

    if (p(0, 1) <= 0d0) then
      fks_sij_impl = 0d0
      return
    end if

    if (is_aorg_state(jj_fks) .and. .not. is_aorg_state(ii_fks) .and. &
        jj_fks > incoming_count) then
      write (*, *) 'Error #0 in fks_Sij', ii_fks, jj_fks, &
        is_aorg_state(ii_fks), is_aorg_state(jj_fks)
      stop
    end if

    if (incoming_count == 2) then
      shattmp = 2d0*dot(p(0:3, 1), p(0:3, 2))
    else
      shattmp = p(0, 1)**2
    end if
    if (abs(shattmp/shat_state - 1d0) > 1d-5) then
      write (*, *) 'Error in fks_Sij: inconsistent shat #1'
      write (*, *) shattmp, shat_state
      stop
    end if

    setsijzero = .false.
    e_ii_resc = -one

    if (xi_i_fks > tiny .or. ii_fks /= i_fks_state) then
      do i = 0, 3
        phat_ii(i) = p(i, ii_fks)
      end do
      e_ii_resc = one
    else if (xi_i_fks <= tiny .and. ii_fks == i_fks_state) then
      isorsc = 0
      if (1d0 - y_ij_fks < tiny .and. jj_fks == j_fks_state .and. &
          particle_mass_state(j_fks_state) == 0d0) isorsc = 2
      if (p_i_fks_cnt_state(0, isorsc) < 0d0) then
        if (xi_i_fks == 0d0) then
          write (*, *) 'Error #7 in fks_Sij', isorsc, xi_i_fks, y_ij_fks
          stop
        end if
        if (p(0, ii_fks) /= 0d0) then
          write (*, *) 'WARNING in fks_Sij: no cnt momenta', &
            isorsc, xi_i_fks, y_ij_fks
          do i = 0, 3
            phat_ii(i) = p(i, ii_fks)
          end do
          e_ii_resc = one
        else
          write (*, *) 'Error #8 in fks_Sij', isorsc, xi_i_fks, y_ij_fks
          stop
        end if
      else
        do i = 0, 3
          phat_ii(i) = p_i_fks_cnt_state(i, isorsc)
        end do
        e_ii_resc = xi_i_fks
      end if
    else
      write (*, *) 'fks_Sij: do not know what to do', &
        ii_fks, i_fks_state, xi_i_fks
    end if

    ! The legacy FIRSTTIME flag was intentionally never reset, so this
    ! configuration-dependent table must be rebuilt on every invocation.
    ijskip = 0
    do i = 1, particle_count
      if (fks_j_from_i_state(i, 0) /= 0) then
        do j = 1, fks_j_from_i_state(i, 0)
          kk = i
          ll = fks_j_from_i_state(i, j)
          if (incoming_count /= 2 .and. ll <= incoming_count) cycle
          if (ijskip(kk, ll) == 0 .and. ijskip(ll, kk) == 0) then
            ijskip(kk, ll) = 1
          else if (ijskip(kk, ll) == 0 .and. ijskip(ll, kk) == 1) then
            ijskip(kk, ll) = 2
            if (.not. is_aorg_state(kk) .or. &
                .not. is_aorg_state(ll)) then
              write (*, *) 'Error #1 in fks_Sij', kk, ll, &
                is_aorg_state(kk), is_aorg_state(ll)
              do k = 1, particle_count
                write (*, *) k, (ijskip(k, l), l=1, particle_count)
              end do
              stop
            end if
          else
            write (*, *) 'Error #2 in fks_Sij', kk, ll
            stop
          end if
        end do
      end if
    end do

    inverse_sij = 0d0
    ihdamp = 0
    hfact = 1d0

    do i = 1, particle_count
      if (fks_j_from_i_state(i, 0) /= 0) then
        do j = 1, fks_j_from_i_state(i, 0)
          kk = i
          ll = fks_j_from_i_state(i, j)
          if (incoming_count /= 2 .and. ll <= incoming_count) cycle
          if (ijskip(kk, ll) /= 1) cycle
          if (is_aorg_state(ll) .and. .not. is_aorg_state(kk) .and. &
              ll > incoming_count) then
            write (*, *) 'Error #3 in fks_Sij', kk, ll, &
              is_aorg_state(kk), is_aorg_state(ll)
            stop
          end if
          if ((.not. is_aorg_state(ll) .and. &
               .not. is_aorg_state(kk) .and. &
               particle_mass_state(ll) /= zero) .or. &
              particle_mass_state(kk) /= zero) then
            write (*, *) 'Error #4 in fks_Sij', kk, ll, &
              is_aorg_state(kk), is_aorg_state(ll), &
              particle_mass_state(kk), particle_mass_state(ll)
            stop
          end if
          if ((kk == ii_fks .and. ll == jj_fks) .or. &
              (ll == ii_fks .and. kk == jj_fks)) then
            inverse_sij = inverse_sij + 1d0
            if (is_aorg_state(ll) .and. is_aorg_state(kk) .and. &
                jj_fks > incoming_count) then
              z = p(0, ii_fks)/(p(0, ii_fks) + p(0, jj_fks))
              hfact = hfact*h_damp_impl(z)
              ihdamp = ihdamp + 1
            end if
          else if (kk == ii_fks .and. ll /= jj_fks) then
            call get_dkl_sij_impl(phat_ii, p(0:3, jj_fks), &
                                  p(0:3, ione), p(0:3, itwo), one, ii_fks, &
                                  particle_type_state(ii_fks), jj_fks, &
                                  particle_type_state(jj_fks), ione, &
                                  xnum, setsijzero)
            call get_dkl_sij_impl(phat_ii, p(0:3, ll), &
                                  p(0:3, ione), p(0:3, itwo), one, ii_fks, &
                                  particle_type_state(ii_fks), ll, &
                                  particle_type_state(ll), ione, &
                                  xden, setsijzero)
            if (setsijzero) then
              fks_sij_impl = 0d0
              return
            end if
            inverse_sij = inverse_sij + xnum/xden
          else if (ll == ii_fks .and. kk /= jj_fks) then
            call get_dkl_sij_impl(phat_ii, p(0:3, jj_fks), &
                                  p(0:3, ione), p(0:3, itwo), one, ii_fks, &
                                  particle_type_state(ii_fks), jj_fks, &
                                  particle_type_state(jj_fks), ione, &
                                  xnum, setsijzero)
            call get_dkl_sij_impl(phat_ii, p(0:3, kk), &
                                  p(0:3, ione), p(0:3, itwo), one, ii_fks, &
                                  particle_type_state(ii_fks), kk, &
                                  particle_type_state(kk), ione, &
                                  xden, setsijzero)
            if (setsijzero) then
              fks_sij_impl = 0d0
              return
            end if
            inverse_sij = inverse_sij + xnum/xden
          else
            call get_dkl_sij_impl(phat_ii, p(0:3, jj_fks), &
                                  p(0:3, ione), p(0:3, itwo), e_ii_resc, ii_fks, &
                                  particle_type_state(ii_fks), jj_fks, &
                                  particle_type_state(jj_fks), itwo, &
                                  xnum, setsijzero)
            call get_dkl_sij_impl(p(0:3, kk), p(0:3, ll), &
                                  p(0:3, ione), p(0:3, itwo), one, kk, &
                                  particle_type_state(kk), ll, particle_type_state(ll), &
                                  itwo, xden, setsijzero)
            if (setsijzero) then
              fks_sij_impl = 0d0
              return
            end if
            inverse_sij = inverse_sij + xnum/xden
          end if
        end do
      end if
    end do

    if (ihdamp /= 0 .and. ihdamp /= 1) then
      write (*, *) 'Error #5 in fks_Sij', ihdamp
      stop
    end if

    fks_sij_impl = 1d0/inverse_sij*hfact
    if (fks_sij_impl < 0d0 .or. fks_sij_impl > 1d0) then
      write (*, *) 'Error #6 in fks_Sij', fks_sij_impl
      stop
    end if
  end function fks_sij_impl

  subroutine get_dkl_sij_impl(p1, p2, ka, kb, e1resc, i1, itype1, &
                              i2, itype2, ioneortwo, dkl_sij, setsijzero)
    implicit none
    double precision, intent(in) :: p1(0:3), p2(0:3)
    double precision, intent(in) :: ka(0:3), kb(0:3), e1resc
    integer, intent(in) :: i1, itype1, i2, itype2
    integer, intent(in) :: ioneortwo
    double precision, intent(out) :: dkl_sij
    logical, intent(out) :: setsijzero

    double precision, parameter :: vtiny = 1d-8
    double precision, parameter :: zero = 0d0
    double precision :: energy, e1, e2, beta, beta1, beta2
    double precision :: angle, costhfks

    call require_kinematic_state()
    call check_particle_index(i1)
    call check_particle_index(i2)

    if (e1resc < 0d0) then
      write (*, *) 'Error #0 in dkl_Sij', e1resc, i1, i2
    end if

    setsijzero = .false.
    dkl_sij = 0d0
    energy = 1d0
    e1 = -1d0
    e2 = -1d0
    if (ioneortwo == 2) then
      e1 = get_cms_energy_impl(p1, ka, kb, ybst_til_tocm_state, &
                               shat_state)
      energy = energy*e1*e1resc/(sqrtshat_state/2d0)
      setsijzero = setsijzero .or. &
                   (e1*e1resc) < (vtiny*sqrtshat_state/2d0)
    else if (ioneortwo /= 1) then
      write (*, *) 'Error in dkl_Sij: unknown option', ioneortwo
      stop
    end if
    if (setsijzero) return

    e2 = get_cms_energy_impl(p2, ka, kb, ybst_til_tocm_state, shat_state)
    energy = energy*e2/(sqrtshat_state/2d0)
    setsijzero = setsijzero .or. e2 < (vtiny*sqrtshat_state/2d0)
    if (setsijzero) return

    if (energy > 0d0) then
      energy = energy**fks_energy_power
    else if (energy < 0d0) then
      write (*, *) 'Error #3 in dkl_Sij', energy
      stop
    end if

    beta = 1d0
    call get_cms_costh_fks_impl(p1, p2, ka, kb, e1, e2, &
                                particle_mass_state(i1), particle_mass_state(i2), beta1, &
                                beta2, costhfks, ybst_til_tocm_state, shat_state)
    if (itype1 /= 8 .and. particle_mass_state(i1) /= zero) &
      beta = beta*beta1
    if (itype2 /= 8 .and. particle_mass_state(i2) /= zero) &
      beta = beta*beta2
    angle = 1d0 - beta*costhfks
    setsijzero = setsijzero .or. angle < vtiny
    if (angle > 0d0) then
      angle = angle**fks_angle_power
    else if (angle < 0d0) then
      write (*, *) 'Error #4 in dkl_Sij', angle
      stop
    end if

    dkl_sij = energy*angle
    if (dkl_sij == 0d0 .and. .not. setsijzero) then
      write (*, *) 'Error #5 in dkl_Sij'
      stop
    end if
  end subroutine get_dkl_sij_impl

  double precision function get_cms_energy_impl(p, ka, kb, &
                                                ybst_til_tocm, shat)
    implicit none
    double precision, intent(in) :: p(0:3), ka(0:3), kb(0:3)
    double precision, intent(in) :: ybst_til_tocm, shat
    double precision :: xden, xnum

    if (ybst_til_tocm == 0d0) then
      get_cms_energy_impl = p(0)
    else
      xden = dot(p, ka) + dot(p, kb)
      xnum = 2*dot(ka, kb)
      if (abs(xnum/shat - 1d0) > 1d-6) then
        write (*, *) 'Inconsistency in get_cms_energy'
        stop
      end if
      get_cms_energy_impl = xden/sqrt(xnum)
    end if
  end function get_cms_energy_impl

  subroutine get_cms_costh_fks_impl(p1, p2, ka, kb, e1, e2, xm1, xm2, &
                                    beta1, beta2, costhfks, ybst_til_tocm, shat)
    implicit none
    double precision, intent(in) :: p1(0:3), p2(0:3)
    double precision, intent(in) :: ka(0:3), kb(0:3)
    double precision, intent(inout) :: e1, e2
    double precision, intent(in) :: xm1, xm2
    double precision, intent(out) :: beta1, beta2, costhfks
    double precision, intent(in) :: ybst_til_tocm, shat
    double precision, parameter :: tiny = 1d-6
    double precision :: tmp

    if (ybst_til_tocm == 0d0) then
      tmp = costh_fks_impl(p1, p2)
      beta1 = sqrt(1d0 - (xm1/p1(0))**2)
      beta2 = sqrt(1d0 - (xm2/p2(0))**2)
    else
      if (e1 < 0d0) then
        e1 = get_cms_energy_impl(p1, ka, kb, ybst_til_tocm, shat)
      end if
      beta1 = sqrt(1d0 - (xm1/e1)**2)
      if (e2 < 0d0) then
        e2 = get_cms_energy_impl(p2, ka, kb, ybst_til_tocm, shat)
      end if
      beta2 = sqrt(1d0 - (xm2/e2)**2)
      tmp = (1d0 - dot(p1, p2)/(e1*e2))/(beta1*beta2)
      if ((abs(tmp) - 1d0) > tiny) then
        write (*, *) 'Warning in get_cms_costh_fks', tmp
        tmp = sign(1d0, tmp)
      else if ((abs(tmp) - 1d0) <= tiny .and. &
               (abs(tmp) - 1d0) >= 0d0) then
        tmp = sign(1d0, tmp)
      end if
    end if
    costhfks = tmp
  end subroutine get_cms_costh_fks_impl

  double precision function costh_fks_impl(p1, p2)
    implicit none
    double precision, intent(in) :: p1(0:3), p2(0:3)
    double precision, parameter :: tiny = 1d-6
    double precision :: length1, length2

    length1 = sqrt(p1(1)**2 + p1(2)**2 + p1(3)**2)
    length2 = sqrt(p2(1)**2 + p2(2)**2 + p2(3)**2)
    if (length1 /= 0d0 .and. length2 /= 0d0) then
      costh_fks_impl = (p1(1)*p2(1) + p1(2)*p2(2) + &
                        p1(3)*p2(3))/length1/length2
      if ((abs(costh_fks_impl) - 1d0) > tiny) then
        write (*, *) 'Error in costh_fks', costh_fks_impl
        stop
      else if ((abs(costh_fks_impl) - 1d0) <= tiny .and. &
               (abs(costh_fks_impl) - 1d0) >= 0d0) then
        costh_fks_impl = sign(1d0, costh_fks_impl)
      end if
    else
      costh_fks_impl = 1d0
    end if
  end function costh_fks_impl

  double precision function h_damp_impl(x)
    implicit none
    double precision, intent(in) :: x
    double precision :: y

    call require_initialized()
    if (x < 0d0 .or. x > 1d0) then
      write (*, *) 'ERROR in h_damp', x
      stop
    end if

    if (x <= h_damp_edge) then
      h_damp_impl = 1d0
    else if (x >= 1d0 - h_damp_edge) then
      h_damp_impl = 0d0
    else
      y = (x - h_damp_edge)/(1d0 - 2d0*h_damp_edge)
      if (y < 0d0 .or. y > 1d0) then
        write (*, *) 'ERROR in h_damp', x, y, h_damp_edge
        stop
      end if
      h_damp_impl = (1d0 - y)**(2d0*h_damp_power)/ &
                    ((1d0 - y)**(2d0*h_damp_power) + y**(2d0*h_damp_power))
    end if
  end function h_damp_impl

  subroutine check_momenta_shape(p)
    implicit none
    double precision, intent(in) :: p(0:, :)

    if (size(p, 1) /= 4 .or. size(p, 2) /= particle_count) then
      call fail_validation('the momentum array has the wrong shape')
    end if
  end subroutine check_momenta_shape

  subroutine check_particle_index(index)
    implicit none
    integer, intent(in) :: index

    if (index < 1 .or. index > particle_count) then
      call fail_validation('a particle index is outside the process')
    end if
  end subroutine check_particle_index

  subroutine require_initialized()
    implicit none

    if (.not. module_initialized) then
      call fail_validation('the module has not been initialized')
    end if
  end subroutine require_initialized

  subroutine require_partition_state()
    implicit none

    call require_initialized()
    if (.not. partition_state_ready) then
      call fail_validation('the FKS partition state has not been set')
    end if
    call require_kinematic_state()
  end subroutine require_partition_state

  subroutine require_kinematic_state()
    implicit none

    call require_initialized()
    if (.not. kinematic_state_ready) then
      call fail_validation('the FKS kinematic state has not been set')
    end if
  end subroutine require_kinematic_state

  subroutine fail_validation(message)
    implicit none
    character(len=*), intent(in) :: message

    write (*, *) 'fks_sij_module: ', trim(message)
    stop 1
  end subroutine fail_validation

end module fks_sij_module
