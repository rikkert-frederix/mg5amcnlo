module multiplicative_product_kinematics
  use iso_fortran_env, only: real64
  use ieee_arithmetic, only: ieee_is_finite
  use kin_functions_module, only: dot => dot_impl
  use fks_qcd_splitting, only: AP_reduced, Qterms_reduced_timelike, &
       Qterms_reduced_spacelike
  use fks_model_state_module, only: g => strong_coupling
  use fks_sij_module, only: initialize_fks_sij_module, &
       set_fks_sij_partition_state, fks_sij_impl
  use decay_chain_kinematics, only: boost_from_rest, minkowski_square
  use boostwdir2_module, only: boostwdir2_in_place
  use multiplicative_product, only: product_stage_event
  implicit none
  private

  integer, save :: stored_stages = 0
  integer, save :: stored_local_max = 0
  logical, allocatable, save :: stage_ready(:)
  integer, allocatable, save :: local_count(:), incoming_count(:)
  integer, allocatable, save :: local_i(:), local_j(:), local_ij(:)
  integer, allocatable, save :: local_color(:, :)
  real(real64), allocatable, save :: local_momenta(:, :, :)
  real(real64), allocatable, save :: local_masses(:, :)
  real(real64), allocatable, save :: reference_momentum(:, :)
  real(real64), allocatable, save :: local_sqrt_shat(:)
  real(real64), allocatable, save :: local_xi(:), local_y(:), local_phi(:)

  public :: reset_product_kinematics
  public :: store_product_stage_kinematics
  public :: product_eikonal_factor
  public :: product_collinear_factors
  public :: product_sij_factor
  public :: product_stage_local_momenta
  public :: map_product_final_state
  public :: map_product_initial_state

contains

  subroutine map_product_final_state(stage, emitted, sister, mother, &
       emission_slot, initial_count, group_sizes, group_slots, &
       real_to_born, momenta, &
       jacobian, real_momenta, real_masses, phat, sqrt_shat, pass)
    type(product_stage_event), intent(in) :: stage
    integer, intent(in) :: emitted, sister, mother, emission_slot
    integer, intent(in) :: initial_count
    integer, intent(in) :: group_sizes(:), group_slots(:, :)
    integer, intent(in) :: real_to_born(:)
    real(real64), intent(inout) :: momenta(0:, :)
    real(real64), intent(out) :: jacobian
    real(real64), intent(out) :: real_momenta(0:, :), real_masses(:)
    real(real64), intent(out) :: phat(0:), sqrt_shat
    logical, intent(out) :: pass
    real(real64) :: system_lab(0:3), inverse_system(0:3)
    real(real64) :: system_rest(0:3), born_mother(0:3), recoil(0:3)
    real(real64) :: emitted_p(0:3), sister_p(0:3), mother_p(0:3)
    real(real64) :: work(0:3, size(momenta, 2))
    real(real64) :: recoil_mass2, sister_mass2, sister_mass
    real(real64) :: emitted_energy, sister_length, mother_length
    real(real64) :: cos_i, sin_i, direction(3), denominator
    real(real64) :: sumrec, sumrec2, shy, chy, chymo
    real(real64) :: boost_direction(3), old_sister(0:3)
    real(real64) :: inverse_sister(0:3), child_rest(0:3)
    real(real64) :: cffa, cffb, cffc, root_expression
    real(real64) :: rechat, mjhat, sister_num, sister_den
    real(real64) :: mass_square
    integer :: canonical, born_leg, real_leg, position, solution_sign
    logical :: used(size(momenta, 2))
    real(real64), parameter :: pi = 3.14159265358979323846_real64

    pass = .false.
    jacobian = 0._real64
    real_momenta = 0._real64
    real_masses = 0._real64
    phat = 0._real64
    sqrt_shat = 0._real64
    if (.not. valid_map_shapes(group_sizes, group_slots, real_to_born, &
                               momenta, real_momenta, real_masses, phat)) return
    if (emitted < 1 .or. emitted > size(real_to_born) .or. &
        sister < 1 .or. sister > size(real_to_born) .or. &
        mother < 1 .or. mother > size(group_sizes) .or. &
        emission_slot < 1 .or. emission_slot > size(momenta, 2)) return
    if (stage%xi < -1._real64 .or. abs(stage%xi) >= 1._real64 .or. &
        stage%y < -1._real64 .or. stage%y > 1._real64) return

    system_lab = 0._real64
    if (initial_count < 1 .or. initial_count > 2 .or. &
        initial_count > size(group_sizes)) return
    do born_leg = 1, initial_count
      do position = 1, group_sizes(born_leg)
        system_lab = system_lab + &
             momenta(:, group_slots(position, born_leg))
      end do
    end do
    sqrt_shat = sqrt(max(0._real64, minkowski_square(system_lab)))
    if (sqrt_shat <= 0._real64) return
    system_rest = 0._real64
    system_rest(0) = sqrt_shat
    inverse_system = system_lab
    inverse_system(1:3) = -inverse_system(1:3)
    work = momenta
    used = .false.
    do born_leg = 1, size(group_sizes)
      do position = 1, group_sizes(born_leg)
        canonical = group_slots(position, born_leg)
        if (used(canonical)) cycle
        call boost_from_rest(momenta(:, canonical), inverse_system, &
                             sqrt_shat, work(:, canonical))
        used(canonical) = .true.
      end do
    end do

    born_mother = 0._real64
    do position = 1, group_sizes(mother)
      born_mother = born_mother + work(:, group_slots(position, mother))
    end do
    recoil = system_rest - born_mother
    recoil_mass2 = max(0._real64, minkowski_square(recoil))
    sister_mass2 = max(0._real64, minkowski_square(born_mother))
    sister_mass = sqrt(sister_mass2)
    emitted_energy = abs(stage%xi)*sqrt_shat/2._real64
    solution_sign = merge(1, -1, stage%xi >= 0._real64)

    if (sister_mass <= 1.e-10_real64*sqrt_shat) then
      sister_length = (sqrt_shat**2 - recoil_mass2 - &
           2._real64*sqrt_shat*emitted_energy)/ &
           (2._real64*(sqrt_shat - emitted_energy*(1._real64-stage%y)))
      if (sister_length <= 0._real64) return
    else
      mjhat = sister_mass/sqrt_shat
      rechat = sqrt(recoil_mass2)/sqrt_shat
      cffa = 1._real64 - mjhat**2*(1._real64-stage%y**2)
      cffb = -2._real64*(1._real64-rechat**2-mjhat**2)
      cffc = (1._real64-(rechat-mjhat)**2)* &
             (1._real64-(rechat+mjhat)**2)
      root_expression = abs(stage%xi)**2*cffa + &
                        abs(stage%xi)*cffb + cffc
      if (root_expression < -1.e-10_real64) return
      root_expression = max(0._real64, root_expression)
      sister_num = -abs(stage%xi)*stage%y* &
           (1._real64-rechat**2+mjhat**2-abs(stage%xi)) + &
           (2._real64-abs(stage%xi))*sqrt(root_expression)* &
           real(solution_sign, real64)
      sister_den = (2._real64-abs(stage%xi)*(1._real64-stage%y))* &
                   (2._real64-abs(stage%xi)*(1._real64+stage%y))
      sister_length = sqrt_shat*sister_num/sister_den
      if (sister_length < 0._real64) return
    end if
    mother_length = sqrt(max(0._real64, emitted_energy**2 + &
         sister_length**2 + 2._real64*emitted_energy*sister_length*stage%y))
    if (emitted_energy <= 1.e-14_real64*sqrt_shat) then
      cos_i = stage%y
    else
      cos_i = (mother_length**2-sister_length**2+emitted_energy**2)/ &
              (2._real64*mother_length*emitted_energy)
    end if
    if (abs(cos_i) > 1._real64 + 1.e-8_real64) return
    cos_i = max(-1._real64, min(1._real64, cos_i))
    sin_i = sqrt(max(0._real64, 1._real64-cos_i**2))
    direction = (/sin_i*cos(stage%phi), sin_i*sin(stage%phi), cos_i/)
    emitted_p(0) = emitted_energy
    emitted_p(1:3) = emitted_energy*direction
    sister_p(0) = sqrt(sister_length**2+sister_mass2)
    sister_p(1) = -emitted_p(1)
    sister_p(2) = -emitted_p(2)
    sister_p(3) = mother_length-emitted_p(3)
    phat(0) = sqrt_shat/2._real64
    phat(1:3) = sqrt_shat/2._real64*direction
    call rotate_to_direction(emitted_p, born_mother)
    call rotate_to_direction(sister_p, born_mother)
    call rotate_to_direction(phat, born_mother)

    mother_p = emitted_p+sister_p
    recoil = system_rest-mother_p
    sumrec = recoil(0)+spatial_norm(recoil)
    if (sumrec <= 0._real64 .or. mother_length <= 0._real64) return
    boost_direction = mother_p(1:3)/mother_length
    if (sister_mass <= 1.e-10_real64*sqrt_shat) then
      sumrec2 = sumrec**2
      shy = -(sqrt_shat**2-sumrec2)/(2._real64*sumrec*sqrt_shat)
      chy = (sqrt_shat**2+sumrec2)/(2._real64*sumrec*sqrt_shat)
      chymo = (sqrt_shat-sumrec)**2/(2._real64*sumrec*sqrt_shat)
    else
      if (recoil_mass2 < 1.e-16_real64*sqrt_shat**2) then
        denominator = sqrt_shat*sumrec/(sqrt_shat**2-sister_mass2)
      else
        denominator = sumrec/(2._real64*sqrt_shat*recoil_mass2)* &
             (sqrt_shat**2+recoil_mass2-sister_mass2- &
              sqrt_shat**2*sqrt(max(0._real64, cffc)))
      end if
      if (denominator <= 0._real64) return
      shy = (denominator-1._real64/denominator)/2._real64
      chy = (denominator+1._real64/denominator)/2._real64
      chymo = chy-1._real64
    end if
    do born_leg = initial_count + 1, size(group_sizes)
      if (born_leg == mother) cycle
      do position = 1, group_sizes(born_leg)
        canonical = group_slots(position, born_leg)
        if (shy /= 0._real64) call boostwdir2_in_place(&
             chy, shy, chymo, boost_direction, work(:, canonical))
      end do
    end do

    old_sister = born_mother
    if (group_sizes(mother) == 1) then
      work(:, group_slots(1, mother)) = sister_p
    else
      if (sister_mass <= 0._real64) return
      inverse_sister = old_sister
      inverse_sister(1:3) = -inverse_sister(1:3)
      do position = 1, group_sizes(mother)
        canonical = group_slots(position, mother)
        call boost_from_rest(work(:, canonical), inverse_sister, &
                             sister_mass, child_rest)
        call boost_from_rest(child_rest, sister_p, sister_mass, &
                             work(:, canonical))
      end do
    end if
    work(:, emission_slot) = emitted_p

    used = .false.
    do born_leg = initial_count + 1, size(group_sizes)
      do position = 1, group_sizes(born_leg)
        canonical = group_slots(position, born_leg)
        if (used(canonical)) cycle
        call boost_from_rest(work(:, canonical), system_lab, sqrt_shat, &
                             momenta(:, canonical))
        used(canonical) = .true.
      end do
    end do
    call boost_from_rest(work(:, emission_slot), system_lab, sqrt_shat, &
                         momenta(:, emission_slot))

    do real_leg = 1, size(real_to_born)
      if (real_leg == emitted) then
        real_momenta(:, real_leg) = emitted_p
      else
        born_leg = real_to_born(real_leg)
        if (initial_count == 1 .and. born_leg == 1) then
          ! The incoming decay parent contains the radiation emitted by this
          ! stage.  Its Born descendant group deliberately does not: that
          ! group is also used as the mutable canonical storage below.
          real_momenta(:, real_leg) = system_rest
        else
          do position = 1, group_sizes(born_leg)
            real_momenta(:, real_leg) = real_momenta(:, real_leg) + &
                 work(:, group_slots(position, born_leg))
          end do
        end if
      end if
      mass_square = minkowski_square(real_momenta(:, real_leg))
      if (abs(mass_square) <= 1.e-12_real64*sqrt_shat**2) &
           mass_square = 0._real64
      real_masses(real_leg) = sqrt(max(0._real64, mass_square))
    end do
    if (sister_mass <= 1.e-10_real64*sqrt_shat) then
      denominator = 2._real64-abs(stage%xi)*(1._real64-stage%y)
    else
      denominator = 2._real64-abs(stage%xi)*(1._real64- &
           sister_p(0)/max(sister_length, 1.e-30_real64)*stage%y)
    end if
    if (denominator <= 0._real64 .or. spatial_norm(born_mother) <= 0._real64) return
    jacobian = 2._real64*sqrt_shat**2/(4._real64*pi)**3* &
         sister_length/spatial_norm(born_mother)/denominator
    pass = jacobian >= 0._real64 .and. ieee_is_finite(jacobian) .and. &
           all(ieee_is_finite(momenta))
  end subroutine map_product_final_state


  subroutine map_product_initial_state(stage, emitted, sister, &
       emission_slot, initial_count, group_sizes, group_slots, &
       real_to_born, momenta, jacobian, real_momenta, real_masses, &
       phat, sqrt_shat, pass)
    type(product_stage_event), intent(in) :: stage
    integer, intent(in) :: emitted, sister, emission_slot, initial_count
    integer, intent(in) :: group_sizes(:), group_slots(:, :)
    integer, intent(in) :: real_to_born(:)
    real(real64), intent(inout) :: momenta(0:, :)
    real(real64), intent(out) :: jacobian
    real(real64), intent(out) :: real_momenta(0:, :), real_masses(:)
    real(real64), intent(out) :: phat(0:), sqrt_shat
    logical, intent(out) :: pass
    real(real64) :: born_system(0:3), inverse_born(0:3)
    real(real64) :: real_system(0:3), inverse_real(0:3)
    real(real64) :: work(0:3, size(momenta, 2))
    real(real64) :: p1(0:3), p2(0:3), emitted_p(0:3), phat_tilde(0:3)
    real(real64) :: boosted(0:3)
    real(real64) :: shat_born, xi, y_directed, bstfact
    real(real64) :: shy_t, chy_t, chymo_t, transverse_direction(3)
    real(real64) :: shy_l, chy_l, encmso2, emitted_energy, sin_theta
    real(real64) :: mass_square
    integer :: direction, born_leg, position, canonical, real_leg
    logical :: used(size(momenta, 2))
    real(real64), parameter :: pi = 3.14159265358979323846_real64

    pass = .false.
    jacobian = 0._real64
    real_momenta = 0._real64
    real_masses = 0._real64
    phat = 0._real64
    sqrt_shat = 0._real64
    if (.not. valid_map_shapes(group_sizes, group_slots, real_to_born, &
                               momenta, real_momenta, real_masses, phat)) return
    if (initial_count /= 2 .or. sister < 1 .or. sister > 2 .or. &
        emitted < 1 .or. emitted > size(real_to_born) .or. &
        emission_slot < 1 .or. emission_slot > size(momenta, 2) .or. &
        group_sizes(1) /= 1 .or. group_sizes(2) /= 1) return
    xi = abs(stage%xi)
    if (xi >= 1._real64 .or. stage%y < -1._real64 .or. &
        stage%y > 1._real64) return
    direction = merge(1, -1, sister == 1)
    y_directed = real(direction, real64)*stage%y

    born_system = momenta(:, group_slots(1, 1)) + &
                  momenta(:, group_slots(1, 2))
    shat_born = minkowski_square(born_system)
    if (shat_born <= 0._real64) return
    inverse_born = born_system
    inverse_born(1:3) = -inverse_born(1:3)
    work = momenta
    used = .false.
    do born_leg = 3, size(group_sizes)
      do position = 1, group_sizes(born_leg)
        canonical = group_slots(position, born_leg)
        if (used(canonical)) cycle
        call boost_from_rest(momenta(:, canonical), inverse_born, &
                             sqrt(shat_born), work(:, canonical))
        used(canonical) = .true.
      end do
    end do

    sqrt_shat = sqrt(shat_born/(1._real64-xi))
    bstfact = sqrt(max(0._real64, &
         (2._real64-xi*(1._real64-y_directed))* &
         (2._real64-xi*(1._real64+y_directed))))
    if (bstfact <= 0._real64) return
    shy_t = -xi*sqrt(max(0._real64, 1._real64-y_directed**2))/ &
            (2._real64*sqrt(1._real64-xi))
    chy_t = bstfact/(2._real64*sqrt(1._real64-xi))
    chymo_t = chy_t-1._real64
    transverse_direction = (/-cos(stage%phi), -sin(stage%phi), 0._real64/)
    do born_leg = 3, size(group_sizes)
      do position = 1, group_sizes(born_leg)
        canonical = group_slots(position, born_leg)
        if (shy_t /= 0._real64) call boostwdir2_in_place(&
             chy_t, shy_t, chymo_t, transverse_direction, work(:, canonical))
      end do
    end do

    shy_l = -xi*y_directed/bstfact
    chy_l = (2._real64-xi)/bstfact
    encmso2 = sqrt_shat/2._real64
    p1 = 0._real64
    p2 = 0._real64
    p1(0) = encmso2*(chy_l-shy_l)
    p1(3) = p1(0)
    p2(0) = encmso2*(chy_l+shy_l)
    p2(3) = -p2(0)
    emitted_energy = xi*encmso2
    sin_theta = sqrt(max(0._real64, 1._real64-stage%y**2))
    emitted_p(0) = emitted_energy*(chy_l-shy_l*y_directed)
    emitted_p(1) = emitted_energy*sin_theta*cos(stage%phi)
    emitted_p(2) = emitted_energy*sin_theta*sin(stage%phi)
    emitted_p(3) = emitted_energy*(chy_l*y_directed-shy_l)
    phat_tilde(0) = encmso2*(chy_l-shy_l*y_directed)
    phat_tilde(1) = encmso2*sin_theta*cos(stage%phi)
    phat_tilde(2) = encmso2*sin_theta*sin(stage%phi)
    phat_tilde(3) = encmso2*(chy_l*y_directed-shy_l)

    work(:, group_slots(1, 1)) = p1
    work(:, group_slots(1, 2)) = p2
    work(:, emission_slot) = emitted_p
    do born_leg = 1, size(group_sizes)
      do position = 1, group_sizes(born_leg)
        canonical = group_slots(position, born_leg)
        call boost_from_rest(work(:, canonical), born_system, &
                             sqrt(shat_born), momenta(:, canonical))
      end do
    end do
    call boost_from_rest(work(:, emission_slot), born_system, &
                         sqrt(shat_born), momenta(:, emission_slot))

    do real_leg = 1, size(real_to_born)
      if (real_leg == emitted) then
        real_momenta(:, real_leg) = emitted_p
      else
        born_leg = real_to_born(real_leg)
        do position = 1, group_sizes(born_leg)
          real_momenta(:, real_leg) = real_momenta(:, real_leg) + &
               work(:, group_slots(position, born_leg))
        end do
      end if
    end do
    real_system = p1+p2
    inverse_real = real_system
    inverse_real(1:3) = -inverse_real(1:3)
    do real_leg = 1, size(real_to_born)
      call boost_from_rest(real_momenta(:, real_leg), inverse_real, &
                           sqrt_shat, boosted)
      real_momenta(:, real_leg) = boosted
      mass_square = minkowski_square(real_momenta(:, real_leg))
      if (abs(mass_square) <= 1.e-12_real64*sqrt_shat**2) &
           mass_square = 0._real64
      real_masses(real_leg) = sqrt(max(0._real64, mass_square))
    end do
    call boost_from_rest(phat_tilde, inverse_real, sqrt_shat, phat)
    jacobian = sqrt_shat**2/(4._real64*pi)**3/(1._real64-xi)
    pass = ieee_is_finite(jacobian) .and. jacobian >= 0._real64 .and. &
           all(ieee_is_finite(momenta))
  end subroutine map_product_initial_state


  logical function valid_map_shapes(group_sizes, group_slots, &
       real_to_born, momenta, real_momenta, real_masses, phat)
    integer, intent(in) :: group_sizes(:), group_slots(:, :)
    integer, intent(in) :: real_to_born(:)
    real(real64), intent(in) :: momenta(0:, :)
    real(real64), intent(in) :: real_momenta(0:, :), real_masses(:)
    real(real64), intent(in) :: phat(0:)
    integer :: leg, position

    valid_map_shapes = .false.
    if (size(momenta, 1) /= 4 .or. size(real_momenta, 1) /= 4 .or. &
        size(real_momenta, 2) /= size(real_to_born) .or. &
        size(real_masses) /= size(real_to_born) .or. size(phat) /= 4 .or. &
        size(group_slots, 2) /= size(group_sizes)) return
    do leg = 1, size(group_sizes)
      if (group_sizes(leg) < 1 .or. &
          group_sizes(leg) > size(group_slots, 1)) return
      do position = 1, group_sizes(leg)
        if (group_slots(position, leg) < 1 .or. &
            group_slots(position, leg) > size(momenta, 2)) return
      end do
    end do
    do leg = 1, size(real_to_born)
      if (real_to_born(leg) < 0 .or. &
          real_to_born(leg) > size(group_sizes)) return
    end do
    valid_map_shapes = .true.
  end function valid_map_shapes


  real(real64) function spatial_norm(momentum)
    real(real64), intent(in) :: momentum(0:)
    spatial_norm = sqrt(max(0._real64, sum(momentum(1:3)**2)))
  end function spatial_norm


  subroutine rotate_to_direction(momentum, direction)
    real(real64), intent(inout) :: momentum(0:)
    real(real64), intent(in) :: direction(0:)
    real(real64) :: length, cosine, sine, phi, cphi, sphi
    real(real64) :: old(3)

    length = spatial_norm(direction)
    if (length <= 0._real64) return
    cosine = max(-1._real64, min(1._real64, direction(3)/length))
    sine = sqrt(max(0._real64, 1._real64-cosine**2))
    phi = atan2(direction(2), direction(1))
    cphi = cos(phi)
    sphi = sin(phi)
    old = momentum(1:3)
    momentum(1) = cosine*cphi*old(1)-sphi*old(2)+sine*cphi*old(3)
    momentum(2) = cosine*sphi*old(1)+cphi*old(2)+sine*sphi*old(3)
    momentum(3) = -sine*old(1)+cosine*old(3)
  end subroutine rotate_to_direction

  subroutine reset_product_kinematics(number_of_stages, maximum_local)
    integer, intent(in) :: number_of_stages, maximum_local

    if (number_of_stages < 1 .or. maximum_local < 2) &
         call fail_product_kinematics('invalid storage dimensions')
    if (allocated(stage_ready)) then
      deallocate(stage_ready, local_count, incoming_count, local_i, &
           local_j, local_ij, local_color, local_momenta, local_masses, &
           reference_momentum, local_sqrt_shat, local_xi, local_y, &
           local_phi)
    end if
    stored_stages = number_of_stages
    stored_local_max = maximum_local
    allocate(stage_ready(stored_stages))
    allocate(local_count(stored_stages), incoming_count(stored_stages))
    allocate(local_i(stored_stages), local_j(stored_stages))
    allocate(local_ij(stored_stages))
    allocate(local_color(stored_local_max, stored_stages))
    allocate(local_momenta(0:3, stored_local_max, stored_stages))
    allocate(local_masses(stored_local_max, stored_stages))
    allocate(reference_momentum(0:3, stored_stages))
    allocate(local_sqrt_shat(stored_stages), local_xi(stored_stages))
    allocate(local_y(stored_stages), local_phi(stored_stages))
    stage_ready = .false.
    local_count = 0
    incoming_count = 0
    local_i = 0
    local_j = 0
    local_ij = 0
    local_color = 1
    local_momenta = 0._real64
    local_masses = 0._real64
    reference_momentum = -1._real64
    local_sqrt_shat = 0._real64
    local_xi = 0._real64
    local_y = 0._real64
    local_phi = 0._real64
  end subroutine reset_product_kinematics


  subroutine store_product_stage_kinematics(stage, count, initial_count, &
       i_fks, j_fks, ij_fks, momenta, masses, colors, phat, sqrt_shat, &
       xi, y, phi, pass)
    integer, intent(in) :: stage, count, initial_count
    integer, intent(in) :: i_fks, j_fks, ij_fks
    real(real64), intent(in) :: momenta(0:, :), masses(:), phat(0:)
    integer, intent(in) :: colors(:)
    real(real64), intent(in) :: sqrt_shat, xi, y, phi
    logical, intent(out) :: pass

    pass = .false.
    call check_storage_stage(stage)
    if (count < 2 .or. count > stored_local_max .or. &
        initial_count < 1 .or. initial_count > 2 .or. &
        initial_count > count .or. i_fks < 1 .or. i_fks > count .or. &
        j_fks < 1 .or. j_fks > count .or. ij_fks < 1 .or. &
        ij_fks >= count) return
    if (size(momenta, 1) /= 4 .or. size(momenta, 2) /= count .or. &
        size(masses) /= count .or. size(colors) /= count .or. &
        size(phat) /= 4 .or. sqrt_shat <= 0._real64) return
    if (.not. all(ieee_is_finite(momenta)) .or. &
        .not. all(ieee_is_finite(masses)) .or. &
        .not. all(ieee_is_finite(phat))) return

    local_count(stage) = count
    incoming_count(stage) = initial_count
    local_i(stage) = i_fks
    local_j(stage) = j_fks
    local_ij(stage) = ij_fks
    local_momenta(:, :, stage) = 0._real64
    local_momenta(:, 1:count, stage) = momenta
    local_masses(:, stage) = 0._real64
    local_masses(1:count, stage) = masses
    local_color(:, stage) = 1
    local_color(1:count, stage) = colors
    reference_momentum(:, stage) = phat
    local_sqrt_shat(stage) = sqrt_shat
    local_xi(stage) = xi
    local_y(stage) = y
    local_phi(stage) = phi
    stage_ready(stage) = .true.
    pass = .true.
  end subroutine store_product_stage_kinematics


  subroutine product_stage_local_momenta(stage, momenta, masses, count, pass)
    integer, intent(in) :: stage
    real(real64), intent(out) :: momenta(0:, :), masses(:)
    integer, intent(out) :: count
    logical, intent(out) :: pass

    pass = .false.
    count = 0
    call check_ready_stage(stage)
    count = local_count(stage)
    if (size(momenta, 1) /= 4 .or. size(momenta, 2) < count .or. &
        size(masses) < count) return
    momenta = 0._real64
    masses = 0._real64
    momenta(:, 1:count) = local_momenta(:, 1:count, stage)
    masses(1:count) = local_masses(1:count, stage)
    pass = .true.
  end subroutine product_stage_local_momenta


  subroutine product_eikonal_factor(stage, first, second, value, pass)
    integer, intent(in) :: stage, first, second
    real(real64), intent(out) :: value
    logical, intent(out) :: pass
    real(real64) :: dot_nm, dot_ni, dot_mi, factor
    integer :: sister
    real(real64), parameter :: tiny = 1.e-12_real64

    value = 0._real64
    pass = .false.
    call check_ready_stage(stage)
    if (first < 1 .or. first > local_count(stage) .or. &
        second < 1 .or. second > local_count(stage)) return
    sister = local_j(stage)
    dot_nm = dot(local_momenta(:, second, stage), &
                 local_momenta(:, first, stage))
    if ((first /= sister .and. second /= sister) .or. &
        local_masses(sister, stage) /= 0._real64) then
      dot_mi = dot(local_momenta(:, first, stage), &
                   reference_momentum(:, stage))
      dot_ni = dot(local_momenta(:, second, stage), &
                   reference_momentum(:, stage))
      factor = 1._real64 - local_y(stage)
    else if (first == sister .and. second /= sister) then
      dot_ni = dot(local_momenta(:, second, stage), &
                   reference_momentum(:, stage))
      dot_mi = local_sqrt_shat(stage)/2._real64* &
               local_momenta(0, sister, stage)
      factor = 1._real64
    else if (first /= sister .and. second == sister) then
      dot_ni = local_sqrt_shat(stage)/2._real64* &
               local_momenta(0, sister, stage)
      dot_mi = dot(local_momenta(:, first, stage), &
                   reference_momentum(:, stage))
      factor = 1._real64
    else
      return
    end if
    if (abs(dot_ni*dot_mi) <= tiny) return
    value = dot_nm/(dot_ni*dot_mi)*factor
    pass = ieee_is_finite(value)
    if (.not. pass) value = 0._real64
  end subroutine product_eikonal_factor


  subroutine product_collinear_factors(stage, ap_value, q_value, phase, pass)
    integer, intent(in) :: stage
    real(real64), intent(out) :: ap_value, q_value
    complex(real64), intent(out) :: phase
    logical, intent(out) :: pass
    real(real64) :: ap(2), qterm(2), z, t
    integer :: emitted_color, sister_color, mother_color, direction
    complex(real64), parameter :: imaginary = (0._real64, 1._real64)
    ap_value = 0._real64
    q_value = 0._real64
    phase = (0._real64, 0._real64)
    pass = .false.
    call check_ready_stage(stage)
    emitted_color = local_color(local_i(stage), stage)
    sister_color = local_color(local_j(stage), stage)
    call mother_color_representation(emitted_color, sister_color, &
         local_j(stage) <= incoming_count(stage), mother_color, pass)
    if (.not. pass) return
    if (local_j(stage) > incoming_count(stage)) then
      if (local_momenta(0, local_i(stage), stage) + &
          local_momenta(0, local_j(stage), stage) <= 0._real64) return
      z = 1._real64 - local_momenta(0, local_i(stage), stage)/ &
          (local_momenta(0, local_i(stage), stage) + &
           local_momenta(0, local_j(stage), stage))
      t = z*local_sqrt_shat(stage)**2/4._real64
      if (t <= 0._real64) return
      call AP_reduced(sister_color, emitted_color, t, z, g, ap)
      call Qterms_reduced_timelike(&
           sister_color, emitted_color, t, z, g, qterm)
      if (abs(sister_color) == 3 .and. emitted_color == 8) &
           qterm(1) = 0._real64
      phase = exp(2._real64*imaginary*local_phi(stage))
    else
      z = 1._real64 - local_xi(stage)
      t = z*local_sqrt_shat(stage)**2/4._real64
      if (t <= 0._real64) return
      call AP_reduced(mother_color, emitted_color, t, z, g, ap)
      call Qterms_reduced_spacelike(&
           mother_color, emitted_color, t, z, g, qterm)
      if (local_j(stage) == 1) then
        direction = 1
      else
        direction = -1
      end if
      phase = exp(-2._real64*real(direction, real64)*imaginary* &
                  local_phi(stage))
    end if
    ap_value = ap(1)
    q_value = qterm(1)
    pass = ieee_is_finite(ap_value) .and. ieee_is_finite(q_value) .and. &
           ieee_is_finite(real(phase)) .and. ieee_is_finite(aimag(phase))
  end subroutine product_collinear_factors


  subroutine product_sij_factor(stage, partners, particle_types, is_aorg, &
       value, pass)
    integer, intent(in) :: stage
    integer, intent(in) :: partners(:, 0:), particle_types(:)
    logical, intent(in) :: is_aorg(:)
    real(real64), intent(out) :: value
    logical, intent(out) :: pass
    real(real64) :: counterevent_momenta(0:3, 0:2)

    value = 0._real64
    pass = .false.
    call check_ready_stage(stage)
    if (size(partners, 1) /= local_count(stage) .or. &
        size(partners, 2) /= local_count(stage) + 1 .or. &
        size(particle_types) /= local_count(stage) .or. &
        size(is_aorg) /= local_count(stage)) return
    call initialize_fks_sij_module(local_count(stage), &
         incoming_count(stage), 1.5_real64, 1.5_real64, 1._real64, &
         1.e-2_real64, .true.)
    counterevent_momenta = -1._real64
    counterevent_momenta(:, 0) = reference_momentum(:, stage)
    counterevent_momenta(:, 2) = reference_momentum(:, stage)
    call set_fks_sij_partition_state(partners, particle_types, is_aorg, &
         local_i(stage), local_j(stage), 0._real64, &
         local_sqrt_shat(stage), local_sqrt_shat(stage)**2, &
         counterevent_momenta, local_masses(1:local_count(stage), stage))
    value = fks_sij_impl(&
         local_momenta(:, 1:local_count(stage), stage), local_i(stage), &
         local_j(stage), local_xi(stage), local_y(stage))
    pass = ieee_is_finite(value) .and. value >= 0._real64
    if (.not. pass) value = 0._real64
  end subroutine product_sij_factor


  subroutine mother_color_representation(emitted, sister, initial_sister, &
       mother, pass)
    integer, intent(in) :: emitted, sister
    logical, intent(in) :: initial_sister
    integer, intent(out) :: mother
    logical, intent(out) :: pass

    mother = 0
    pass = .true.
    if (abs(emitted) == abs(sister) .and. abs(emitted) > 1) then
      mother = 8
      if ((initial_sister .and. abs(emitted) == 3 .and. &
           sister /= emitted) .or. &
          (.not. initial_sister .and. abs(emitted) == 3 .and. &
           sister /= -emitted)) pass = .false.
    else if (abs(emitted) == 3 .and. sister == 8 .and. initial_sister) then
      mother = -emitted
    else if (emitted == 8 .and. abs(sister) == 3) then
      mother = sister
    else
      pass = .false.
    end if
  end subroutine mother_color_representation


  subroutine check_storage_stage(stage)
    integer, intent(in) :: stage
    if (.not. allocated(stage_ready)) &
         call fail_product_kinematics('stage storage is not initialized')
    if (stage < 1 .or. stage > stored_stages) &
         call fail_product_kinematics('stage index is outside storage')
  end subroutine check_storage_stage


  subroutine check_ready_stage(stage)
    integer, intent(in) :: stage
    call check_storage_stage(stage)
    if (.not. stage_ready(stage)) &
         call fail_product_kinematics('stage kinematics are unavailable')
  end subroutine check_ready_stage


  subroutine fail_product_kinematics(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in multiplicative product kinematics: '// &
         trim(message)
    stop 1
  end subroutine fail_product_kinematics

end module multiplicative_product_kinematics
