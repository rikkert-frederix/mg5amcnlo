module nlo_decay_kinematics
  use process_dimensions, only: nexternal, nincoming, &
                                validate_process_dimensions
  use phase_space_kinematics, only: phase_space_lambda, rotate_invar
  use boostwdir2_module, only: boostwdir2_in_place
  use decay_chain_kinematics, only: generate_nbody, boost_from_rest, &
                                    minkowski_square
  use decay_chain_parameters, only: decay_dummy_width_ratio, &
                                    decay_physical_width
  use nlo_decay_metadata, only: initialize_nlo_decay_metadata, &
       has_nlo_decay, nlo_decay_metadata_revision, corrected_parent_pdg, &
       nlo_decay_born_context, nlo_decay_context_for_fks, &
       nlo_decay_production_count, nlo_decay_production_pdg, &
       nlo_decay_production_target_kind, &
       nlo_decay_production_target_id, nlo_decay_local_count, &
       nlo_decay_local_pdg, nlo_decay_local_is_final, &
       nlo_decay_local_target_kind, nlo_decay_local_target_id, &
       nlo_decay_fks_i, nlo_decay_fks_j, nlo_decay_fks_ij, &
       nlo_decay_real_to_born, nlo_decay_leg_target, &
       nlo_decay_node_target, nlo_decay_corrected_node, &
       nlo_decay_node_count, nlo_decay_leaf_count, &
       nlo_decay_node_pdg, nlo_decay_node_child_count, &
       nlo_decay_node_child_kind, nlo_decay_node_child_id, &
       nlo_decay_leaf_pdg, nlo_decay_leaf_visible, &
       nlo_decay_leaf_child, nlo_decay_node_child
  use fnlo_process_common, only: soft_counterevent, collinear_counterevent, &
       soft_collinear_counterevent, real_event, softtest, colltest, &
       xi_i_fks_fix, y_ij_fks_fix, xij_aor, &
       event_from_decay => from_decay
  implicit none
  private

  logical, save :: initialized = .false.
  integer, save :: loaded_metadata_revision = -1
  double precision, save :: parent_mass = 0d0
  double precision, save :: massive_xjac_cache = 1d0
  double precision, allocatable, save :: production_masses(:)
  double precision, allocatable, save :: born_local_masses(:)
  double precision, allocatable, save :: node_masses(:)
  double precision, allocatable, save :: leaf_masses(:)
  double precision, allocatable, save :: production_born(:, :)
  double precision, allocatable, save :: born_local(:, :)
  integer, allocatable, save :: node_random_start(:)
  double precision, allocatable, save :: local_event_cache(:, :, :)
  double precision, allocatable, save :: local_fks_momentum_cache(:, :)
  integer, allocatable, save :: local_event_configuration(:)
  logical, allocatable, save :: local_event_valid(:)
  double precision, save :: parent_born(0:3) = 0d0

  public :: initialize_nlo_decay_kinematics
  public :: nlo_decay_minimum_production_mass
  public :: nlo_decay_production_mass
  public :: get_nlo_decay_production_momenta
  public :: fill_nlo_decay_born_masses
  public :: fill_nlo_decay_event_masses
  public :: generate_nlo_decay_born_momenta
  public :: generate_nlo_decay_fks_event
  public :: nlo_decay_fks_sister_mass
  public :: get_nlo_decay_event_momenta
  public :: get_nlo_decay_born_kernel
  public :: get_nlo_decay_counterevent_fks_momenta
  public :: get_nlo_decay_mass_buffer
  public :: nlo_decay_parent_mass
  public :: set_nlo_decay_cut_mask

  interface
    double precision function get_mass_from_id(id)
      integer, intent(in) :: id
    end function get_mass_from_id
  end interface

contains

  subroutine initialize_nlo_decay_kinematics()
    integer :: context, leg, node, target, current_revision

    call initialize_nlo_decay_metadata()
    current_revision = nlo_decay_metadata_revision()
    if (initialized .and. loaded_metadata_revision == current_revision) return
    call clear_nlo_decay_kinematics()
    initialized = .false.
    loaded_metadata_revision = current_revision
    call validate_process_dimensions()
    if (.not. has_nlo_decay()) then
      initialized = .true.
      return
    end if

    allocate(production_masses(nlo_decay_production_count()))
    allocate(born_local_masses(nexternal))
    allocate(node_masses(nlo_decay_node_count()))
    allocate(leaf_masses(nlo_decay_leaf_count()))
    allocate(production_born(0:3, nexternal))
    allocate(born_local(0:3, nexternal))
    allocate(node_random_start(nlo_decay_node_count()))
    allocate(local_event_cache(0:3, nexternal, 0:real_event))
    allocate(local_fks_momentum_cache(0:3, 0:real_event))
    allocate(local_event_configuration(0:real_event))
    allocate(local_event_valid(0:real_event))
    production_masses = 0d0
    born_local_masses = 0d0
    node_masses = 0d0
    leaf_masses = 0d0
    production_born = 0d0
    born_local = 0d0
    node_random_start = 0
    local_event_cache = 0d0
    local_fks_momentum_cache = -1d0
    local_event_configuration = 0
    local_event_valid = .false.
    parent_born = 0d0

    do node = 1, nlo_decay_node_count()
      node_masses(node) = abs(get_mass_from_id(nlo_decay_node_pdg(node)))
      if (node_masses(node) <= 0d0) then
        call fail_kinematics('a forced decay parent has zero model mass')
      end if
    end do
    do leg = 1, nlo_decay_leaf_count()
      leaf_masses(leg) = abs(get_mass_from_id(nlo_decay_leaf_pdg(leg)))
    end do

    parent_mass = abs(get_mass_from_id(corrected_parent_pdg()))
    if (parent_mass <= 0d0) then
      call fail_kinematics('the corrected parent has zero model mass')
    end if
    do leg = 1, nlo_decay_production_count()
      if (nlo_decay_production_target_kind( &
              nlo_decay_born_context(), leg) == nlo_decay_node_target) then
        target = nlo_decay_production_target_id( &
             nlo_decay_born_context(), leg)
        production_masses(leg) = node_masses(target)
      else
        production_masses(leg) = &
             abs(get_mass_from_id(nlo_decay_production_pdg(leg)))
      end if
    end do
    context = nlo_decay_born_context()
    do leg = 1, nlo_decay_local_count(context)
      if (nlo_decay_local_target_kind(context, leg) == &
          nlo_decay_node_target .and. &
          nlo_decay_local_is_final(context, leg)) then
        target = nlo_decay_local_target_id(context, leg)
        born_local_masses(leg) = node_masses(target)
      else
        born_local_masses(leg) = &
             abs(get_mass_from_id(nlo_decay_local_pdg(context, leg)))
      end if
      if (.not. nlo_decay_local_is_final(context, leg) .and. &
          abs(born_local_masses(leg) - parent_mass) > &
          1d-10*max(1d0, parent_mass)) then
        call fail_kinematics('local incoming mass differs from parent mass')
      end if
    end do
    initialized = .true.
  end subroutine initialize_nlo_decay_kinematics


  subroutine clear_nlo_decay_kinematics()
    if (allocated(production_masses)) deallocate(production_masses)
    if (allocated(born_local_masses)) deallocate(born_local_masses)
    if (allocated(node_masses)) deallocate(node_masses)
    if (allocated(leaf_masses)) deallocate(leaf_masses)
    if (allocated(production_born)) deallocate(production_born)
    if (allocated(born_local)) deallocate(born_local)
    if (allocated(node_random_start)) deallocate(node_random_start)
    if (allocated(local_event_cache)) deallocate(local_event_cache)
    if (allocated(local_fks_momentum_cache)) &
         deallocate(local_fks_momentum_cache)
    if (allocated(local_event_configuration)) &
         deallocate(local_event_configuration)
    if (allocated(local_event_valid)) deallocate(local_event_valid)
    parent_mass = 0d0
    massive_xjac_cache = 1d0
    parent_born = 0d0
  end subroutine clear_nlo_decay_kinematics


  integer function production_random_dimension()
    integer :: final_count
    call require_enabled()
    final_count = nlo_decay_production_count() - nincoming
    if (final_count == 1) then
      production_random_dimension = 0
    else
      production_random_dimension = 3*final_count - 4
    end if
  end function production_random_dimension


  integer function decay_born_random_dimension()
    integer :: node
    call require_enabled()
    decay_born_random_dimension = 0
    do node = 1, nlo_decay_node_count()
      decay_born_random_dimension = decay_born_random_dimension + &
           3*nlo_decay_node_child_count(node) - 4
    end do
  end function decay_born_random_dimension


  double precision function nlo_decay_minimum_production_mass()
    integer :: leg
    call require_enabled()
    nlo_decay_minimum_production_mass = 0d0
    do leg = nincoming + 1, nlo_decay_production_count()
      nlo_decay_minimum_production_mass = &
           nlo_decay_minimum_production_mass + production_masses(leg)
    end do
  end function nlo_decay_minimum_production_mass


  double precision function nlo_decay_production_mass(leg)
    integer, intent(in) :: leg
    call require_enabled()
    if (leg < 1 .or. leg > nlo_decay_production_count()) then
      call fail_kinematics('production mass index is out of range')
    end if
    nlo_decay_production_mass = production_masses(leg)
  end function nlo_decay_production_mass


  subroutine get_nlo_decay_production_momenta(momenta)
    double precision, intent(out) :: momenta(0:3, nexternal)
    call require_enabled()
    momenta = 0d0
    momenta(:, 1:nlo_decay_production_count()) = &
         production_born(:, 1:nlo_decay_production_count())
  end subroutine get_nlo_decay_production_momenta


  double precision function nlo_decay_fks_sister_mass(configuration)
    integer, intent(in) :: configuration
    integer :: context, local_j
    call require_enabled()
    context = nlo_decay_context_for_fks(configuration)
    local_j = nlo_decay_fks_j(configuration)
    nlo_decay_fks_sister_mass = &
         abs(get_mass_from_id(nlo_decay_local_pdg(context, local_j)))
  end function nlo_decay_fks_sister_mass


  double precision function nlo_decay_parent_mass()
    call require_enabled()
    nlo_decay_parent_mass = parent_mass
  end function nlo_decay_parent_mass


  subroutine get_nlo_decay_mass_buffer(configuration, masses)
    integer, intent(in) :: configuration
    double precision, intent(out) :: masses(nexternal)
    integer :: context, leg

    call require_enabled()
    context = nlo_decay_context_for_fks(configuration)
    masses = 0d0
    do leg = 1, nlo_decay_local_count(context)
      masses(leg) = &
           abs(get_mass_from_id(nlo_decay_local_pdg(context, leg)))
    end do
  end subroutine get_nlo_decay_mass_buffer


  subroutine get_nlo_decay_event_momenta(configuration, event_slot, momenta)
    integer, intent(in) :: configuration, event_slot
    double precision, intent(out) :: momenta(0:3, nexternal)

    call require_enabled()
    call check_event_slot(event_slot)
    if (.not. local_event_valid(event_slot) .or. &
        local_event_configuration(event_slot) /= configuration) then
      write (*, *) 'NLO-decay event cache request:', configuration, &
           event_slot, local_event_valid(event_slot), &
           local_event_configuration(event_slot)
      call fail_kinematics('decay-local event momenta are unavailable')
    end if
    momenta = local_event_cache(:, :, event_slot)
  end subroutine get_nlo_decay_event_momenta


  subroutine get_nlo_decay_born_kernel(configuration, momenta, masses, &
                                       particle_count)
    integer, intent(in) :: configuration
    double precision, intent(out) :: momenta(0:3, nexternal)
    double precision, intent(out) :: masses(nexternal)
    integer, intent(out) :: particle_count
    integer :: context, leg, born_leg
    double precision :: inverse_parent(0:3), parent_rest(0:3)

    call require_enabled()
    context = nlo_decay_context_for_fks(configuration)
    particle_count = nlo_decay_local_count(context)
    momenta = 0d0
    masses = 0d0
    parent_rest = 0d0
    parent_rest(0) = parent_mass
    inverse_parent = parent_born
    inverse_parent(1:3) = -inverse_parent(1:3)

    do leg = 1, particle_count
      masses(leg) = &
           abs(get_mass_from_id(nlo_decay_local_pdg(context, leg)))
      if (.not. nlo_decay_local_is_final(context, leg)) then
        momenta(:, leg) = parent_rest
      else if (leg /= nlo_decay_fks_i(configuration)) then
        born_leg = nlo_decay_real_to_born(configuration, leg)
        if (born_leg == 0) then
          call fail_kinematics( &
               'a non-radiated decay leg has no underlying-Born image')
        end if
        call boost_from_rest(born_local(:, born_leg), inverse_parent, &
                             parent_mass, momenta(:, leg))
      end if
    end do
  end subroutine get_nlo_decay_born_kernel


  subroutine get_nlo_decay_counterevent_fks_momenta(configuration, momenta)
    integer, intent(in) :: configuration
    double precision, intent(out) :: momenta(0:3, 0:2)
    integer :: event_slot

    call require_enabled()
    momenta = -1d0
    do event_slot = soft_counterevent, soft_collinear_counterevent
      if (local_event_valid(event_slot) .and. &
          local_event_configuration(event_slot) == configuration) then
        momenta(:, event_slot) = local_fks_momentum_cache(:, event_slot)
      end if
    end do
  end subroutine get_nlo_decay_counterevent_fks_momenta


  subroutine fill_nlo_decay_born_masses(masses)
    double precision, intent(out) :: masses(nexternal - 1)
    integer :: context, leg, target, leaf

    call require_enabled()
    masses = 0d0
    context = nlo_decay_born_context()
    do leg = 1, nlo_decay_production_count()
      if (nlo_decay_production_target_kind(context, leg) == &
          nlo_decay_leg_target) then
        target = nlo_decay_production_target_id(context, leg)
        masses(target) = production_masses(leg)
      end if
    end do
    do leg = 1, nlo_decay_local_count(context)
      if (nlo_decay_local_target_kind(context, leg) == &
          nlo_decay_leg_target) then
        target = nlo_decay_local_target_id(context, leg)
        masses(target) = born_local_masses(leg)
      end if
    end do
    do leaf = 1, nlo_decay_leaf_count()
      target = nlo_decay_leaf_visible(context, leaf)
      if (target /= 0) masses(target) = leaf_masses(leaf)
    end do
  end subroutine fill_nlo_decay_born_masses


  subroutine fill_nlo_decay_event_masses(configuration, masses)
    integer, intent(in) :: configuration
    double precision, intent(out) :: masses(nexternal)
    integer :: context, leg, target, leaf
    logical :: covered(nexternal)

    call require_enabled()
    masses = 0d0
    covered = .false.
    context = nlo_decay_context_for_fks(configuration)
    do leg = 1, nlo_decay_production_count()
      if (nlo_decay_production_target_kind(context, leg) == &
          nlo_decay_leg_target) then
        target = nlo_decay_production_target_id(context, leg)
        masses(target) = production_masses(leg)
        covered(target) = .true.
      end if
    end do
    do leg = 1, nlo_decay_local_count(context)
      if (nlo_decay_local_target_kind(context, leg) == &
          nlo_decay_leg_target) then
        target = nlo_decay_local_target_id(context, leg)
        masses(target) = &
             abs(get_mass_from_id(nlo_decay_local_pdg(context, leg)))
        covered(target) = .true.
      end if
    end do
    do leaf = 1, nlo_decay_leaf_count()
      target = nlo_decay_leaf_visible(context, leaf)
      if (target /= 0) then
        masses(target) = leaf_masses(leaf)
        covered(target) = .true.
      end if
    end do
    if (.not. all(covered)) then
      call fail_kinematics('NLO-decay event masses do not cover all legs')
    end if
  end subroutine fill_nlo_decay_event_masses


  subroutine generate_nlo_decay_born_momenta(x, shat, sqrtshat, xjac, &
                                               xpswgt, visible, pass)
    double precision, intent(in) :: x(99), shat, sqrtshat
    double precision, intent(inout) :: xjac, xpswgt
    double precision, intent(out) :: visible(0:3, nexternal - 1)
    logical, intent(out) :: pass

    integer :: context, final_count, leg, decay_index
    double precision :: system(0:3), final_masses(nexternal)
    double precision :: final_momenta(0:3, nexternal)

    call require_enabled()
    pass = .false.
    visible = 0d0
    production_born = 0d0
    born_local = 0d0
    node_random_start = 0
    parent_born = 0d0
    local_event_cache = 0d0
    local_fks_momentum_cache = -1d0
    local_event_configuration = 0
    local_event_valid = .false.
    context = nlo_decay_born_context()

    system = 0d0
    system(0) = sqrtshat
    final_count = nlo_decay_production_count() - nincoming
    do leg = 1, final_count
      final_masses(leg) = production_masses(nincoming + leg)
    end do
    call generate_nbody(system, final_count, final_masses, x, 1, &
                        final_momenta, xjac, xpswgt, pass)
    if (.not. pass) return
    call fill_production_incoming(shat, sqrtshat, production_born, pass)
    if (.not. pass) return
    do leg = 1, final_count
      production_born(:, nincoming + leg) = final_momenta(:, leg)
    end do

    decay_index = production_random_dimension() + 1
    call expand_born_context(context, x, decay_index, visible, xjac, &
                             xpswgt, pass)
    if (.not. pass) return
    if (abs(minkowski_square(parent_born) - parent_mass**2) > &
        1d-8*max(1d0, parent_mass**2)) then
      call fail_kinematics('generated corrected parent is off shell')
    end if
    if (decay_index /= production_random_dimension() + &
        decay_born_random_dimension() + 1) then
      call fail_kinematics('Born random-variable accounting is inconsistent')
    end if

    xpswgt = xpswgt*nlo_decay_nwa_weight()
    pass = .true.
  end subroutine generate_nlo_decay_born_momenta


  subroutine expand_born_context(context, x, index, visible, xjac, &
                                 xpswgt, pass)
    integer, intent(in) :: context
    double precision, intent(in) :: x(99)
    integer, intent(inout) :: index
    double precision, intent(out) :: visible(0:3, nexternal - 1)
    double precision, intent(inout) :: xjac, xpswgt
    logical, intent(out) :: pass
    integer :: leg, target

    call set_nlo_decay_cut_mask(context)
    visible = 0d0
    pass = .true.
    do leg = 1, nlo_decay_production_count()
      target = nlo_decay_production_target_id(context, leg)
      if (nlo_decay_production_target_kind(context, leg) == &
          nlo_decay_leg_target) then
        visible(:, target) = production_born(:, leg)
      else
        call expand_born_node(context, target, production_born(:, leg), x, &
             index, visible, xjac, xpswgt, pass)
        if (.not. pass) return
      end if
    end do
  end subroutine expand_born_context


  subroutine set_nlo_decay_cut_mask(context)
    integer, intent(in) :: context
    integer :: leaf, leg, visible_leg

    call require_enabled()
    event_from_decay = .false.

    ! The topology leaves cover every ordinary LO decay in the complete
    ! forest, including ancestors, siblings and descendants of the corrected
    ! node.
    do leaf = 1, nlo_decay_leaf_count()
      visible_leg = nlo_decay_leaf_visible(context, leaf)
      if (visible_leg > nexternal) then
        call fail_kinematics('a decay leaf has an invalid visible target')
      end if
      if (visible_leg > nincoming) then
        event_from_decay(visible_leg) = .true.
      end if
    end do

    ! Direct final legs of the corrected local process are not necessarily
    ! topology leaves.  In particular, the real-emission parton exists only
    ! in a real context and must also disappear from generation cuts when
    ! cut_decays is disabled.
    do leg = 1, nlo_decay_local_count(context)
      if (.not. nlo_decay_local_is_final(context, leg)) cycle
      if (nlo_decay_local_target_kind(context, leg) /= &
          nlo_decay_leg_target) cycle
      visible_leg = nlo_decay_local_target_id(context, leg)
      if (visible_leg <= nincoming .or. visible_leg > nexternal) then
        call fail_kinematics( &
             'a corrected-decay leg has an invalid visible target')
      end if
      event_from_decay(visible_leg) = .true.
    end do
  end subroutine set_nlo_decay_cut_mask


  recursive subroutine expand_born_node(context, node, parent, x, index, &
                                        visible, xjac, xpswgt, pass)
    integer, intent(in) :: context, node
    double precision, intent(in) :: parent(0:3), x(99)
    integer, intent(inout) :: index
    double precision, intent(inout) :: visible(0:3, nexternal - 1)
    double precision, intent(inout) :: xjac, xpswgt
    logical, intent(out) :: pass
    integer :: child_count, child, child_kind, identifier
    integer :: leg, final_count, final_index, target
    double precision :: child_masses(nexternal)
    double precision :: child_momenta(0:3, nexternal)

    node_random_start(node) = index
    if (node == nlo_decay_corrected_node()) then
      final_count = 0
      do leg = 1, nlo_decay_local_count(context)
        if (nlo_decay_local_is_final(context, leg)) then
          final_count = final_count + 1
          child_masses(final_count) = born_local_masses(leg)
        else
          born_local(:, leg) = parent
        end if
      end do
      if (final_count /= nlo_decay_node_child_count(node)) then
        call fail_kinematics('corrected-node child count is inconsistent')
      end if
      call generate_nbody(parent, final_count, child_masses, x, index, &
           child_momenta, xjac, xpswgt, pass)
      if (.not. pass) return
      index = index + 3*final_count - 4
      final_index = 0
      do leg = 1, nlo_decay_local_count(context)
        if (.not. nlo_decay_local_is_final(context, leg)) cycle
        final_index = final_index + 1
        born_local(:, leg) = child_momenta(:, final_index)
        target = nlo_decay_local_target_id(context, leg)
        if (nlo_decay_local_target_kind(context, leg) == &
            nlo_decay_node_target) then
          call expand_born_node(context, target, born_local(:, leg), x, &
               index, visible, xjac, xpswgt, pass)
          if (.not. pass) return
        else
          visible(:, target) = born_local(:, leg)
        end if
      end do
      parent_born = parent
      return
    end if

    child_count = nlo_decay_node_child_count(node)
    do child = 1, child_count
      identifier = nlo_decay_node_child_id(node, child)
      child_kind = nlo_decay_node_child_kind(node, child)
      if (child_kind == nlo_decay_node_child) then
        child_masses(child) = node_masses(identifier)
      else
        child_masses(child) = leaf_masses(identifier)
      end if
    end do
    call generate_nbody(parent, child_count, child_masses, x, index, &
         child_momenta, xjac, xpswgt, pass)
    if (.not. pass) return
    index = index + 3*child_count - 4
    do child = 1, child_count
      identifier = nlo_decay_node_child_id(node, child)
      child_kind = nlo_decay_node_child_kind(node, child)
      if (child_kind == nlo_decay_node_child) then
        call expand_born_node(context, identifier, child_momenta(:, child), &
             x, index, visible, xjac, xpswgt, pass)
        if (.not. pass) return
      else
        target = nlo_decay_leaf_visible(context, identifier)
        visible(:, target) = child_momenta(:, child)
      end if
    end do
  end subroutine expand_born_node


  subroutine generate_nlo_decay_fks_event(configuration, event_slot, x, ndim, &
       xjac, xpswgt, visible, xiimax, xinorm, xi_i, xi_hat, y_ij, &
       p_i_hat, solution_sign, pass)
    integer, intent(in) :: configuration, event_slot, ndim
    double precision, intent(in) :: x(99)
    double precision, intent(inout) :: xjac, xpswgt, xi_hat
    double precision, intent(out) :: visible(0:3, nexternal)
    double precision, intent(out) :: xiimax, xinorm, xi_i, y_ij
    double precision, intent(out) :: p_i_hat(0:3)
    integer, intent(out) :: solution_sign
    logical, intent(out) :: pass

    integer :: context, local_i, local_j, local_ij, leg, born_leg
    double precision :: inverse_parent(0:3), parent_rest(0:3)
    double precision :: local_rest(0:3, nexternal)
    double precision :: local_event(0:3, nexternal)
    double precision :: born_emitter(0:3), recoil(0:3), recoil_mass2
    double precision :: sister_mass, phi, reference_rest(0:3)

    call require_enabled()
    call check_event_slot(event_slot)
    local_event_valid(event_slot) = .false.
    local_event_configuration(event_slot) = 0
    local_event_cache(:, :, event_slot) = 0d0
    local_fks_momentum_cache(:, event_slot) = -1d0
    pass = .false.
    visible = 0d0
    local_rest = 0d0
    local_event = 0d0
    context = nlo_decay_context_for_fks(configuration)
    local_i = nlo_decay_fks_i(configuration)
    local_j = nlo_decay_fks_j(configuration)
    local_ij = nlo_decay_fks_ij(configuration)
    sister_mass = nlo_decay_fks_sister_mass(configuration)
    if (abs(get_mass_from_id(nlo_decay_local_pdg(context, local_i))) > 0d0) then
      call fail_kinematics('the emitted NLO-decay FKS leg must be massless')
    end if

    ! The three FKS variables live in the corrected-parent rest frame.
    ! xi=2 E_i*/M, y=cos(theta_ij*), and phi
    ! is the azimuth around the underlying-Born ij direction in that frame.
    ! No production momentum participates in this construction.
    phi = 2d0*pi_value()*x(ndim)
    xjac = xjac*2d0*pi_value()
    parent_rest = 0d0
    parent_rest(0) = parent_mass
    inverse_parent = parent_born
    inverse_parent(1:3) = -inverse_parent(1:3)

    do leg = 1, nlo_decay_local_count(context)
      if (.not. nlo_decay_local_is_final(context, leg)) then
        local_rest(:, leg) = parent_rest
      else if (leg /= local_i) then
        born_leg = nlo_decay_real_to_born(configuration, leg)
        call boost_from_rest(born_local(:, born_leg), inverse_parent, &
                             parent_mass, local_rest(:, leg))
      end if
    end do
    call boost_from_rest(born_local(:, local_ij), inverse_parent, &
                         parent_mass, born_emitter)
    recoil = parent_rest - born_emitter
    recoil_mass2 = minkowski_square(recoil)
    if (recoil_mass2 < -1d-10*parent_mass**2) then
      call fail_kinematics('decay-local recoil has negative invariant mass')
    end if
    recoil_mass2 = max(0d0, recoil_mass2)

    if (sister_mass == 0d0) then
      call map_massless_decay_fks(context, event_slot, local_i, local_j, &
           born_emitter, recoil_mass2, x(ndim - 2), x(ndim - 1), phi, &
           local_rest, xiimax, xinorm, xi_i, xi_hat, y_ij, p_i_hat, &
           xjac, xpswgt, solution_sign, pass)
    else
      call map_massive_decay_fks(context, event_slot, local_i, local_j, &
           born_emitter, recoil_mass2, sister_mass, x(ndim - 2), &
           x(ndim - 1), phi, local_rest, xiimax, xinorm, xi_i, xi_hat, &
           y_ij, p_i_hat, xjac, xpswgt, solution_sign, pass)
    end if
    if (.not. pass) return

    do leg = 1, nlo_decay_local_count(context)
      if (nlo_decay_local_is_final(context, leg)) then
        call boost_from_rest(local_rest(:, leg), parent_born, parent_mass, &
                             local_event(:, leg))
      else
        local_event(:, leg) = parent_born
      end if
    end do
    reference_rest = p_i_hat
    call boost_from_rest(reference_rest, parent_born, parent_mass, p_i_hat)
    call validate_local_recoil(context, local_event, pass)
    if (.not. pass) return
    local_event_cache(:, :, event_slot) = local_rest
    local_fks_momentum_cache(:, event_slot) = reference_rest
    local_event_configuration(event_slot) = configuration
    local_event_valid(event_slot) = .true.
    call expand_event_context(context, local_event, x, visible, pass)
    if (.not. pass) return
    pass = .true.
  end subroutine generate_nlo_decay_fks_event


  subroutine map_massless_decay_fks(context, event_slot, local_i, local_j, &
       born_emitter, recoil_mass2, random_xi, random_y, phi, momenta, &
       xiimax, xinorm, xi_i, xi_hat, y_ij, p_i_hat, xjac, xpswgt, &
       solution_sign, pass)
    integer, intent(in) :: context, event_slot, local_i, local_j
    double precision, intent(in) :: born_emitter(0:3), recoil_mass2
    double precision, intent(in) :: random_xi, random_y, phi
    double precision, intent(inout) :: momenta(0:3, nexternal)
    double precision, intent(out) :: xiimax, xinorm, xi_i, y_ij
    double precision, intent(inout) :: xi_hat, xjac, xpswgt
    double precision, intent(out) :: p_i_hat(0:3)
    integer, intent(out) :: solution_sign
    logical, intent(out) :: pass

    double precision :: sstiny, cctiny, emitted_energy, sister_length
    double precision :: mother_length, cos_i, sin_i, direction(0:3)
    double precision :: mother_theta, mother_cos, mother_sin, mother_phi
    double precision :: mother_cosphi, mother_sinphi, mother(0:3)
    double precision :: recoil(0:3), sumrec, sumrec2, shy, chy, chymo
    double precision :: boost_direction(3), sister_norm, born_norm
    integer :: leg
    double precision, parameter :: stiny = 1d-6, ctiny = 5d-7
    complex(kind=kind(0d0)), parameter :: imaginary = (0d0, 1d0)

    pass = .false.
    solution_sign = 1
    sstiny = merge(0d0, stiny, softtest)
    cctiny = merge(0d0, ctiny, colltest)
    call choose_y(event_slot, random_y, cctiny, y_ij, xjac, pass)
    if (.not. pass) return
    call get_angles(born_emitter, mother_theta, mother_cos, mother_sin, &
                    mother_phi, mother_cosphi, mother_sinphi)

    xiimax = 1d0 - recoil_mass2/parent_mass**2
    xinorm = xiimax
    call choose_massless_xi(event_slot, random_xi, sstiny, xiimax, &
                            xi_hat, xi_i, xjac, pass)
    if (.not. pass) return

    emitted_energy = xi_i*parent_mass/2d0
    sister_length = (parent_mass**2 - recoil_mass2 - &
         2d0*parent_mass*emitted_energy)/(2d0*(parent_mass - &
         emitted_energy*(1d0 - y_ij)))
    mother_length = sqrt(max(0d0, emitted_energy**2 + sister_length**2 + &
                         2d0*emitted_energy*sister_length*y_ij))
    call emitted_angle(xi_i, y_ij, recoil_mass2, emitted_energy, &
                       sister_length, mother_length, cos_i)
    sin_i = sqrt(max(0d0, 1d0 - cos_i**2))
    direction(0) = 1d0
    direction(1) = sin_i*cos(phi)
    direction(2) = sin_i*sin(phi)
    direction(3) = cos_i

    momenta(0, local_i) = emitted_energy
    momenta(1:3, local_i) = emitted_energy*direction(1:3)
    momenta(0, local_j) = sister_length
    momenta(1, local_j) = -momenta(1, local_i)
    momenta(2, local_j) = -momenta(2, local_i)
    momenta(3, local_j) = mother_length - momenta(3, local_i)
    p_i_hat(0) = parent_mass/2d0
    p_i_hat(1:3) = parent_mass/2d0*direction(1:3)
    call rotate_invar(momenta(:, local_i), momenta(:, local_i), &
                      mother_cos, mother_sin, mother_cosphi, mother_sinphi)
    call rotate_invar(momenta(:, local_j), momenta(:, local_j), &
                      mother_cos, mother_sin, mother_cosphi, mother_sinphi)
    call rotate_invar(p_i_hat, p_i_hat, mother_cos, mother_sin, &
                      mother_cosphi, mother_sinphi)

    mother = momenta(:, local_i) + momenta(:, local_j)
    recoil = 0d0
    recoil(0) = parent_mass
    recoil = recoil - mother
    sumrec = recoil(0) + norm3(recoil)
    sumrec2 = sumrec**2
    if (sumrec <= 0d0) return
    shy = -(parent_mass**2 - sumrec2)/(2d0*sumrec*parent_mass)
    chy = (parent_mass**2 + sumrec2)/(2d0*sumrec*parent_mass)
    chymo = (parent_mass - sumrec)**2/(2d0*sumrec*parent_mass)
    if (mother_length > 0d0) then
      boost_direction = mother(1:3)/mother_length
    else
      boost_direction = (/0d0, 0d0, 1d0/)
    end if
    do leg = 1, nlo_decay_local_count(context)
      if (leg /= local_i .and. leg /= local_j .and. &
          nlo_decay_local_is_final(context, leg) .and. &
          shy /= 0d0) then
        call boostwdir2_in_place(chy, shy, chymo, boost_direction, &
                                 momenta(:, leg))
      end if
    end do

    if (event_slot == real_event .or. &
        (event_slot == collinear_counterevent .and. xij_aor == 0)) then
      xij_aor = -exp(2d0*imaginary*(mother_phi + phi))
    end if
    sister_norm = norm3(momenta(:, local_j))
    born_norm = norm3(born_emitter)
    if (sister_norm <= 0d0 .or. born_norm <= 0d0) return
    xpswgt = xpswgt*2d0*parent_mass**2/(4d0*pi_value())**3* &
         sister_norm/born_norm/(2d0 - xi_i*(1d0 - y_ij))
    xpswgt = abs(xpswgt)
    pass = .true.
  end subroutine map_massless_decay_fks


  subroutine map_massive_decay_fks(context, event_slot, local_i, local_j, &
       born_emitter, recoil_mass2, sister_mass, random_xi, random_y, phi, &
       momenta, xiimax, xinorm, xi_i, xi_hat, y_ij, p_i_hat, xjac, &
       xpswgt, solution_sign, pass)
    integer, intent(in) :: context, event_slot, local_i, local_j
    double precision, intent(in) :: born_emitter(0:3), recoil_mass2
    double precision, intent(in) :: sister_mass, random_xi, random_y, phi
    double precision, intent(inout) :: momenta(0:3, nexternal)
    double precision, intent(out) :: xiimax, xinorm, xi_i, y_ij
    double precision, intent(inout) :: xi_hat, xjac, xpswgt
    double precision, intent(out) :: p_i_hat(0:3)
    integer, intent(out) :: solution_sign
    logical, intent(out) :: pass

    double precision :: sstiny, cctiny, mjhat, rechat, xim, xi_boundary
    double precision :: xi_absolute_max, xi_minus, ratio, temp
    double precision :: cffa, cffb, cffc, discriminant, root_expression
    double precision :: emitted_energy, sister_num, sister_den, sister_length
    double precision :: mother_length, cos_i, sin_i, direction(0:3)
    double precision :: mother_theta, mother_cos, mother_sin, mother_phi
    double precision :: mother_cosphi, mother_sinphi, mother(0:3)
    double precision :: recoil(0:3), sumrec, expy, shy, chy, chymo
    double precision :: boost_direction(3), sister_norm, born_norm
    integer :: leg
    double precision, parameter :: stiny = 1d-6, ctiny = 5d-7
    double precision, parameter :: qtiny = 1d-7

    pass = .false.
    if (colltest .or. event_slot == collinear_counterevent .or. &
        event_slot == soft_collinear_counterevent) return
    sstiny = merge(0d0, stiny, softtest)
    cctiny = merge(0d0, ctiny, colltest)
    call choose_y(event_slot, random_y, cctiny, y_ij, xjac, pass)
    if (.not. pass) return
    call get_angles(born_emitter, mother_theta, mother_cos, mother_sin, &
                    mother_phi, mother_cosphi, mother_sinphi)

    mjhat = sister_mass/parent_mass
    rechat = sqrt(recoil_mass2)/parent_mass
    xim = (1d0 - rechat**2 - 2d0*mjhat + mjhat**2)/(1d0 - mjhat)
    cffa = 1d0 - mjhat**2*(1d0 - y_ij**2)
    cffb = -2d0*(1d0 - rechat**2 - mjhat**2)
    cffc = (1d0 - (rechat - mjhat)**2)*(1d0 - (rechat + mjhat)**2)
    discriminant = max(0d0, cffb**2 - 4d0*cffa*cffc)
    xi_boundary = (-cffb - sqrt(discriminant))/(2d0*cffa)
    xi_absolute_max = 1d0 - (rechat + mjhat)**2
    if (xim < 0d0 .or. xi_boundary < -1d-10 .or. &
        xi_boundary > xi_absolute_max + 1d-8) return
    if (y_ij >= 0d0) then
      xiimax = xim
      xi_minus = 0d0
    else
      xiimax = xi_boundary
      xi_minus = xi_boundary - xim
    end if
    xinorm = xiimax + xi_minus
    if (xinorm <= 0d0) return
    ratio = xiimax/xinorm
    call choose_massive_xi(event_slot, random_xi, sstiny, ratio, xiimax, &
                           xinorm, xi_hat, xi_i, solution_sign, xjac, pass)
    if (.not. pass) return

    emitted_energy = xi_i*parent_mass/2d0
    root_expression = xi_i**2*cffa + xi_i*cffb + cffc
    if (root_expression < -1d-8) return
    root_expression = max(0d0, root_expression)
    sister_num = -xi_i*y_ij*(1d0 - rechat**2 + mjhat**2 - xi_i) + &
         (2d0 - xi_i)*sqrt(root_expression)*solution_sign
    sister_den = (2d0 - xi_i*(1d0 - y_ij))* &
                 (2d0 - xi_i*(1d0 + y_ij))
    sister_length = parent_mass*sister_num/sister_den
    if (sister_length < 0d0) return
    mother_length = sqrt(max(0d0, emitted_energy**2 + sister_length**2 + &
                         2d0*emitted_energy*sister_length*y_ij))
    if (xi_i < qtiny) then
      cos_i = y_ij + (1d0 - y_ij**2)*xi_i/sqrt(max(qtiny, cffc))
    else
      cos_i = (mother_length**2 - sister_length**2 + emitted_energy**2)/ &
              (2d0*mother_length*emitted_energy)
    end if
    if (abs(cos_i) > 1d0 + 1d-7) return
    cos_i = max(-1d0, min(1d0, cos_i))
    sin_i = sqrt(max(0d0, 1d0 - cos_i**2))
    direction(0) = 1d0
    direction(1) = sin_i*cos(phi)
    direction(2) = sin_i*sin(phi)
    direction(3) = cos_i

    momenta(0, local_i) = emitted_energy
    momenta(1:3, local_i) = emitted_energy*direction(1:3)
    momenta(0, local_j) = sqrt(sister_length**2 + sister_mass**2)
    momenta(1, local_j) = -momenta(1, local_i)
    momenta(2, local_j) = -momenta(2, local_i)
    momenta(3, local_j) = mother_length - momenta(3, local_i)
    p_i_hat(0) = parent_mass/2d0
    p_i_hat(1:3) = parent_mass/2d0*direction(1:3)
    call rotate_invar(momenta(:, local_i), momenta(:, local_i), &
                      mother_cos, mother_sin, mother_cosphi, mother_sinphi)
    call rotate_invar(momenta(:, local_j), momenta(:, local_j), &
                      mother_cos, mother_sin, mother_cosphi, mother_sinphi)
    call rotate_invar(p_i_hat, p_i_hat, mother_cos, mother_sin, &
                      mother_cosphi, mother_sinphi)

    mother = momenta(:, local_i) + momenta(:, local_j)
    recoil = 0d0
    recoil(0) = parent_mass
    recoil = recoil - mother
    sumrec = recoil(0) + norm3(recoil)
    if (sumrec <= 0d0) return
    if (recoil_mass2 < 1d-16*parent_mass**2) then
      expy = parent_mass*sumrec/(parent_mass**2 - sister_mass**2)* &
           (1d0 + sister_mass**2*recoil_mass2/ &
            (parent_mass**2 - sister_mass**2)**2)
    else
      expy = sumrec/(2d0*parent_mass*recoil_mass2)*(parent_mass**2 + &
           recoil_mass2 - sister_mass**2 - parent_mass**2*sqrt(cffc))
    end if
    if (expy <= 0d0) return
    shy = (expy - 1d0/expy)/2d0
    chy = (expy + 1d0/expy)/2d0
    chymo = chy - 1d0
    if (mother_length > 0d0) then
      boost_direction = mother(1:3)/mother_length
    else
      boost_direction = (/0d0, 0d0, 1d0/)
    end if
    do leg = 1, nlo_decay_local_count(context)
      if (leg /= local_i .and. leg /= local_j .and. &
          nlo_decay_local_is_final(context, leg) .and. shy /= 0d0) then
        call boostwdir2_in_place(chy, shy, chymo, boost_direction, &
                                 momenta(:, leg))
      end if
    end do

    sister_norm = norm3(momenta(:, local_j))
    born_norm = norm3(born_emitter)
    if (sister_norm <= 0d0 .or. born_norm <= 0d0) return
    temp = momenta(0, local_j)/sister_norm
    xpswgt = xpswgt*2d0*parent_mass**2/(4d0*pi_value())**3* &
         sister_norm/born_norm/(2d0 - xi_i*(1d0 - temp*y_ij))
    xpswgt = abs(xpswgt)
    pass = .true.
  end subroutine map_massive_decay_fks


  subroutine choose_y(event_slot, random_y, cutoff, y_ij, xjac, pass)
    integer, intent(in) :: event_slot
    double precision, intent(in) :: random_y, cutoff
    double precision, intent(out) :: y_ij
    double precision, intent(inout) :: xjac
    logical, intent(out) :: pass

    pass = .true.
    if ((event_slot == real_event .or. event_slot == soft_counterevent) .and. &
        ((.not. softtest) .or. y_ij_fks_fix == -2d0) .and. &
        (.not. colltest)) then
      y_ij = 1d0 - 2d0*(cutoff + (1d0 - cutoff)*random_y**2)
    else if ((event_slot == real_event .or. &
              event_slot == soft_counterevent) .and. &
             ((softtest .and. y_ij_fks_fix /= -2d0) .or. colltest)) then
      y_ij = y_ij_fks_fix
    else if (event_slot == collinear_counterevent .or. &
             event_slot == soft_collinear_counterevent) then
      y_ij = 1d0
    else
      pass = .false.
      return
    end if
    xjac = xjac*4d0*random_y
  end subroutine choose_y


  subroutine choose_massless_xi(event_slot, random_xi, cutoff, xiimax, &
                                xi_hat, xi_i, xjac, pass)
    integer, intent(in) :: event_slot
    double precision, intent(in) :: random_xi, cutoff, xiimax
    double precision, intent(inout) :: xi_hat, xjac
    double precision, intent(out) :: xi_i
    logical, intent(out) :: pass

    pass = .true.
    if ((event_slot == real_event .or. &
         event_slot == collinear_counterevent) .and. &
        ((.not. colltest) .or. xi_i_fks_fix == -2d0) .and. &
        (.not. softtest)) then
      if (event_slot == real_event) then
        xi_hat = cutoff + (1d0 - cutoff)*random_xi**2
      end if
      xi_i = xi_hat*xiimax
    else if ((event_slot == real_event .or. &
              event_slot == collinear_counterevent) .and. &
             colltest .and. xi_i_fks_fix /= -2d0 .and. &
             (.not. softtest)) then
      xi_i = xi_i_fks_fix*xiimax
    else if ((event_slot == real_event .or. &
              event_slot == collinear_counterevent) .and. softtest) then
      if (xi_i_fks_fix >= 1d0) then
        pass = .false.
        return
      end if
      xi_i = xi_i_fks_fix*xiimax
    else if (event_slot == soft_counterevent .or. &
             event_slot == soft_collinear_counterevent) then
      xi_i = 0d0
    else
      pass = .false.
      return
    end if
    if (xi_i < 0d0 .or. xi_i > xiimax + 1d-12) then
      pass = .false.
      return
    end if
    xjac = xjac*2d0*random_xi
  end subroutine choose_massless_xi


  subroutine choose_massive_xi(event_slot, random_xi, cutoff, ratio, &
       xiimax, xinorm, xi_hat, xi_i, solution_sign, xjac, pass)
    integer, intent(in) :: event_slot
    double precision, intent(in) :: random_xi, cutoff, ratio, xiimax, xinorm
    double precision, intent(inout) :: xi_hat, xjac
    double precision, intent(out) :: xi_i
    integer, intent(out) :: solution_sign
    logical, intent(out) :: pass
    double precision :: mapped

    pass = .true.
    solution_sign = 1
    if (event_slot == real_event .and. (.not. softtest)) then
      massive_xjac_cache = 1d0
      if (random_xi <= ratio) then
        if (ratio <= 0d0) then
          pass = .false.
          return
        end if
        mapped = random_xi/ratio
        massive_xjac_cache = massive_xjac_cache/ratio
        xi_hat = (cutoff + (1d0 - cutoff)*mapped**2)*ratio
        massive_xjac_cache = massive_xjac_cache*2d0*mapped*ratio
        xi_i = xinorm*xi_hat
        solution_sign = 1
      else
        xi_hat = random_xi
        xi_i = -xinorm*xi_hat + 2d0*xiimax
        solution_sign = -1
      end if
    else if (event_slot == real_event .and. softtest) then
      massive_xjac_cache = 1d0
      xi_i = xi_i_fks_fix
      if (xi_i >= xiimax) then
        pass = .false.
        return
      end if
    else if (event_slot == soft_counterevent) then
      xi_i = 0d0
      solution_sign = 1
    else
      pass = .false.
      return
    end if
    xjac = xjac*massive_xjac_cache
  end subroutine choose_massive_xi


  subroutine emitted_angle(xi_i, y_ij, recoil_mass2, emitted_energy, &
                           sister_length, mother_length, cos_i)
    double precision, intent(in) :: xi_i, y_ij, recoil_mass2
    double precision, intent(in) :: emitted_energy, sister_length
    double precision, intent(in) :: mother_length
    double precision, intent(out) :: cos_i
    double precision, parameter :: qtiny = 1d-7

    if (xi_i < qtiny) then
      cos_i = y_ij + parent_mass**2*(1d0 - y_ij**2)*xi_i/ &
              (parent_mass**2 - recoil_mass2)
    else if (1d0 - y_ij < qtiny) then
      cos_i = 1d0 - (parent_mass**2*(1d0 - xi_i) - recoil_mass2)**2* &
              (1d0 - y_ij)/(parent_mass**2 - recoil_mass2)**2
    else
      cos_i = (mother_length**2 - sister_length**2 + emitted_energy**2)/ &
              (2d0*mother_length*emitted_energy)
    end if
    if (abs(cos_i) <= 1d0 + 1d-6) then
      cos_i = max(-1d0, min(1d0, cos_i))
    else
      call fail_kinematics('invalid emitted-parton polar angle')
    end if
  end subroutine emitted_angle


  subroutine expand_event_context(context, local_momenta, x, visible, pass)
    integer, intent(in) :: context
    double precision, intent(in) :: local_momenta(0:3, nexternal)
    double precision, intent(in) :: x(99)
    double precision, intent(out) :: visible(0:, :)
    integer :: leg, target
    logical, intent(out) :: pass

    call set_nlo_decay_cut_mask(context)
    visible = 0d0
    pass = .true.
    do leg = 1, nlo_decay_production_count()
      if (nlo_decay_production_target_kind(context, leg) == &
          nlo_decay_leg_target) then
        target = nlo_decay_production_target_id(context, leg)
        visible(:, target) = production_born(:, leg)
      else
        target = nlo_decay_production_target_id(context, leg)
        call expand_event_node(context, target, production_born(:, leg), &
             local_momenta, x, visible, pass)
        if (.not. pass) return
      end if
    end do
  end subroutine expand_event_context


  recursive subroutine expand_event_node(context, node, parent, &
                                         local_momenta, x, visible, pass)
    integer, intent(in) :: context, node
    double precision, intent(in) :: parent(0:3)
    double precision, intent(in) :: local_momenta(0:3, nexternal), x(99)
    double precision, intent(inout) :: visible(0:, :)
    logical, intent(out) :: pass
    integer :: leg, target, child, child_count, child_kind, identifier
    double precision :: child_masses(nexternal)
    double precision :: child_momenta(0:3, nexternal)
    double precision :: unused_jacobian, unused_weight

    pass = .true.
    if (node == nlo_decay_corrected_node()) then
      do leg = 1, nlo_decay_local_count(context)
        if (.not. nlo_decay_local_is_final(context, leg)) cycle
        target = nlo_decay_local_target_id(context, leg)
        if (nlo_decay_local_target_kind(context, leg) == &
            nlo_decay_node_target) then
          call expand_event_node(context, target, local_momenta(:, leg), &
               local_momenta, x, visible, pass)
          if (.not. pass) return
        else
          visible(:, target) = local_momenta(:, leg)
        end if
      end do
      return
    end if

    child_count = nlo_decay_node_child_count(node)
    do child = 1, child_count
      identifier = nlo_decay_node_child_id(node, child)
      child_kind = nlo_decay_node_child_kind(node, child)
      if (child_kind == nlo_decay_node_child) then
        child_masses(child) = node_masses(identifier)
      else
        child_masses(child) = leaf_masses(identifier)
      end if
    end do
    unused_jacobian = 1d0
    unused_weight = 1d0
    call generate_nbody(parent, child_count, child_masses, x, &
         node_random_start(node), child_momenta, unused_jacobian, &
         unused_weight, pass)
    if (.not. pass) return
    do child = 1, child_count
      identifier = nlo_decay_node_child_id(node, child)
      child_kind = nlo_decay_node_child_kind(node, child)
      if (child_kind == nlo_decay_node_child) then
        call expand_event_node(context, identifier, child_momenta(:, child), &
             local_momenta, x, visible, pass)
        if (.not. pass) return
      else
        target = nlo_decay_leaf_visible(context, identifier)
        visible(:, target) = child_momenta(:, child)
      end if
    end do
  end subroutine expand_event_node


  subroutine fill_production_incoming(shat, sqrtshat, momenta, pass)
    double precision, intent(in) :: shat, sqrtshat
    double precision, intent(inout) :: momenta(0:3, nexternal)
    logical, intent(out) :: pass
    double precision :: lambda_value, momentum_length, mass1, mass2

    pass = .false.
    if (nincoming == 1) then
      momenta(:, 1) = 0d0
      momenta(0, 1) = sqrtshat
      pass = .true.
      return
    end if
    mass1 = production_masses(1)
    mass2 = production_masses(2)
    lambda_value = phase_space_lambda(shat, mass1**2, mass2**2)
    if (lambda_value < 0d0) return
    momentum_length = sqrt(lambda_value)/(2d0*sqrtshat)
    momenta(:, 1) = 0d0
    momenta(:, 2) = 0d0
    momenta(0, 1) = (shat + mass1**2 - mass2**2)/(2d0*sqrtshat)
    momenta(0, 2) = (shat + mass2**2 - mass1**2)/(2d0*sqrtshat)
    momenta(3, 1) = momentum_length
    momenta(3, 2) = -momentum_length
    pass = .true.
  end subroutine fill_production_incoming


  subroutine validate_local_recoil(context, local_momenta, pass)
    integer, intent(in) :: context
    double precision, intent(in) :: local_momenta(0:3, nexternal)
    logical, intent(out) :: pass
    double precision :: total(0:3), scale
    integer :: leg

    total = 0d0
    do leg = 1, nlo_decay_local_count(context)
      if (nlo_decay_local_is_final(context, leg)) then
        total = total + local_momenta(:, leg)
      end if
    end do
    scale = max(1d0, abs(parent_born(0)))
    pass = all(abs(total - parent_born) <= 2d-9*scale)
  end subroutine validate_local_recoil


  subroutine get_angles(momentum, theta, costheta, sintheta, phi, &
                        cosphi, sinphi)
    double precision, intent(in) :: momentum(0:3)
    double precision, intent(out) :: theta, costheta, sintheta, phi
    double precision, intent(out) :: cosphi, sinphi
    double precision :: length

    length = norm3(momentum)
    if (length == 0d0) then
      theta = 0d0
      costheta = 1d0
      sintheta = 0d0
      phi = 0d0
      cosphi = 1d0
      sinphi = 0d0
      return
    end if
    costheta = max(-1d0, min(1d0, momentum(3)/length))
    theta = acos(costheta)
    sintheta = sqrt(max(0d0, 1d0 - costheta**2))
    if (sintheta == 0d0) then
      phi = 0d0
      cosphi = 1d0
      sinphi = 0d0
    else
      phi = atan2(momentum(2), momentum(1))
      cosphi = cos(phi)
      sinphi = sin(phi)
    end if
  end subroutine get_angles


  double precision function nlo_decay_nwa_weight()
    integer :: node
    double precision :: denominator_scale
    nlo_decay_nwa_weight = 1d0
    do node = 1, nlo_decay_node_count()
      denominator_scale = decay_dummy_width_ratio()*node_masses(node)**2
      nlo_decay_nwa_weight = nlo_decay_nwa_weight*denominator_scale**2/ &
           (2d0*node_masses(node)* &
            decay_physical_width(&
                 nlo_decay_node_pdg(node), &
                 node == nlo_decay_corrected_node()))
    end do
  end function nlo_decay_nwa_weight


  double precision function norm3(momentum)
    double precision, intent(in) :: momentum(0:3)
    norm3 = sqrt(max(0d0, sum(momentum(1:3)**2)))
  end function norm3


  double precision function pi_value()
    pi_value = 3.141592653589793238462643d0
  end function pi_value


  subroutine check_event_slot(event_slot)
    integer, intent(in) :: event_slot
    if (event_slot < soft_counterevent .or. event_slot > real_event) then
      call fail_kinematics('FKS event slot is out of range')
    end if
  end subroutine check_event_slot


  subroutine require_enabled()
    call initialize_nlo_decay_kinematics()
    if (.not. has_nlo_decay()) then
      call fail_kinematics('no NLO-decay metadata are present')
    end if
  end subroutine require_enabled


  subroutine fail_kinematics(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in nlo_decay_kinematics: '//trim(message)
    stop 1
  end subroutine fail_kinematics

end module nlo_decay_kinematics
