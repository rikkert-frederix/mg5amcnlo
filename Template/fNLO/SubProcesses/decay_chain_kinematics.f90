module decay_chain_kinematics
  use process_dimensions, only: nexternal, nincoming, validate_process_dimensions
  use phase_space_kinematics, only: phase_space_lambda
  use decay_chain_metadata, only: has_decay_chains, decay_node_count, &
       decay_leaf_count, decay_random_dimension, born_context, &
       context_for_fks, context_core_count, &
       core_target_kind, core_target_id, core_leg_pdg, &
       direct_leg_target, decay_node_target, leaf_visible_leg, &
       node_pdg, node_child_count, node_child_kind, node_child_id, &
       decay_leaf_child, decay_node_child, leaf_pdg, &
       visible_color_pair
  use decay_chain_parameters, only: decay_dummy_width_ratio, &
       decay_physical_width
  implicit none
  private

  logical, save :: initialized = .false.
  double precision, allocatable, save :: node_masses(:)
  double precision, allocatable, save :: leaf_masses(:)
  double precision, allocatable, save :: core_born_storage(:, :)
  double precision, allocatable, save :: core_event_storage(:, :, :)

  public :: initialize_decay_chain_kinematics
  public :: minimum_core_final_mass, core_mass
  public :: generate_core_born_and_decays, expand_real_decay_momenta
  public :: store_core_event_momenta, get_core_event_momenta
  public :: get_core_born_momenta, get_core_mass_buffer
  public :: active_core_count, fks_leg_mass
  public :: map_core_color_pair
  public :: contract_visible_momenta

  interface
    double precision function get_mass_from_id(id)
      integer, intent(in) :: id
    end function get_mass_from_id
  end interface

contains

  subroutine initialize_decay_chain_kinematics()
    integer :: node, leaf

    if (initialized) return
    call validate_process_dimensions()
    if (.not. has_decay_chains()) then
      initialized = .true.
      return
    end if

    allocate(node_masses(decay_node_count()))
    allocate(leaf_masses(decay_leaf_count()))
    allocate(core_born_storage(0:3, nexternal - 1))
    allocate(core_event_storage(0:3, nexternal, 0:3))
    core_born_storage = 0d0
    core_event_storage = 0d0

    do node = 1, decay_node_count()
      node_masses(node) = abs(get_mass_from_id(node_pdg(node)))
      if (node_masses(node) <= 0d0) then
        call fail_kinematics('a forced decay parent has zero model mass')
      end if
    end do
    do leaf = 1, decay_leaf_count()
      leaf_masses(leaf) = abs(get_mass_from_id(leaf_pdg(leaf)))
    end do
    initialized = .true.
  end subroutine initialize_decay_chain_kinematics


  integer function production_random_dimension()
    integer :: final_count
    call require_enabled()
    final_count = context_core_count(born_context()) - nincoming
    if (final_count == 1) then
      production_random_dimension = 0
    else
      production_random_dimension = 3*final_count - 4
    end if
  end function production_random_dimension


  integer function decay_variable_start()
    call require_enabled()
    decay_variable_start = production_random_dimension() + 1
  end function decay_variable_start


  double precision function minimum_core_final_mass()
    integer :: context, leg
    call require_enabled()
    context = born_context()
    minimum_core_final_mass = 0d0
    do leg = nincoming + 1, context_core_count(context)
      minimum_core_final_mass = minimum_core_final_mass + &
           core_mass(context, leg)
    end do
  end function minimum_core_final_mass


  double precision function core_mass(context, leg)
    integer, intent(in) :: context, leg
    integer :: target
    call require_enabled()
    if (core_target_kind(context, leg) == decay_node_target) then
      target = core_target_id(context, leg)
      core_mass = node_masses(target)
    else
      core_mass = abs(get_mass_from_id(core_leg_pdg(context, leg)))
    end if
  end function core_mass


  subroutine generate_core_born_and_decays(x, shat, sqrtshat, xjac, &
                                            xpswgt, visible_momenta, pass)
    double precision, intent(in) :: x(99), shat, sqrtshat
    double precision, intent(inout) :: xjac, xpswgt
    double precision, intent(out) :: visible_momenta(0:3, nexternal - 1)
    logical, intent(out) :: pass

    integer :: context, core_count, final_count, leg, index
    double precision :: parent(0:3)
    double precision :: final_masses(nexternal)
    double precision :: final_momenta(0:3, nexternal)

    call require_enabled()
    context = born_context()
    core_count = context_core_count(context)
    final_count = core_count - nincoming
    parent = 0d0
    parent(0) = sqrtshat
    do leg = 1, final_count
      final_masses(leg) = core_mass(context, nincoming + leg)
    end do

    call generate_nbody(parent, final_count, final_masses, x, 1, &
                        final_momenta, xjac, xpswgt, pass)
    if (.not. pass) return

    core_born_storage = 0d0
    call fill_incoming_momenta(shat, sqrtshat, context, &
                               core_born_storage, pass)
    if (.not. pass) return
    do leg = 1, final_count
      core_born_storage(:, nincoming + leg) = final_momenta(:, leg)
    end do

    visible_momenta = 0d0
    index = decay_variable_start()
    call expand_context(context, core_born_storage, x, index, &
                        visible_momenta, xjac, xpswgt, .true., pass)
    if (.not. pass) return
    if (index /= decay_variable_start() + decay_random_dimension()) then
      call fail_kinematics('decay random-variable accounting is inconsistent')
    end if
    xpswgt = xpswgt*decay_nwa_weight()
  end subroutine generate_core_born_and_decays


  subroutine expand_real_decay_momenta(configuration, x, core_momenta, &
                                       visible_momenta, pass)
    integer, intent(in) :: configuration
    double precision, intent(in) :: x(99)
    double precision, intent(in) :: core_momenta(0:3, nexternal)
    double precision, intent(out) :: visible_momenta(0:3, nexternal)
    logical, intent(out) :: pass
    integer :: context, index
    double precision :: unused_jacobian, unused_weight

    call require_enabled()
    context = context_for_fks(configuration)
    visible_momenta = 0d0
    index = decay_variable_start()
    unused_jacobian = 1d0
    unused_weight = 1d0
    call expand_context(context, core_momenta, x, index, visible_momenta, &
                        unused_jacobian, unused_weight, .false., pass)
    if (.not. pass) return
    if (index /= decay_variable_start() + decay_random_dimension()) then
      call fail_kinematics('real decay random-variable accounting is inconsistent')
    end if
  end subroutine expand_real_decay_momenta


  subroutine expand_context(context, core_momenta, x, index, &
                            visible_momenta, xjac, xpswgt, keep_weight, pass)
    integer, intent(in) :: context
    double precision, intent(in) :: core_momenta(0:, :), x(99)
    integer, intent(inout) :: index
    double precision, intent(out) :: visible_momenta(0:, :)
    double precision, intent(inout) :: xjac, xpswgt
    logical, intent(in) :: keep_weight
    logical, intent(out) :: pass

    integer :: leg, target
    pass = .true.
    do leg = 1, context_core_count(context)
      target = core_target_id(context, leg)
      if (core_target_kind(context, leg) == direct_leg_target) then
        visible_momenta(:, target) = core_momenta(:, leg)
      else
        call expand_node(context, target, core_momenta(:, leg), x, index, &
                         visible_momenta, xjac, xpswgt, keep_weight, pass)
        if (.not. pass) return
      end if
    end do
  end subroutine expand_context


  recursive subroutine expand_node(context, node, parent_momentum, x, &
                                   index, visible_momenta, xjac, xpswgt, &
                                   keep_weight, pass)
    integer, intent(in) :: context, node
    double precision, intent(in) :: parent_momentum(0:3), x(99)
    integer, intent(inout) :: index
    double precision, intent(inout) :: visible_momenta(0:, :)
    double precision, intent(inout) :: xjac, xpswgt
    logical, intent(in) :: keep_weight
    logical, intent(out) :: pass

    integer :: child_count, child, identifier, child_kind
    double precision :: child_masses(nexternal)
    double precision :: child_momenta(0:3, nexternal)
    double precision :: local_jacobian, local_weight

    child_count = node_child_count(node)
    do child = 1, child_count
      identifier = node_child_id(node, child)
      child_kind = node_child_kind(node, child)
      if (child_kind == decay_node_child) then
        child_masses(child) = node_masses(identifier)
      else
        child_masses(child) = leaf_masses(identifier)
      end if
    end do

    local_jacobian = 1d0
    local_weight = 1d0
    call generate_nbody(parent_momentum, child_count, child_masses, x, &
                        index, child_momenta, local_jacobian, local_weight, &
                        pass)
    if (.not. pass) return
    index = index + 3*child_count - 4
    if (keep_weight) then
      xjac = xjac*local_jacobian
      xpswgt = xpswgt*local_weight
    end if

    do child = 1, child_count
      identifier = node_child_id(node, child)
      child_kind = node_child_kind(node, child)
      if (child_kind == decay_node_child) then
        call expand_node(context, identifier, child_momenta(:, child), x, &
                         index, visible_momenta, xjac, xpswgt, keep_weight, &
                         pass)
        if (.not. pass) return
      else
        visible_momenta(:, leaf_visible_leg(context, identifier)) = &
             child_momenta(:, child)
      end if
    end do
  end subroutine expand_node


  subroutine generate_nbody(parent_momentum, particle_count, masses, x, &
                            first_index, momenta, xjac, xpswgt, pass)
    double precision, intent(in) :: parent_momentum(0:3)
    integer, intent(in) :: particle_count, first_index
    double precision, intent(in) :: masses(:), x(99)
    double precision, intent(out) :: momenta(0:, :)
    double precision, intent(inout) :: xjac, xpswgt
    logical, intent(out) :: pass

    integer :: particle, random_index
    double precision :: parent_mass2, parent_mass, lower, upper, invariant
    double precision :: remainder_masses(nexternal)
    double precision :: remainder(0:3), child(0:3), next_remainder(0:3)
    double precision :: costheta, phi, lambda_value, momentum_length
    double precision :: rest_child(0:3), rest_remainder(0:3)
    double precision, parameter :: pi = 3.141592653589793238462643d0
    double precision, parameter :: tolerance = 1d-10

    pass = .false.
    momenta = 0d0
    if (particle_count < 1 .or. particle_count > size(masses) .or. &
        particle_count > size(momenta, 2)) return
    parent_mass2 = minkowski_square(parent_momentum)
    if (parent_mass2 <= 0d0) return
    parent_mass = sqrt(parent_mass2)
    if (parent_mass + tolerance < sum(masses(1:particle_count))) return

    if (particle_count == 1) then
      if (abs(parent_mass - masses(1)) > &
          tolerance*max(1d0, parent_mass)) return
      momenta(:, 1) = parent_momentum
      pass = .true.
      return
    end if

    remainder_masses = 0d0
    remainder_masses(1) = parent_mass
    remainder_masses(particle_count) = masses(particle_count)
    random_index = first_index
    do particle = particle_count - 1, 2, -1
      lower = masses(particle) + remainder_masses(particle + 1)
      upper = parent_mass - sum(masses(1:particle - 1))
      if (upper + tolerance < lower) return
      lower = lower**2
      upper = max(lower, upper**2)
      invariant = lower + (upper - lower)*x(random_index)
      remainder_masses(particle) = sqrt(max(0d0, invariant))
      xjac = xjac*(upper - lower)
      random_index = random_index + 1
    end do

    remainder = parent_momentum
    do particle = 1, particle_count - 1
      parent_mass = remainder_masses(particle)
      lambda_value = phase_space_lambda(parent_mass**2, &
           masses(particle)**2, remainder_masses(particle + 1)**2)
      if (lambda_value < -tolerance*max(1d0, parent_mass**4)) return
      lambda_value = max(0d0, lambda_value)
      momentum_length = sqrt(lambda_value)/(2d0*parent_mass)
      costheta = 2d0*x(random_index) - 1d0
      phi = 2d0*pi*x(random_index + 1)
      random_index = random_index + 2
      rest_child(0) = sqrt(masses(particle)**2 + momentum_length**2)
      rest_child(1) = momentum_length*sqrt(max(0d0, 1d0 - costheta**2))*cos(phi)
      rest_child(2) = momentum_length*sqrt(max(0d0, 1d0 - costheta**2))*sin(phi)
      rest_child(3) = momentum_length*costheta
      rest_remainder(0) = sqrt(remainder_masses(particle + 1)**2 + &
                               momentum_length**2)
      rest_remainder(1:3) = -rest_child(1:3)
      call boost_from_rest(rest_child, remainder, parent_mass, child)
      call boost_from_rest(rest_remainder, remainder, parent_mass, &
                           next_remainder)
      momenta(:, particle) = child
      remainder = next_remainder
      xjac = xjac*4d0*pi
      xpswgt = xpswgt*sqrt(lambda_value)/(8d0*parent_mass**2)
    end do
    momenta(:, particle_count) = remainder
    pass = .true.
  end subroutine generate_nbody


  subroutine boost_from_rest(rest_momentum, parent_momentum, parent_mass, &
                             lab_momentum)
    double precision, intent(in) :: rest_momentum(0:3)
    double precision, intent(in) :: parent_momentum(0:3), parent_mass
    double precision, intent(out) :: lab_momentum(0:3)
    double precision :: spatial_product, denominator

    spatial_product = dot_product(parent_momentum(1:3), rest_momentum(1:3))
    denominator = parent_mass*(parent_momentum(0) + parent_mass)
    lab_momentum(0) = (parent_momentum(0)*rest_momentum(0) + &
                       spatial_product)/parent_mass
    lab_momentum(1:3) = rest_momentum(1:3) + parent_momentum(1:3)*( &
         rest_momentum(0)/parent_mass + spatial_product/denominator)
  end subroutine boost_from_rest


  subroutine fill_incoming_momenta(shat, sqrtshat, context, momenta, pass)
    double precision, intent(in) :: shat, sqrtshat
    integer, intent(in) :: context
    double precision, intent(inout) :: momenta(0:, :)
    logical, intent(out) :: pass
    double precision :: mass1, mass2, lambda_value, momentum_length

    pass = .false.
    if (nincoming == 1) then
      momenta(:, 1) = 0d0
      momenta(0, 1) = sqrtshat
      pass = .true.
      return
    end if
    mass1 = core_mass(context, 1)
    mass2 = core_mass(context, 2)
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
  end subroutine fill_incoming_momenta


  subroutine store_core_event_momenta(event_slot, momenta)
    integer, intent(in) :: event_slot
    double precision, intent(in) :: momenta(0:, :)
    integer :: count
    call require_enabled()
    if (event_slot < 0 .or. event_slot > 3) then
      call fail_kinematics('event slot is out of range')
    end if
    count = min(size(momenta, 2), nexternal)
    core_event_storage(:, :, event_slot) = 0d0
    core_event_storage(:, 1:count, event_slot) = momenta(:, 1:count)
  end subroutine store_core_event_momenta


  subroutine get_core_event_momenta(event_slot, momenta)
    integer, intent(in) :: event_slot
    double precision, intent(out) :: momenta(0:3, nexternal)
    call require_enabled()
    if (event_slot < 0 .or. event_slot > 3) then
      call fail_kinematics('event slot is out of range')
    end if
    momenta = core_event_storage(:, :, event_slot)
  end subroutine get_core_event_momenta


  subroutine get_core_born_momenta(momenta)
    double precision, intent(out) :: momenta(0:3, nexternal - 1)
    call require_enabled()
    momenta = core_born_storage
  end subroutine get_core_born_momenta


  subroutine get_core_mass_buffer(context, masses)
    integer, intent(in) :: context
    double precision, intent(out) :: masses(nexternal)
    integer :: leg
    call require_enabled()
    masses = 0d0
    do leg = 1, context_core_count(context)
      masses(leg) = core_mass(context, leg)
    end do
  end subroutine get_core_mass_buffer


  integer function active_core_count(configuration)
    integer, intent(in) :: configuration
    call require_enabled()
    active_core_count = context_core_count(context_for_fks(configuration))
  end function active_core_count


  double precision function fks_leg_mass(configuration, leg)
    integer, intent(in) :: configuration, leg
    call require_enabled()
    fks_leg_mass = core_mass(context_for_fks(configuration), leg)
  end function fks_leg_mass


  subroutine map_core_color_pair(core_first, core_second, visible_first, &
                                 visible_second)
    integer, intent(in) :: core_first, core_second
    integer, intent(out) :: visible_first, visible_second
    call require_enabled()
    call visible_color_pair(core_first, core_second, visible_first, &
                            visible_second)
  end subroutine map_core_color_pair


  subroutine contract_visible_momenta(context, visible_momenta, &
                                      core_momenta)
    integer, intent(in) :: context
    double precision, intent(in) :: visible_momenta(0:3, nexternal)
    double precision, intent(out) :: core_momenta(0:3, nexternal)
    integer :: leg, target

    call require_enabled()
    core_momenta = 0d0
    do leg = 1, context_core_count(context)
      target = core_target_id(context, leg)
      if (core_target_kind(context, leg) == direct_leg_target) then
        core_momenta(:, leg) = visible_momenta(:, target)
      else
        call sum_decay_node_momenta(context, target, visible_momenta, &
                                    core_momenta(:, leg))
      end if
    end do
  end subroutine contract_visible_momenta


  recursive subroutine sum_decay_node_momenta(context, node, &
                                               visible_momenta, momentum)
    integer, intent(in) :: context, node
    double precision, intent(in) :: visible_momenta(0:3, nexternal)
    double precision, intent(out) :: momentum(0:3)
    double precision :: child_momentum(0:3)
    integer :: child, identifier

    momentum = 0d0
    do child = 1, node_child_count(node)
      identifier = node_child_id(node, child)
      if (node_child_kind(node, child) == decay_node_child) then
        call sum_decay_node_momenta(context, identifier, visible_momenta, &
                                    child_momentum)
      else
        child_momentum = &
             visible_momenta(:, leaf_visible_leg(context, identifier))
      end if
      momentum = momentum + child_momentum
    end do
  end subroutine sum_decay_node_momenta


  double precision function decay_nwa_weight()
    integer :: node
    double precision :: denominator_scale, physical_width
    call require_enabled()
    decay_nwa_weight = 1d0
    do node = 1, decay_node_count()
      denominator_scale = decay_dummy_width_ratio()*node_masses(node)**2
      physical_width = decay_physical_width(node_pdg(node))
      decay_nwa_weight = decay_nwa_weight*denominator_scale**2/ &
           (2d0*node_masses(node)*physical_width)
    end do
  end function decay_nwa_weight


  double precision function minkowski_square(momentum)
    double precision, intent(in) :: momentum(0:3)
    minkowski_square = momentum(0)**2 - sum(momentum(1:3)**2)
  end function minkowski_square


  subroutine require_enabled()
    if (.not. initialized) call initialize_decay_chain_kinematics()
    if (.not. has_decay_chains()) then
      call fail_kinematics('no decay-chain metadata are present')
    end if
  end subroutine require_enabled


  subroutine fail_kinematics(message)
    character(len=*), intent(in) :: message
    write (*, '(a)') 'ERROR in decay_chain_kinematics: '//trim(message)
    stop 1
  end subroutine fail_kinematics

end module decay_chain_kinematics
