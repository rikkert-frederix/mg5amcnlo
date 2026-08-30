module factorized_block_kinematics
  use process_dimensions, only: nexternal
  use phase_space_kinematics, only: phase_space_lambda
  implicit none
  private

  public :: generate_factorized_nbody
  public :: generate_factorized_nbody_rest
  public :: generate_factorized_decay_tree_rest
  public :: boost_factorized_block_from_rest
  public :: boost_factorized_momentum_from_rest
  public :: factorized_minkowski_square

contains

  subroutine generate_factorized_nbody(parent_momentum, particle_count, &
       masses, x, first_index, momenta, xjac, xpswgt, pass)
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
    parent_mass2 = factorized_minkowski_square(parent_momentum)
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
      rest_child(1) = &
           momentum_length*sqrt(max(0d0, 1d0 - costheta**2))*cos(phi)
      rest_child(2) = &
           momentum_length*sqrt(max(0d0, 1d0 - costheta**2))*sin(phi)
      rest_child(3) = momentum_length*costheta
      rest_remainder(0) = sqrt(remainder_masses(particle + 1)**2 + &
                               momentum_length**2)
      rest_remainder(1:3) = -rest_child(1:3)
      call boost_factorized_momentum_from_rest( &
           rest_child, remainder, parent_mass, child)
      call boost_factorized_momentum_from_rest( &
           rest_remainder, remainder, parent_mass, next_remainder)
      momenta(:, particle) = child
      remainder = next_remainder
      xjac = xjac*4d0*pi
      xpswgt = xpswgt*sqrt(lambda_value)/(8d0*parent_mass**2)
    end do
    momenta(:, particle_count) = remainder
    pass = .true.
  end subroutine generate_factorized_nbody


  subroutine generate_factorized_nbody_rest(parent_mass, particle_count, &
       masses, x, first_index, rest_momenta, xjac, xpswgt, pass)
    double precision, intent(in) :: parent_mass
    integer, intent(in) :: particle_count, first_index
    double precision, intent(in) :: masses(:), x(99)
    double precision, intent(out) :: rest_momenta(0:, :)
    double precision, intent(inout) :: xjac, xpswgt
    logical, intent(out) :: pass
    double precision :: rest_parent(0:3)

    rest_parent = 0d0
    rest_parent(0) = parent_mass
    call generate_factorized_nbody( &
         rest_parent, particle_count, masses, x, first_index, &
         rest_momenta, xjac, xpswgt, pass)
  end subroutine generate_factorized_nbody_rest


  subroutine generate_factorized_decay_tree_rest( &
       parent_mass, particle_count, masses, tree, propagator_masses, &
       propagator_widths, propagator_ids, x, first_index, rest_momenta, &
       xjac, xpswgt, pass)
    double precision, intent(in) :: parent_mass
    integer, intent(in) :: particle_count, first_index
    double precision, intent(in) :: masses(:), x(99)
    integer, intent(in) :: tree(2, -nexternal:-1)
    double precision, intent(in) :: propagator_masses(-nexternal:-1)
    double precision, intent(in) :: propagator_widths(-nexternal:-1)
    integer, intent(in) :: propagator_ids(-nexternal:-1)
    double precision, intent(out) :: rest_momenta(0:, :)
    double precision, intent(inout) :: xjac, xpswgt
    logical, intent(out) :: pass

    integer :: branch, child, internal_count, random_index
    integer :: first_child, second_child
    double precision :: branch_masses(-nexternal:nexternal)
    double precision :: branch_momenta(0:3, -nexternal:nexternal)
    double precision :: lower, upper, invariant, invariant_jacobian
    double precision :: total_mass, lambda_value, momentum_length
    double precision :: costheta, phi
    double precision :: first_rest(0:3), second_rest(0:3)
    double precision, parameter :: pi = 3.141592653589793238462643d0
    double precision, parameter :: tolerance = 1d-10

    pass = .false.
    rest_momenta = 0d0
    if (particle_count < 2 .or. particle_count > size(masses) .or. &
        particle_count > size(rest_momenta, 2) .or. &
        particle_count > nexternal) return
    if (parent_mass <= 0d0 .or. parent_mass + tolerance < &
        sum(masses(1:particle_count))) return

    internal_count = particle_count - 1
    branch_masses = 0d0
    branch_momenta = 0d0
    branch_masses(1:particle_count) = masses(1:particle_count)
    total_mass = sum(masses(1:particle_count))
    random_index = first_index

    ! Internal labels are postordered: -1 is the first completed subtree and
    ! -(n-1) is the fixed-mass decay parent.  Generate every proper
    ! s-channel invariant before constructing any momenta.
    do branch = -1, -(internal_count - 1), -1
      first_child = tree(1, branch)
      second_child = tree(2, branch)
      if (.not. decay_tree_child_is_ready(first_child, branch, &
           particle_count) .or. &
          .not. decay_tree_child_is_ready(second_child, branch, &
           particle_count)) return
      lower = branch_masses(first_child) + branch_masses(second_child)
      upper = parent_mass - total_mass + lower
      if (upper + tolerance < lower) return
      lower = lower**2
      upper = max(lower, upper**2)
      call sample_decay_invariant( &
           x(random_index), lower, upper, propagator_masses(branch), &
           propagator_widths(branch), propagator_ids(branch), invariant, &
           invariant_jacobian)
      branch_masses(branch) = sqrt(max(0d0, invariant))
      xjac = xjac*invariant_jacobian
      total_mass = total_mass + branch_masses(branch) - &
           branch_masses(first_child) - branch_masses(second_child)
      random_index = random_index + 1
    end do

    branch = -internal_count
    first_child = tree(1, branch)
    second_child = tree(2, branch)
    if (.not. decay_tree_child_is_ready(first_child, branch, &
         particle_count) .or. &
        .not. decay_tree_child_is_ready(second_child, branch, &
         particle_count)) return
    if (branch_masses(first_child) + branch_masses(second_child) > &
        parent_mass + tolerance) return
    branch_masses(branch) = parent_mass
    branch_momenta(0, branch) = parent_mass

    ! Split the tree from its root towards the external children.  The same
    ! two-body measure is used as in generate_factorized_nbody; only the
    ! invariant map and clustering topology differ.
    do branch = -internal_count, -1
      first_child = tree(1, branch)
      second_child = tree(2, branch)
      lambda_value = phase_space_lambda(branch_masses(branch)**2, &
           branch_masses(first_child)**2, branch_masses(second_child)**2)
      if (lambda_value < -tolerance*max( &
           1d0, branch_masses(branch)**4)) return
      lambda_value = max(0d0, lambda_value)
      momentum_length = sqrt(lambda_value)/(2d0*branch_masses(branch))
      costheta = 2d0*x(random_index) - 1d0
      phi = 2d0*pi*x(random_index + 1)
      random_index = random_index + 2

      first_rest(0) = sqrt( &
           branch_masses(first_child)**2 + momentum_length**2)
      first_rest(1) = momentum_length*sqrt(max( &
           0d0, 1d0 - costheta**2))*cos(phi)
      first_rest(2) = momentum_length*sqrt(max( &
           0d0, 1d0 - costheta**2))*sin(phi)
      first_rest(3) = momentum_length*costheta
      second_rest(0) = sqrt( &
           branch_masses(second_child)**2 + momentum_length**2)
      second_rest(1:3) = -first_rest(1:3)
      call boost_factorized_momentum_from_rest( &
           first_rest, branch_momenta(:, branch), &
           branch_masses(branch), branch_momenta(:, first_child))
      call boost_factorized_momentum_from_rest( &
           second_rest, branch_momenta(:, branch), &
           branch_masses(branch), branch_momenta(:, second_child))
      xjac = xjac*4d0*pi
      xpswgt = xpswgt*sqrt(lambda_value)/( &
           8d0*branch_masses(branch)**2)
    end do

    do child = 1, particle_count
      rest_momenta(:, child) = branch_momenta(:, child)
    end do
    pass = .true.
  end subroutine generate_factorized_decay_tree_rest


  logical function decay_tree_child_is_ready(label, parent, particle_count)
    integer, intent(in) :: label, parent, particle_count

    decay_tree_child_is_ready = label /= 0
    if (label > 0) then
      decay_tree_child_is_ready = &
           decay_tree_child_is_ready .and. label <= particle_count
    else
      decay_tree_child_is_ready = decay_tree_child_is_ready .and. &
           label >= -nexternal .and. label > parent
    end if
  end function decay_tree_child_is_ready


  subroutine sample_decay_invariant( &
       random_number, lower, upper, mass, width, propagator_id, &
       invariant, jacobian)
    double precision, intent(in) :: random_number, lower, upper
    double precision, intent(in) :: mass, width
    integer, intent(in) :: propagator_id
    double precision, intent(out) :: invariant, jacobian

    if (upper <= lower) then
      invariant = lower
      jacobian = 0d0
    else if (mass > 0d0 .and. width > 0d0) then
      call sample_breit_wigner_invariant( &
           random_number, lower, upper, mass, width, invariant, jacobian)
    else if (propagator_id /= 0 .and. lower > 0d0) then
      call sample_logarithmic_invariant( &
           random_number, lower, upper, invariant, jacobian)
    else
      invariant = lower + (upper - lower)*random_number
      jacobian = upper - lower
    end if
  end subroutine sample_decay_invariant


  subroutine sample_breit_wigner_invariant( &
       random_number, lower, upper, mass, width, invariant, jacobian)
    double precision, intent(in) :: random_number, lower, upper
    double precision, intent(in) :: mass, width
    double precision, intent(out) :: invariant, jacobian
    double precision :: angle, angle_lower, angle_upper, mass_width

    mass_width = mass*width
    angle_lower = atan((lower - mass**2)/mass_width)
    angle_upper = atan((upper - mass**2)/mass_width)
    angle = angle_lower + (angle_upper - angle_lower)*random_number
    invariant = mass**2 + mass_width*tan(angle)
    invariant = min(upper, max(lower, invariant))
    jacobian = mass_width*(angle_upper - angle_lower)/cos(angle)**2
  end subroutine sample_breit_wigner_invariant


  subroutine sample_logarithmic_invariant( &
       random_number, lower, upper, invariant, jacobian)
    double precision, intent(in) :: random_number, lower, upper
    double precision, intent(out) :: invariant, jacobian
    double precision :: logarithmic_range

    logarithmic_range = log(upper/lower)
    invariant = lower*exp(logarithmic_range*random_number)
    invariant = min(upper, max(lower, invariant))
    jacobian = invariant*logarithmic_range
  end subroutine sample_logarithmic_invariant


  subroutine boost_factorized_block_from_rest(rest_momenta, particle_count, &
       parent_momentum, parent_mass, momenta)
    double precision, intent(in) :: rest_momenta(0:, :)
    integer, intent(in) :: particle_count
    double precision, intent(in) :: parent_momentum(0:3), parent_mass
    double precision, intent(out) :: momenta(0:, :)
    integer :: particle

    momenta = 0d0
    do particle = 1, particle_count
      call boost_factorized_momentum_from_rest( &
           rest_momenta(0:3, particle), parent_momentum, parent_mass, &
           momenta(0:3, particle))
    end do
  end subroutine boost_factorized_block_from_rest


  subroutine boost_factorized_momentum_from_rest(rest_momentum, &
       parent_momentum, parent_mass, event_momentum)
    double precision, intent(in) :: rest_momentum(0:3)
    double precision, intent(in) :: parent_momentum(0:3), parent_mass
    double precision, intent(out) :: event_momentum(0:3)
    double precision :: spatial_product, denominator

    spatial_product = dot_product(parent_momentum(1:3), rest_momentum(1:3))
    denominator = parent_mass*(parent_momentum(0) + parent_mass)
    event_momentum(0) = (parent_momentum(0)*rest_momentum(0) + &
                         spatial_product)/parent_mass
    event_momentum(1:3) = rest_momentum(1:3) + parent_momentum(1:3)*( &
         rest_momentum(0)/parent_mass + spatial_product/denominator)
  end subroutine boost_factorized_momentum_from_rest


  double precision function factorized_minkowski_square(momentum)
    double precision, intent(in) :: momentum(0:3)
    factorized_minkowski_square = &
         momentum(0)**2 - sum(momentum(1:3)**2)
  end function factorized_minkowski_square

end module factorized_block_kinematics
