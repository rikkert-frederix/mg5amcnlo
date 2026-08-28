module factorized_block_kinematics
  use process_dimensions, only: nexternal
  use phase_space_kinematics, only: phase_space_lambda
  implicit none
  private

  public :: generate_factorized_nbody
  public :: generate_factorized_nbody_rest
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
    call generate_factorized_nbody(rest_parent, particle_count, masses, x, &
         first_index, rest_momenta, xjac, xpswgt, pass)
  end subroutine generate_factorized_nbody_rest


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
