module production_hard_split
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  private

  logical, save :: split_enabled = .false.
  logical, save :: hard_branch = .false.
  double precision, save :: split_scale = 0.1d0
  integer, save :: hard_matrix = 0, hard_configuration = 0
  integer, save :: conditional_branch = 0

  public :: configure_production_hard_split, production_hard_split_enabled
  public :: production_hard_split_scale, set_production_hard_branch
  public :: production_hard_branch_active, production_hard_channel_matrix
  public :: production_hard_channel_config, production_split_factors
  public :: hard_split_partition
  public :: set_production_conditional_branch, production_conditional_branch
  public :: production_conditional_real_branch

contains

  subroutine configure_production_hard_split(enabled, scale)
    logical, intent(in) :: enabled
    double precision, intent(in) :: scale

    if (.not. ieee_is_finite(scale)) then
      write (*, '(a)') 'ERROR: production hard split scale must be finite and positive'
      stop 1
    else if (scale <= 0d0) then
      write (*, '(a)') 'ERROR: production hard split scale must be finite and positive'
      stop 1
    end if
    split_enabled = enabled
    split_scale = scale
    call set_production_hard_branch(.false.)
    conditional_branch = 0
  end subroutine configure_production_hard_split


  logical function production_hard_split_enabled()
    production_hard_split_enabled = split_enabled
  end function production_hard_split_enabled


  double precision function production_hard_split_scale()
    production_hard_split_scale = split_scale
  end function production_hard_split_scale


  subroutine set_production_hard_branch(active, real_matrix, real_config)
    logical, intent(in) :: active
    integer, intent(in), optional :: real_matrix, real_config

    hard_branch = active
    hard_matrix = 0
    hard_configuration = 0
    if (.not. active) return
    if (.not. split_enabled) then
      write (*, '(a)') 'ERROR: a production hard branch requires the hard split'
      stop 1
    end if
    if (present(real_matrix)) hard_matrix = real_matrix
    if (present(real_config)) hard_configuration = real_config
    if (hard_configuration < 1) then
      write (*, '(a)') 'ERROR: a production hard branch has no real diagram channel'
      stop 1
    end if
  end subroutine set_production_hard_branch


  logical function production_hard_branch_active()
    production_hard_branch_active = hard_branch
  end function production_hard_branch_active


  integer function production_hard_channel_matrix()
    production_hard_channel_matrix = hard_matrix
  end function production_hard_channel_matrix


  integer function production_hard_channel_config()
    production_hard_channel_config = hard_configuration
  end function production_hard_channel_config


  subroutine set_production_conditional_branch(branch)
    integer, intent(in) :: branch
    if (branch < 0 .or. branch > 2) stop 'invalid conditional production branch'
    conditional_branch = branch
  end subroutine set_production_conditional_branch


  integer function production_conditional_branch()
    production_conditional_branch = conditional_branch
  end function production_conditional_branch


  logical function production_conditional_real_branch()
    production_conditional_real_branch = conditional_branch == 2
  end function production_conditional_real_branch


  pure subroutine hard_split_partition(invariant, scale, singular, hard)
    ! MadFKS, arXiv:0908.4272, eqs. (6.26)--(6.29).  Return both terms
    ! independently: computing hard=1-singular loses all significant bits
    ! close to a soft/collinear boundary.  The scaled form also avoids
    ! squaring a large invariant/scale ratio.
    double precision, intent(in) :: invariant, scale
    double precision, intent(out) :: singular, hard
    double precision :: ratio_squared

    if (invariant <= scale) then
      ratio_squared = (invariant/scale)**2
      singular = 1d0/(1d0 + ratio_squared)
      hard = ratio_squared/(1d0 + ratio_squared)
    else
      ratio_squared = (scale/invariant)**2
      singular = ratio_squared/(1d0 + ratio_squared)
      hard = 1d0/(1d0 + ratio_squared)
    end if
  end subroutine hard_split_partition


  subroutine production_split_factors(momenta, emitted, sister, shat, &
                                       singular, hard, sister_mass)
    double precision, intent(in) :: momenta(0:, :), shat
    integer, intent(in) :: emitted, sister
    double precision, intent(out) :: singular, hard
    double precision, intent(in), optional :: sister_mass
    double precision :: dot_product, emitted_norm, sister_norm
    double precision :: invariant, energy_difference, direction_difference(3)

    if (.not. ieee_is_finite(shat)) then
      write (*, '(a)') 'ERROR: invalid production hard split kinematics'
      stop 1
    end if
    if (shat <= 0d0 .or. &
        emitted < 1 .or. sister < 1 .or. &
        emitted > size(momenta, 2) .or. sister > size(momenta, 2)) then
      write (*, '(a)') 'ERROR: invalid production hard split kinematics'
      stop 1
    end if
    if (.not. split_enabled) then
      singular = 1d0
      hard = 0d0
      return
    end if

    if (present(sister_mass)) then
      ! The FKS parton is massless.  With the known on-shell sister mass,
      ! this positive form of k_i.k_j avoids subtracting almost equal
      ! energies in a collinear configuration.  It works for incoming
      ! sisters and massive final-state sisters in any common frame.
      emitted_norm = sqrt(sum(momenta(1:3, emitted)**2))
      sister_norm = sqrt(sum(momenta(1:3, sister)**2))
      if (emitted_norm > 0d0 .and. sister_norm > 0d0) then
        direction_difference = momenta(1:3, emitted)/emitted_norm - &
                               momenta(1:3, sister)/sister_norm
        energy_difference = sister_mass**2/ &
             (momenta(0, sister) + sister_norm)
        dot_product = momenta(0, emitted)*(energy_difference + &
             0.5d0*sister_norm*sum(direction_difference**2))
      else
        dot_product = momenta(0, emitted)*momenta(0, sister)
      end if
    else
      dot_product = momenta(0, emitted)*momenta(0, sister) - &
                    sum(momenta(1:3, emitted)*momenta(1:3, sister))
    end if
    invariant = 2d0*dot_product/shat
    if (.not. ieee_is_finite(invariant)) then
      write (*, '(a)') 'ERROR: nonfinite production hard split invariant'
      stop 1
    else if (invariant < 0d0) then
      write (*, '(a,es24.16)') &
           'ERROR: invalid production hard split invariant ', invariant
      stop 1
    end if
    call hard_split_partition(invariant, split_scale, singular, hard)
  end subroutine production_split_factors

end module production_hard_split
