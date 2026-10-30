program production_hard_split_checks
  use production_hard_split
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  implicit none
  integer, parameter :: n = 96
  double precision :: p(0:3, 3), boosted(0:3, 3)
  double precision :: singular, hard, singular_boosted, hard_boosted
  double precision :: nodes(n), weights(n), x, z, jacobian, real_value
  double precision :: eta, ch, sh, local_weights(3), channel_partition
  complex(kind=8) :: a(2, 2), b(2, 2), c(2, 2), real_density(2, 2)
  complex(kind=8) :: singular_integral(2, 2), hard_integral(2, 2), expected(2, 2)
  complex(kind=8) :: point_singular(2, 2), point_hard(2, 2)
  integer :: i, j, k, channel
  character(len=32) :: mode

  call get_command_argument(1, mode)
  if (mode == 'bad_scale') then
    call configure_production_hard_split(.true., 0d0)
    stop 99
  else if (mode == 'nan_scale') then
    call configure_production_hard_split(.true., ieee_value(1d0, ieee_quiet_nan))
    stop 99
  else if (mode == 'bad_channel') then
    call configure_production_hard_split(.true., 0.1d0)
    call set_production_hard_branch(.true.)
    stop 99
  end if

  if (production_hard_split_enabled() .or. production_hard_branch_active()) stop 1
  call configure_production_hard_split(.true., 0.1d0)
  call set_production_hard_branch(.true., real_matrix=3, real_config=2)
  if (.not. production_hard_branch_active() .or. &
      production_hard_channel_matrix() /= 3 .or. &
      production_hard_channel_config() /= 2) stop 2
  call set_production_hard_branch(.false.)

  ! An ISR collinear limit has exactly zero invariant.  At a merely tiny
  ! invariant the hard term must remain nonzero even when singular rounds
  ! to one, which a numerical subtraction 1-singular cannot achieve.
  call hard_split_partition(0d0, 0.1d0, singular, hard)
  if (singular /= 1d0 .or. hard /= 0d0) stop 3
  call hard_split_partition(1d-120, 0.1d0, singular, hard)
  if (hard <= 0d0 .or. abs(hard/1d-238-1d0) > 1d-14) stop 4
  call hard_split_partition(1d200, 1d-100, singular, hard)
  if (singular /= 0d0 .or. hard /= 1d0) stop 5

  ! A massive sister has no collinear singularity: its hard piece need not
  ! disappear at y=1.  Lorentz boosting the same core must not change it.
  p(:, 1) = [50d0, 0d0, 0d0, 50d0]
  p(:, 2) = [sqrt(100d0**2+175d0**2), 0d0, 0d0, 100d0]
  p(:, 3) = [500d0, 0d0, 0d0, -500d0]
  call production_split_factors(p, 1, 2, 1d6, singular, hard, 175d0)
  if (hard <= 0d0 .or. singular >= 1d0) stop 6
  eta = 1.3d0
  ch = cosh(eta)
  sh = sinh(eta)
  do i = 1, 3
    boosted(0, i) = ch*p(0, i) + sh*p(1, i)
    boosted(1, i) = sh*p(0, i) + ch*p(1, i)
    boosted(2:3, i) = p(2:3, i)
  end do
  call production_split_factors( &
       boosted, 1, 2, 1d6, singular_boosted, hard_boosted, 175d0)
  if (abs(hard_boosted-hard) > 1d-14 .or. &
      abs(singular_boosted-singular) > 1d-14) stop 7
  p(:, 1) = 0d0
  call production_split_factors(p, 1, 2, 1d6, singular, hard, 175d0)
  if (hard /= 0d0 .or. singular /= 1d0) stop 8
  write (*, '(a)') 'PASS: IR limits, massive sister and Lorentz invariance'

  ! Integrate a genuinely subtracted two-variable model.  The singular
  ! contribution uses x=u,z=v; the hard contribution uses a different
  ! direct phase-space proposal x=u**2,z=v**2 with its own Jacobian.
  ! The matrix-valued real density contains 1/(x*z), with an unchanged
  ! local counterterm A/(x*z).  Their known integrated sum is finite.
  a = reshape([(2d0,0d0), (0d0,-1d0), (0d0,1d0), (3d0,0d0)], [2,2])
  b = reshape([(0.3d0,0d0), (0.2d0,-0.1d0), &
               (0.2d0,0.1d0), (0.4d0,0d0)], [2,2])
  c = reshape([(0.7d0,0d0), (0d0,0.1d0), &
               (0d0,-0.1d0), (0.8d0,0d0)], [2,2])
  expected = 6d0*a + 0.5d0*(b+c)
  singular_integral = (0d0,0d0)
  hard_integral = (0d0,0d0)
  call gauss_legendre(nodes, weights)
  do i = 1, n
    do j = 1, n
      x = nodes(i)
      z = nodes(j)
      call hard_split_partition(0.5d0*x*z, 0.1d0, singular, hard)
      real_value = 1d0/(x*z) + 3d0 + 2d0*x + 4d0*z
      real_density = real_value*a + x*b + z*c
      point_singular = singular*real_density - a/(x*z)
      point_hard = hard*real_density
      if (maxval(abs(point_singular+point_hard - &
          ((3d0+2d0*x+4d0*z)*a+x*b+z*c))) > &
          5d-15*maxval(abs(real_density))) stop 9
      singular_integral = singular_integral + weights(i)*weights(j)*point_singular

      x = nodes(i)**2
      z = nodes(j)**2
      jacobian = 4d0*nodes(i)*nodes(j)
      call hard_split_partition(0.5d0*x*z, 0.1d0, singular, hard)
      real_value = 1d0/(x*z) + 3d0 + 2d0*x + 4d0*z
      real_density = real_value*a + x*b + z*c
      local_weights = [1d0, x, z]
      do channel = 1, 3
        channel_partition = local_weights(channel)/sum(local_weights)
        hard_integral = hard_integral + &
             weights(i)*weights(j)*jacobian*channel_partition*hard*real_density
      end do
    end do
  end do
  if (maxval(abs(singular_integral+hard_integral-expected)) > 2d-11) then
    write (*, *) 'Integral mismatch:', singular_integral+hard_integral-expected
    stop 10
  end if
  write (*, '(a)') 'PASS: subtracted scalar/density closure with independent proposals'

  call configure_production_hard_split(.false., 0.1d0)
  call production_split_factors(p, 1, 2, 1d6, singular, hard, 175d0)
  if (singular /= 1d0 .or. hard /= 0d0 .or. production_hard_branch_active()) stop 11
  write (*, '(a)') 'PASS: disabled split'

contains

  subroutine gauss_legendre(nodes, weights)
    double precision, intent(out) :: nodes(:), weights(:)
    double precision :: root, previous, polynomial, p_previous, p_before
    double precision :: derivative, pi
    integer :: order, first, iteration
    order = size(nodes)
    pi = acos(-1d0)
    do first = 1, (order+1)/2
      root = cos(pi*(first-0.25d0)/(order+0.5d0))
      do iteration = 1, 100
        polynomial = 1d0
        p_previous = 0d0
        do k = 1, order
          p_before = p_previous
          p_previous = polynomial
          polynomial = ((2d0*k-1d0)*root*p_previous - &
                        (k-1d0)*p_before)/k
        end do
        derivative = order*(root*polynomial-p_previous)/(root**2-1d0)
        previous = root
        root = previous-polynomial/derivative
        if (abs(root-previous) <= 4d0*epsilon(1d0)) exit
      end do
      nodes(first) = (1d0-root)/2d0
      nodes(order+1-first) = (1d0+root)/2d0
      weights(first) = 1d0/((1d0-root**2)*derivative**2)
      weights(order+1-first) = weights(first)
    end do
  end subroutine gauss_legendre

end program production_hard_split_checks
