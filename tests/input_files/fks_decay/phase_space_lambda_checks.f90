program phase_space_lambda_checks
  use phase_space_kinematics, only: phase_space_lambda
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  double precision :: value, scale, parent, daughter, expected
  integer :: i

  value = phase_space_lambda(0d0, 0d0, 0d0)
  if (.not. ieee_is_finite(value) .or. value /= 0d0) error stop 1
  if (phase_space_lambda(100d0, 0d0, 0d0) /= 10000d0) error stop 2
  if (phase_space_lambda(36d0, 4d0, 9d0) /= 385d0) error stop 3
  if (phase_space_lambda(36d0, 9d0, 4d0) /= 385d0) error stop 4
  ! The expanded polynomial subtracts O(S^2) terms to obtain O(delta^2)
  ! here. A positive result can still have percent-level relative error.
  parent = 172.5d0**2
  daughter = parent*(1d0-1d-7)
  expected = (parent-daughter)**2
  value = phase_space_lambda(parent,0d0,daughter)
  if (abs(value/expected-1d0)>1d-13) error stop 8
  value = phase_space_lambda(parent,daughter,0d0)
  if (abs(value/expected-1d0)>1d-13) error stop 9
  do i = -6, 6, 6
    scale = 10d0**i
    value = phase_space_lambda(25d0*scale, 4d0*scale, 9d0*scale)
    if (abs(value)/scale**2 > 1d-12) error stop 5
    value = phase_space_lambda((25d0-1d-9)*scale, 4d0*scale, 9d0*scale)
    if (value /= 0d0) error stop 6
    ! A real violation must not be mistaken for roundoff in any unit system.
    value = phase_space_lambda(24d0*scale, 4d0*scale, 9d0*scale)
    if (abs(value/scale**2+23d0) > 1d-10) error stop 7
  end do
  print *, 'PASS: zero invariant and dimensionless threshold checks'
end program phase_space_lambda_checks
