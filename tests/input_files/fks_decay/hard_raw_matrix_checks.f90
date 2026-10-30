program hard_raw_matrix_checks
  use hard_matrix_fixture
  use spin_density_fks_matrices
  use production_hard_split
  implicit none
  double precision :: raw_weight, damped_weight, damping, singular, hard
  complex(kind=8) :: raw_density(2,2), damped_density(2,2)
  real(kind=16) :: reference, xi_q, y_q, invariant_q, scale_q
  integer :: point

  call set_spin_density_fks_collection(.true.)
  do point = 1, 3
    select case(point)
    case(1)
      test_xi=0.2d0; test_y=-0.4d0
    case(2)
      test_xi=0.03d0; test_y=0.8d0
    case(3)
      test_xi=0.8d0; test_y=-0.7d0
    end select
    call evaluate_hard_real_matrix(3,raw_weight)
    call get_spin_density_reduced_matrix(raw_density)
    call sreal(3,test_xi,test_y,damped_weight)
    call get_spin_density_reduced_matrix(damped_density)
    damping=test_xi**2*(1d0-test_y)
    if (abs(damped_weight/(raw_weight*damping)-1d0)>1d-14) stop 1
    if (maxval(abs(damped_density-raw_density*damping))> &
        1d-14*maxval(abs(damped_density))) stop 2
    if (abs(amp_split(1)/(raw_weight*damping)-1d0)>1d-14) stop 3
  end do
  if (limit_calls /= 0) stop 4
  print *, 'PASS: raw and damped native real scalar/density agree away from IR'

  ! Raw real evaluation must keep the actual matrix even in the region
  ! where the ordinary subtracted kernel deliberately substitutes a limit.
  ! Quadruple precision supplies an independent finite hard-weight reference.
  test_xi=1d-100; test_y=0.3d0
  call evaluate_hard_real_matrix(3,raw_weight)
  call get_spin_density_reduced_matrix(raw_density)
  if (limit_calls /= 0) stop 5
  xi_q=real(test_xi,kind=16); y_q=real(test_y,kind=16)
  invariant_q=xi_q*(1q0-y_q)/2q0
  scale_q=0.1q0
  reference=(3.1q0/(xi_q**2*(1q0-y_q)))* &
       invariant_q**2/(scale_q**2+invariant_q**2)
  call hard_split_partition(0.5d0*test_xi*(1d0-test_y),0.1d0,singular,hard)
  if (abs(raw_weight*hard/real(reference,kind=8)-1d0)>2d-14) stop 6
  if (maxval(abs(raw_density*hard - &
      real(reference,kind=8)*fixture_density_shape()))> &
      2d-14*real(reference,kind=8)) stop 7
  call sreal(3,test_xi,test_y,damped_weight)
  if (limit_calls /= 1 .or. damped_weight /= 777d0) stop 8
  print *, 'PASS: raw hard IR weight agrees with quadruple precision; native limit unchanged'
end program hard_raw_matrix_checks

integer function sdm_branch_max_open_size()
  sdm_branch_max_open_size=2
end function

integer function sdm_contribution_component_position(contribution)
  integer, intent(in) :: contribution
  sdm_contribution_component_position=1
end function
