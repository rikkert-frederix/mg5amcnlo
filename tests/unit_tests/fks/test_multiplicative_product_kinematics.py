################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of MadGraph5_aMC@NLO.
#
################################################################################

"""Compiled conservation checks for multiplicative product radiation maps."""

from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile

import tests.unit_tests as unittest

from madgraph import MG5DIR


class TestMultiplicativeProductKinematics(unittest.TestCase):

    def test_initial_and_three_simultaneous_final_state_maps(self):
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is unavailable')

        stubs = r"""
module kin_functions_module
  use iso_fortran_env, only: real64
  implicit none
contains
  real(real64) function dot_impl(first, second)
    real(real64), intent(in) :: first(0:), second(0:)
    dot_impl = first(0)*second(0)-sum(first(1:3)*second(1:3))
  end function dot_impl
end module kin_functions_module

module fks_qcd_splitting
  use iso_fortran_env, only: real64
  implicit none
contains
  subroutine AP_reduced(col1, col2, t, z, gs, result)
    integer, intent(in) :: col1, col2
    real(real64), intent(in) :: t, z, gs
    real(real64), intent(out) :: result(2)
    result = 0._real64
  end subroutine AP_reduced
  subroutine Qterms_reduced_timelike(col1, col2, t, z, gs, result)
    integer, intent(in) :: col1, col2
    real(real64), intent(in) :: t, z, gs
    real(real64), intent(out) :: result(2)
    result = 0._real64
  end subroutine Qterms_reduced_timelike
  subroutine Qterms_reduced_spacelike(col1, col2, t, z, gs, result)
    integer, intent(in) :: col1, col2
    real(real64), intent(in) :: t, z, gs
    real(real64), intent(out) :: result(2)
    result = 0._real64
  end subroutine Qterms_reduced_spacelike
end module fks_qcd_splitting

module fks_model_state_module
  use iso_fortran_env, only: real64
  implicit none
  real(real64), target, save :: coupling = 1._real64
  real(real64), pointer :: strong_coupling => coupling
end module fks_model_state_module

module fks_sij_module
  use iso_fortran_env, only: real64
  implicit none
contains
  subroutine initialize_fks_sij_module(count, incoming, a, b, c, d, local)
    integer, intent(in) :: count, incoming
    real(real64), intent(in) :: a, b, c, d
    logical, intent(in), optional :: local
  end subroutine initialize_fks_sij_module
  subroutine set_fks_sij_partition_state(partners, types, aorg, i, j, &
       rapidity, sqrtshat, shat, phat, masses)
    integer, intent(in) :: partners(:,0:), types(:), i, j
    logical, intent(in) :: aorg(:)
    real(real64), intent(in) :: rapidity, sqrtshat, shat
    real(real64), intent(in) :: phat(0:,0:), masses(:)
  end subroutine set_fks_sij_partition_state
  real(real64) function fks_sij_impl(p, i, j, xi, y)
    real(real64), intent(in) :: p(0:,:), xi, y
    integer, intent(in) :: i, j
    fks_sij_impl = 1._real64
  end function fks_sij_impl
end module fks_sij_module

module decay_chain_kinematics
  use iso_fortran_env, only: real64
  implicit none
contains
  subroutine boost_from_rest(rest_momentum, parent_momentum, parent_mass, &
       lab_momentum)
    real(real64), intent(in) :: rest_momentum(0:3), parent_momentum(0:3)
    real(real64), intent(in) :: parent_mass
    real(real64), intent(out) :: lab_momentum(0:3)
    real(real64) :: spatial_product, denominator
    spatial_product = dot_product(parent_momentum(1:3),rest_momentum(1:3))
    denominator = parent_mass*(parent_momentum(0)+parent_mass)
    lab_momentum(0) = (parent_momentum(0)*rest_momentum(0)+ &
         spatial_product)/parent_mass
    lab_momentum(1:3) = rest_momentum(1:3)+parent_momentum(1:3)*( &
         rest_momentum(0)/parent_mass+spatial_product/denominator)
  end subroutine boost_from_rest
  real(real64) function minkowski_square(momentum)
    real(real64), intent(in) :: momentum(0:3)
    minkowski_square = momentum(0)**2-sum(momentum(1:3)**2)
  end function minkowski_square
end module decay_chain_kinematics

module multiplicative_product
  use iso_fortran_env, only: real64
  implicit none
  type :: product_stage_event
    integer :: stage_id = 0
    real(real64) :: xi = 0._real64
    real(real64) :: y = 0._real64
    real(real64) :: phi = 0._real64
  end type product_stage_event
end module multiplicative_product
"""

        driver = r"""
program check_product_maps
  use iso_fortran_env, only: real64
  use multiplicative_product, only: product_stage_event
  use multiplicative_product_kinematics, only: map_product_final_state, &
       map_product_initial_state
  implicit none
  real(real64), parameter :: mt=173._real64, mw=80.4_real64, mb=4.7_real64
  real(real64), parameter :: energy=250._real64, tolerance=2.e-9_real64
  type(product_stage_event) :: stage
  real(real64) :: p(0:3,9), seed(0:3,9), mass(9)
  real(real64) :: local_p(0:3,5), local_m(5), phat(0:3)
  real(real64) :: sqrtshat, jacobian, q, eb, ew, beta, gamma, ptop
  real(real64) :: initial(0:3), top_before(0:3), antitop_before(0:3)
  integer :: group_sizes(4), group_slots(2,4), real_to_born(5)
  logical :: pass

  p = 0._real64
  mass = 0._real64
  p(0,1)=energy; p(3,1)=energy
  p(0,2)=energy; p(3,2)=-energy
  ptop=sqrt(energy**2-mt**2); beta=ptop/energy; gamma=energy/mt
  q=sqrt((mt**2-(mw+mb)**2)*(mt**2-(mw-mb)**2))/(2._real64*mt)
  eb=sqrt(q**2+mb**2); ew=sqrt(q**2+mw**2)
  p(:,3)=(/gamma*eb,q,0._real64,gamma*beta*eb/)
  p(:,4)=(/gamma*ew,-q,0._real64,gamma*beta*ew/)
  p(:,5)=(/gamma*eb,0._real64,q,-gamma*beta*eb/)
  p(:,6)=(/gamma*ew,0._real64,-q,-gamma*beta*ew/)
  mass(3)=mb; mass(4)=mw; mass(5)=mb; mass(6)=mw
  initial=p(:,1)+p(:,2)
  seed=p

  group_sizes=(/1,1,2,2/)
  group_slots=0
  group_slots(1,1)=1; group_slots(1,2)=2
  group_slots(:,3)=(/3,4/); group_slots(:,4)=(/5,6/)
  real_to_born=(/1,2,3,4,0/)
  stage%stage_id=1; stage%xi=.12_real64; stage%y=-.3_real64
  stage%phi=.4_real64
  call map_product_initial_state(stage,5,1,7,2,group_sizes,group_slots, &
       real_to_born,p,jacobian,local_p,local_m,phat,sqrtshat,pass)
  if (.not.pass .or. jacobian<=0._real64) error stop 1
  call check_close(p(:,1)+p(:,2),sum_columns(p,(/3,4,5,6,7/)))
  call check_close(local_p(:,1)+local_p(:,2), &
       local_p(:,3)+local_p(:,4)+local_p(:,5))
  call check_mass(p(:,1),0._real64); call check_mass(p(:,2),0._real64)
  call check_mass(p(:,3),mb); call check_mass(p(:,4),mw)
  call check_mass(p(:,5),mb); call check_mass(p(:,6),mw)
  call check_mass(p(:,7),0._real64)
  call check_mass(p(:,3)+p(:,4),mt)
  call check_mass(p(:,5)+p(:,6),mt)

  p=seed
  stage%stage_id=1; stage%xi=.08_real64; stage%y=.2_real64
  stage%phi=.7_real64
  call map_product_final_state(stage,5,3,3,7,2,group_sizes,group_slots, &
       real_to_born,p,jacobian,local_p,local_m,phat,sqrtshat,pass)
  if (.not.pass .or. jacobian<=0._real64) error stop 2
  call check_close(sum_columns(p,(/3,4,5,6,7/)),initial)
  call check_mass(p(:,3),mb); call check_mass(p(:,4),mw)
  call check_mass(p(:,5),mb); call check_mass(p(:,6),mw)
  call check_mass(p(:,7),0._real64)
  call check_mass(p(:,3)+p(:,4),mt)
  call check_mass(p(:,5)+p(:,6),mt)
  call check_close(local_p(:,1)+local_p(:,2), &
       local_p(:,3)+local_p(:,4)+local_p(:,5))
  top_before=p(:,3)+p(:,4)

  group_sizes(1:3)=(/2,1,1/); group_slots=0
  group_slots(:,1)=(/3,4/); group_slots(1,2)=3; group_slots(1,3)=4
  real_to_born(1:4)=(/1,2,3,0/)
  stage%stage_id=2; stage%xi=.1_real64; stage%y=.3_real64
  stage%phi=1.1_real64
  call map_product_final_state(stage,4,2,2,8,1,group_sizes(1:3), &
       group_slots(:,1:3),real_to_born(1:4),p,jacobian,local_p(:,1:4), &
       local_m(1:4),phat,sqrtshat,pass)
  if (.not.pass .or. jacobian<=0._real64) error stop 3
  call check_close(local_p(:,1),local_p(:,2)+local_p(:,3)+local_p(:,4))
  call check_mass(local_p(:,1),mt)
  call check_close(p(:,3)+p(:,4)+p(:,8),top_before)
  call check_close(sum_columns(p,(/3,4,5,6,7,8/)),initial)
  call check_mass(p(:,3),mb); call check_mass(p(:,4),mw)
  call check_mass(p(:,8),0._real64)
  antitop_before=p(:,5)+p(:,6)

  group_slots=0; group_slots(:,1)=(/5,6/)
  group_slots(1,2)=5; group_slots(1,3)=6
  stage%stage_id=3; stage%xi=.07_real64; stage%y=-.25_real64
  stage%phi=2.2_real64
  call map_product_final_state(stage,4,2,2,9,1,group_sizes(1:3), &
       group_slots(:,1:3),real_to_born(1:4),p,jacobian,local_p(:,1:4), &
       local_m(1:4),phat,sqrtshat,pass)
  if (.not.pass .or. jacobian<=0._real64) error stop 4
  call check_close(local_p(:,1),local_p(:,2)+local_p(:,3)+local_p(:,4))
  call check_mass(local_p(:,1),mt)
  call check_close(p(:,5)+p(:,6)+p(:,9),antitop_before)
  call check_close(sum_columns(p,(/3,4,5,6,7,8,9/)),initial)
  call check_mass(p(:,5),mb); call check_mass(p(:,6),mw)
  call check_mass(p(:,9),0._real64)
contains
  function sum_columns(values, indices) result(total)
    real(real64), intent(in) :: values(0:3,9)
    integer, intent(in) :: indices(:)
    real(real64) :: total(0:3)
    integer :: index
    total=0._real64
    do index=1,size(indices)
      total=total+values(:,indices(index))
    end do
  end function sum_columns
  subroutine check_close(first,second)
    real(real64), intent(in) :: first(0:3),second(0:3)
    if (maxval(abs(first-second))>tolerance*energy) error stop 5
  end subroutine check_close
  subroutine check_mass(momentum,expected)
    real(real64), intent(in) :: momentum(0:3),expected
    real(real64) :: square
    square=momentum(0)**2-sum(momentum(1:3)**2)
    if (abs(square-expected**2)>tolerance*energy**2) error stop 6
  end subroutine check_mass
end program check_product_maps
"""

        boost = os.path.join(
            MG5DIR, 'Template', 'fNLO', 'SubProcesses', 'boostwdir2.f90')
        kinematics = os.path.join(
            MG5DIR, 'Template', 'fNLO', 'SubProcesses',
            'multiplicative_product_kinematics.f90')
        with tempfile.TemporaryDirectory() as workdir:
            stubs_path = os.path.join(workdir, 'stubs.f90')
            driver_path = os.path.join(workdir, 'driver.f90')
            executable = os.path.join(workdir, 'check_product_maps')
            with open(stubs_path, 'w') as stream:
                stream.write(stubs)
            with open(driver_path, 'w') as stream:
                stream.write(driver)
            subprocess.check_call([
                compiler, '-std=f2008', '-Wall', '-Wextra',
                '-Wno-unused-dummy-argument', '-Wno-compare-reals',
                '-fcheck=all', '-J', workdir, stubs_path, boost,
                kinematics, driver_path, '-o', executable], cwd=workdir)
            subprocess.check_call([executable], cwd=workdir)


if __name__ == '__main__':
    unittest.main()
