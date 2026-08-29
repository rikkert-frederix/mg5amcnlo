from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class MultiplicativeKinematicsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.compiler = shutil.which('gfortran')
        cls.repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        cls.sources = [os.path.join(
            cls.repository, 'Template', 'fNLO', 'SubProcesses', name)
            for name in (
                'factorized_phase_space.f90',
                'factorized_block_kinematics.f90',
                'multiplicative_kinematics.f90')]

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_nested_blocks_are_reboosted_for_each_tuple(self):
        stubs = r'''
module process_dimensions
  implicit none
  integer, parameter :: nexternal = 5
contains
  subroutine validate_process_dimensions()
  end subroutine validate_process_dimensions
end module process_dimensions

module fnlo_process_common
  implicit none
  integer, parameter :: soft_counterevent = 0
  integer, parameter :: real_event = 3
end module fnlo_process_common

module phase_space_kinematics
  implicit none
contains
  double precision function phase_space_lambda(x, y, z)
    double precision, intent(in) :: x, y, z
    phase_space_lambda = x*x + y*y + z*z - 2d0*(x*y + x*z + y*z)
  end function phase_space_lambda
end module phase_space_kinematics

module kin_functions_module
  implicit none
contains
  double precision function et_impl(p)
    double precision, intent(in) :: p(0:3)
    double precision :: pt, spatial
    pt = sqrt(p(1)**2 + p(2)**2)
    spatial = sqrt(pt**2 + p(3)**2)
    if (spatial > 0d0) then
      et_impl = p(0)*pt/spatial
    else
      et_impl = 0d0
    end if
  end function et_impl
end module kin_functions_module
'''
        program = r'''
program test_nested_tuple_boosts
  use factorized_phase_space
  use multiplicative_kinematics
  implicit none
  integer :: slots(0:5), pdgs(5), kinds(5), targets(5)
  integer :: previous_slots(0:5)
  integer(kind=8) :: root_revision, first_revision, second_revision
  integer :: visible_pdgs(5), visible_origins(5), visible_count
  logical :: finals(5), pass, available
  type(factorized_measure_state) :: measure, fetched
  double precision :: root(0:3, 5), first(0:3, 5), second(0:3, 5)
  double precision :: realized_first(0:3, 3), realized_second(0:3, 3)
  double precision :: previous_second_parent(0:3)
  double precision :: visible(0:3, 5)
  double precision :: sum_et, sum_mt

  call reset_factorized_phase_space()
  measure = factorized_measure_state()
  measure%jacobian = 2d0
  measure%phase_space_weight = 3d0
  call store_factorized_event_measure(0, 0, measure)
  call store_factorized_event_measure(0, 2, measure)
  call scale_factorized_event_measures(5d0, 2)
  call fetch_factorized_event_measure(0, 0, fetched, available)
  if (.not. available .or. abs(fetched%jacobian - 2d0) > 1d-12) stop 9
  call fetch_factorized_event_measure(0, 2, fetched, available)
  if (.not. available .or. abs(fetched%jacobian - 10d0) > 1d-12) stop 10
  call reset_factorized_phase_space()
  slots = 0
  slots(0) = 3
  slots(1) = 0
  slots(2) = 3

  root = 0d0
  root(:, 1) = (/5d0, 0d0, 0d0, 5d0/)
  root(:, 2) = (/5d0, 0d0, 0d0, -5d0/)
  root(:, 3) = (/5d0, 0d0, 0d0, 3d0/)
  call store_factorized_local_momenta(3, 0, 3, root)
  pdgs = 0
  pdgs(1:3) = (/2, -2, 6/)
  finals = .false.
  finals(3) = .true.
  kinds = factorized_no_target
  kinds(3) = factorized_block_target
  targets = 0
  targets(3) = 1
  call store_factorized_local_layout( &
       3, 0, 3, pdgs, finals, kinds, targets)

  first = 0d0
  first(:, 1) = (/4d0, 0d0, 0d0, 0d0/)
  first(:, 2) = (/2.5d0, 0d0, 0d0, 1.5d0/)
  first(:, 3) = (/1.5d0, 0d0, 0d0, -1.5d0/)
  call store_factorized_local_momenta(0, 1, 3, first)
  pdgs = 0
  pdgs(1:3) = (/6, 24, 5/)
  finals = .true.
  finals(1) = .false.
  kinds = factorized_visible_target
  kinds(1) = factorized_no_target
  kinds(2) = factorized_block_target
  targets = (/0, 2, 3, 0, 0/)
  call store_factorized_local_layout( &
       0, 1, 3, pdgs, finals, kinds, targets)

  second = 0d0
  second(:, 1) = (/2d0, 0d0, 0d0, 0d0/)
  second(:, 2) = (/1d0, 1d0, 0d0, 0d0/)
  second(:, 3) = (/1d0, -1d0, 0d0, 0d0/)
  call store_factorized_local_momenta(3, 2, 3, second)
  pdgs = 0
  pdgs(1:3) = (/24, 11, -12/)
  finals = .true.
  finals(1) = .false.
  kinds = factorized_visible_target
  kinds(1) = factorized_no_target
  targets = (/0, 4, 5, 0, 0/)
  call store_factorized_local_layout( &
       3, 2, 3, pdgs, finals, kinds, targets)

  call realize_factorized_event_tuple(slots, pass)
  if (.not. pass) stop 11
  call fetch_factorized_block_momenta( &
       0, 1, 3, realized_first, available)
  if (.not. available) stop 12
  if (maxval(abs(realized_first(:, 1) - root(:, 3))) > 1d-12) stop 13
  call fetch_factorized_matrix_momenta( &
       0, 1, 2, realized_first(:, 1:2), available)
  if (.not. available) stop 24
  if (maxval(abs(realized_first(:, 2) - first(:, 2))) < 1d-12) stop 25
  call fetch_factorized_block_momenta( &
       3, 2, 3, realized_second, available)
  if (.not. available) stop 14
  if (maxval(abs(realized_second(:, 1) - &
                 realized_first(:, 2))) > 1d-12) stop 15
  previous_second_parent = realized_second(:, 1)

  root_revision = factorized_block_momentum_revision(3, 0)
  first_revision = factorized_block_momentum_revision(0, 1)
  second_revision = factorized_block_momentum_revision(3, 2)
  previous_slots = slots
  call realize_factorized_event_transition( &
       slots, previous_slots, .true., pass)
  if (.not. pass) stop 26
  if (factorized_block_momentum_revision(3, 0) /= root_revision .or. &
      factorized_block_momentum_revision(0, 1) /= first_revision .or. &
      factorized_block_momentum_revision(3, 2) /= second_revision) stop 27

  ! A leaf-only slot transition must not rewrite either ancestor cache.
  call store_factorized_local_momenta(0, 2, 3, second)
  call store_factorized_local_layout( &
       0, 2, 3, pdgs, finals, kinds, targets)
  slots(2) = 0
  call realize_factorized_event_transition( &
       slots, previous_slots, .true., pass)
  if (.not. pass) stop 28
  if (factorized_block_momentum_revision(3, 0) /= root_revision .or. &
      factorized_block_momentum_revision(0, 1) /= first_revision) stop 29
  if (factorized_block_momentum_revision(0, 2) <= 0_8) stop 30
  call fetch_factorized_block_momenta( &
       0, 2, 3, realized_second, available)
  if (.not. available) stop 31
  previous_second_parent = realized_second(:, 1)

  call factorized_production_scale_sums(slots, sum_et, sum_mt)
  if (abs(sum_et) > 1d-12 .or. abs(sum_mt - 4d0) > 1d-12) stop 21
  call factorized_visible_scale_sums(slots, sum_et, sum_mt)
  if (abs(sum_et - 2d0) > 1d-12 .or. abs(sum_mt - 2d0) > 1d-12) stop 22

  call materialize_factorized_event_tuple( &
       slots, 5, visible_count, visible, visible_pdgs, pass, &
       visible_origins)
  if (.not. pass .or. visible_count /= 5) stop 16
  if (any(visible_pdgs /= (/2, -2, 5, 11, -12/))) stop 17
  if (any(visible_origins /= (/0, 0, 1, 2, 2/))) stop 23

  root(:, 3) = (/sqrt(32d0), 0d0, 0d0, 4d0/)
  call store_factorized_local_momenta(3, 0, 3, root)
  call realize_factorized_event_tuple(slots, pass)
  if (.not. pass) stop 18
  call fetch_factorized_block_momenta( &
       slots(2), 2, 3, realized_second, available)
  if (.not. available) stop 19
  if (maxval(abs(realized_second(:, 1) - &
                 previous_second_parent)) < 1d-6) stop 20
end program test_nested_tuple_boosts
'''
        with tempfile.TemporaryDirectory() as directory:
            stub_path = os.path.join(directory, 'stubs.f90')
            program_path = os.path.join(directory, 'test.f90')
            executable = os.path.join(directory, 'test')
            with open(stub_path, 'w') as stream:
                stream.write(stubs)
            with open(program_path, 'w') as stream:
                stream.write(program)
            subprocess.check_call([
                self.compiler, '-std=f2008', '-ffree-line-length-none',
                '-J', directory, '-I', directory, stub_path] +
                self.sources + [program_path, '-o', executable])
            subprocess.check_call([executable])


if __name__ == '__main__':
    unittest.main()
