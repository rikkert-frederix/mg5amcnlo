from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class MultiplicativeDensityTermsTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.compiler = shutil.which('gfortran')
        cls.repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        cls.source = os.path.join(
            cls.repository, 'Template', 'fNLO', 'SubProcesses',
            'multiplicative_density_terms.f90')

    def compile_program(self, directory, name, program):
        stub = os.path.join(directory, 'stubs.f90')
        with open(stub, 'w') as stream:
            stream.write('''
module process_dimensions
  implicit none
  integer, parameter :: nexternal = 4
contains
  subroutine validate_process_dimensions()
  end subroutine validate_process_dimensions
end module process_dimensions

module fnlo_process_common
  implicit none
  integer, parameter :: soft_counterevent = 0
  integer, parameter :: real_event = 3
end module fnlo_process_common

module spin_density_matrix_results
  implicit none
  integer, parameter :: spin_density_no_insertion = 0
  integer, parameter :: spin_density_color_insertion = 4
  integer, parameter :: spin_density_fast_virtual_insertion = 5
end module spin_density_matrix_results
''')
        program_path = os.path.join(directory, name + '.f90')
        with open(program_path, 'w') as stream:
            stream.write(program)
        executable = os.path.join(directory, name)
        subprocess.check_call([
            self.compiler, '-std=f2008', '-fno-automatic', '-J', directory,
            '-I', directory, stub, self.source, program_path, '-o',
            executable])
        return executable

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_scheduled_terms_retain_independent_event_slots(self):
        program = r'''
program test_density_tuples
  use multiplicative_density_terms
  implicit none
  type(block_nlo_distribution) :: distributions(2)
  type(multiplicative_density_tuple) :: tuple
  type(density_tuple_schedule) :: schedule
  complex(kind=8) :: coefficients(3), evaluated
  integer :: distribution, term, tuple_index
  integer :: expected_signs(4)
  integer :: expected_slots(2, 4)

  coefficients = (0d0, 0d0)
  coefficients(1) = (2d0, 0d0)
  expected_signs = (/1, -1, -1, 1/)
  expected_slots(:, 1) = (/0, 0/)
  expected_slots(:, 2) = (/0, 3/)
  expected_slots(:, 3) = (/3, 0/)
  expected_slots(:, 4) = (/3, 3/)

  do distribution = 1, 2
    call initialize_block_distribution( &
         distributions(distribution), distribution - 1, 2)
    call initialize_block_distribution_term( &
         distributions(distribution), 1, 3, 1, 1, 1)
    call initialize_block_distribution_term( &
         distributions(distribution), 2, 0, -1, 1, 1)
    do term = 1, 2
      call set_density_primitive( &
           distributions(distribution)%terms(term), 1, 2, 1, 1, 0, 1, &
           coefficients, .true.)
      call finalize_block_distribution_term( &
           distributions(distribution)%terms(term))
    end do
    call finalize_block_distribution(distributions(distribution))
  end do

  if (density_cartesian_tuple_count(distributions) /= 4) stop 11
  call initialize_density_tuple_schedule(distributions, schedule)
  do tuple_index = 1, 4
    call decode_scheduled_density_tuple( &
         distributions, schedule, tuple_index, tuple)
    if (tuple%sign /= expected_signs(tuple_index)) stop 12
    if (tuple%nlo_order /= 2) stop 13
    if (tuple%event_slots(0) /= expected_slots(1, tuple_index)) stop 14
    if (tuple%event_slots(1) /= expected_slots(2, tuple_index)) stop 15
  end do
  evaluated = evaluate_density_primitive_coefficient( &
       distributions(1)%terms(1)%primitives(1), 3d0, 5d0)
  if (abs(evaluated - (2d0, 0d0)) > 1d-12) stop 16
end program test_density_tuples
'''
        with tempfile.TemporaryDirectory() as directory:
            executable = self.compile_program(directory, 'tuples', program)
            subprocess.check_call([executable])

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_raw_laurent_primitive_is_rejected(self):
        program = r'''
program test_raw_laurent_rejection
  use multiplicative_density_terms
  implicit none
  type(block_nlo_distribution) :: distribution
  complex(kind=8) :: coefficients(3)

  coefficients = (0d0, 0d0)
  coefficients(1) = (1d0, 0d0)
  call initialize_block_distribution(distribution, 0, 1)
  call initialize_block_distribution_term(distribution, 1, 0, 1, 1, 1)
  call set_density_primitive(distribution%terms(1), 1, 3, 1, 1, 0, 1, &
       coefficients, .false.)
  call finalize_block_distribution_term(distribution%terms(1))
end program test_raw_laurent_rejection
'''
        with tempfile.TemporaryDirectory() as directory:
            executable = self.compile_program(directory, 'raw', program)
            result = subprocess.run(
                [executable], stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, universal_newlines=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('raw Laurent-pole density', result.stdout)

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_exact_kinematic_families_fold_signs_and_orders(self):
        program = r'''
program test_family_coalescing
  use multiplicative_density_terms
  implicit none
  type(block_nlo_distribution) :: distribution
  type(multiplicative_density_tuple) :: tuple
  type(density_tuple_schedule) :: schedule
  complex(kind=8) :: coefficients(3)

  call initialize_block_distribution(distribution, 0, 3)
  coefficients = (0d0, 0d0)
  coefficients(1) = (2d0, 0d0)
  call fill_term(1, 0, 7, 1, 0, 1, coefficients)
  coefficients(1) = (3d0, 0d0)
  call fill_term(2, 0, 7, -1, 1, 2, coefficients)
  coefficients(1) = (5d0, 0d0)
  call fill_term(3, 3, 7, 1, 1, 3, coefficients)
  call finalize_block_distribution(distribution)

  call coalesce_block_kinematic_families(distribution)
  if (.not. distribution%kinematic_families_coalesced) stop 11
  if (distribution%term_count /= 2) stop 12
  if (distribution%terms(1)%event_slot /= 0 .or. &
      distribution%terms(1)%luminosity_configuration /= 7) stop 13
  if (distribution%terms(1)%sign /= 1 .or. &
      distribution%terms(1)%nlo_order /= -1 .or. &
      distribution%terms(1)%primitive_count /= 2) stop 14
  if (abs(distribution%terms(1)%primitives(1)% &
      scale_coefficients(1) - (2d0, 0d0)) > 1d-12) stop 15
  if (abs(distribution%terms(1)%primitives(2)% &
      scale_coefficients(1) + (3d0, 0d0)) > 1d-12) stop 16
  if (distribution%terms(1)%primitives(1)%radiation_group /= 1 .or. &
      distribution%terms(1)%primitives(2)%radiation_group /= 2) stop 17
  if (distribution%terms(2)%event_slot /= 3 .or. &
      distribution%terms(2)%nlo_order /= 1) stop 18
  if (density_cartesian_tuple_count((/distribution/)) /= 2) stop 19
  call initialize_density_tuple_schedule((/distribution/), schedule)
  call decode_scheduled_density_tuple( &
       (/distribution/), schedule, 1, tuple)
  if (tuple%nlo_order /= -1) stop 20

  ! Re-entering the coalescer must be safe with the production compiler's
  ! -fno-automatic flag, and identical descriptors in one radiation group
  ! must load just one matrix with the sum of their coefficients.
  call initialize_block_distribution(distribution, 0, 2)
  coefficients = (0d0, 0d0)
  coefficients(1) = (7d0, 0d0)
  call fill_term(1, 0, 4, 1, 1, 2, coefficients)
  coefficients(1) = (11d0, 0d0)
  call fill_term(2, 0, 4, -1, 1, 2, coefficients)
  call finalize_block_distribution(distribution)
  call coalesce_block_kinematic_families(distribution)
  if (distribution%term_count /= 1) stop 21
  if (distribution%terms(1)%primitive_count /= 1) stop 22
  if (abs(distribution%terms(1)%primitives(1)% &
      scale_coefficients(1) + (4d0, 0d0)) > 1d-12) stop 23

contains

  subroutine fill_term(index, slot, luminosity, sign, order, group, values)
    integer, intent(in) :: index, slot, luminosity, sign, order, group
    complex(kind=8), intent(in) :: values(3)

    call initialize_block_distribution_term( &
         distribution, index, slot, sign, order, 1)
    distribution%terms(index)%luminosity_configuration = luminosity
    call set_density_primitive( &
         distribution%terms(index), 1, 0, 0, 0, 0, order, values, .true.)
    distribution%terms(index)%primitives(1)%radiation_group = group
    call finalize_block_distribution_term(distribution%terms(index))
  end subroutine fill_term
end program test_family_coalescing
'''
        with tempfile.TemporaryDirectory() as directory:
            executable = self.compile_program(
                directory, 'families', program)
            subprocess.check_call([executable])


if __name__ == '__main__':
    unittest.main()
