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
end module spin_density_matrix_results
''')
        program_path = os.path.join(directory, name + '.f90')
        with open(program_path, 'w') as stream:
            stream.write(program)
        executable = os.path.join(directory, name)
        subprocess.check_call([
            self.compiler, '-std=f2008', '-J', directory, '-I', directory,
            stub, self.source, program_path, '-o', executable])
        return executable

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_cartesian_terms_retain_independent_event_slots(self):
        program = r'''
program test_density_tuples
  use multiplicative_density_terms
  implicit none
  type(block_nlo_distribution) :: distributions(2)
  type(multiplicative_density_tuple) :: tuple
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
  do tuple_index = 1, 4
    call decode_density_cartesian_tuple( &
         distributions, tuple_index, tuple)
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


if __name__ == '__main__':
    unittest.main()
