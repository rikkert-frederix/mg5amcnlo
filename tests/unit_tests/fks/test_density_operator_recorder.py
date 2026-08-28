from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class DensityOperatorRecorderTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.compiler = shutil.which('gfortran')
        cls.repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        cls.sources = [os.path.join(
            cls.repository, 'Template', 'fNLO', 'SubProcesses', name)
            for name in (
                'multiplicative_density_terms.f90',
                'density_operator_recorder.f90')]

    def compile_program(self, directory, name, program):
        stubs = os.path.join(directory, 'stubs.f90')
        with open(stubs, 'w') as stream:
            stream.write(r'''
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
            self.compiler, '-std=f2008', '-ffree-line-length-none',
            '-J', directory, '-I', directory, stubs] + self.sources +
            [program_path, '-o', executable])
        return executable

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_one_recording_produces_one_same_momentum_term(self):
        program = r'''
program test_recorder
  use multiplicative_density_terms
  use density_operator_recorder
  implicit none
  type(block_distribution_term) :: term
  complex(kind=8) :: coefficients(3), factor(3)

  coefficients = (0d0, 0d0)
  coefficients(1) = (2d0, 0d0)
  factor = (0d0, 0d0)
  factor(1) = (3d0, 0d0)
  factor(2) = (5d0, 0d0)

  call begin_density_operator_recording(2, 3, 1)
  call record_density_operator_primitive(2, 7, 1, 0, coefficients, .true.)
  call scale_recorded_density_operator(factor)
  if (recorded_density_operator_count() /= 1) stop 11
  call finish_density_operator_recording(term, -1)

  if (density_operator_is_recording()) stop 12
  if (.not. term%finalized) stop 13
  if (term%block /= 2 .or. term%event_slot /= 3) stop 14
  if (term%sign /= -1 .or. term%nlo_order /= 1) stop 15
  if (term%primitive_count /= 1) stop 16
  if (abs(term%primitives(1)%scale_coefficients(1) - &
          (6d0, 0d0)) > 1d-12) stop 17
  if (abs(term%primitives(1)%scale_coefficients(2) - &
          (10d0, 0d0)) > 1d-12) stop 18
  if (term%primitives(1)%scale_coefficients(3) /= (0d0, 0d0)) stop 19
end program test_recorder
'''
        with tempfile.TemporaryDirectory() as directory:
            executable = self.compile_program(directory, 'recorder', program)
            subprocess.check_call([executable])

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_quadratic_scale_logs_are_rejected(self):
        program = r'''
program test_quadratic_logs
  use density_operator_recorder
  implicit none
  complex(kind=8) :: coefficients(3), factor(3)

  coefficients = (0d0, 0d0)
  coefficients(1) = (1d0, 0d0)
  coefficients(2) = (2d0, 0d0)
  factor = (0d0, 0d0)
  factor(1) = (3d0, 0d0)
  factor(2) = (4d0, 0d0)
  call begin_density_operator_recording(0, 0, 1)
  call record_density_operator_primitive(3, 1, 1, 0, coefficients, .true.)
  call scale_recorded_density_operator(factor)
end program test_quadratic_logs
'''
        with tempfile.TemporaryDirectory() as directory:
            executable = self.compile_program(
                directory, 'quadratic_logs', program)
            result = subprocess.run(
                [executable], stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, universal_newlines=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('quadratic scale logs', result.stdout)


if __name__ == '__main__':
    unittest.main()
