from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class MultiplicativeLambdaValidationTest(unittest.TestCase):

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_global_linear_and_mixed_formal_lambda_coefficients(self):
        repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        source = os.path.join(
            repository, 'Template', 'fNLO', 'SubProcesses',
            'multiplicative_lambda_validation.f90')
        program = r'''
program test_lambda_validation
  use multiplicative_lambda_validation
  implicit none
  type(multiplicative_lambda_accumulator) :: accumulator
  integer :: orders(3), width_blocks(2)
  double precision :: lo_widths(2), nlo_widths(2)
  double precision :: epsilon, finite_difference, value

  call initialize_multiplicative_lambda_accumulator(accumulator, 3)
  orders = (/0, 0, 0/)
  call accumulate_multiplicative_lambda_atom(accumulator, orders, 100d0)
  orders = (/1, 0, 0/)
  call accumulate_multiplicative_lambda_atom(accumulator, orders, 10d0)
  orders = (/0, 1, 0/)
  call accumulate_multiplicative_lambda_atom(accumulator, orders, 20d0)
  orders = (/0, 0, 1/)
  call accumulate_multiplicative_lambda_atom(accumulator, orders, 30d0)
  orders = (/1, 1, 0/)
  call accumulate_multiplicative_lambda_atom(accumulator, orders, 2d0)
  orders = (/1, 0, 1/)
  call accumulate_multiplicative_lambda_atom(accumulator, orders, 3d0)
  orders = (/0, 1, 1/)
  call accumulate_multiplicative_lambda_atom(accumulator, orders, 6d0)
  orders = (/1, 1, 1/)
  call accumulate_multiplicative_lambda_atom(accumulator, orders, 1d0)

  lo_widths = (/2d0, 4d0/)
  nlo_widths = (/2.2d0, 3.6d0/)
  width_blocks = (/2, 3/)

  if (abs(accumulator%exact_weight - 172d0) > 1d-12) stop 1
  if (maxval(abs(accumulator%order_coefficients - &
      (/100d0, 60d0, 11d0, 1d0/))) > 1d-12) stop 2
  if (abs(formal_lambda_global_weight( &
      accumulator, 0d0, lo_widths, nlo_widths) - 99d0) > 1d-12) stop 3
  if (abs(formal_lambda_global_weight( &
      accumulator, 1d0, lo_widths, nlo_widths) - 172d0) > 1d-12) stop 4
  if (abs(formal_lambda_lo_weight( &
      accumulator, lo_widths, nlo_widths) - 99d0) > 1d-12) stop 5
  if (abs(formal_lambda_linear_correction( &
      accumulator, lo_widths, nlo_widths) - 59.4d0) > 1d-11) stop 6
  if (abs(formal_lambda_additive_weight( &
      accumulator, lo_widths, nlo_widths) - 158.4d0) > 1d-11) stop 7

  value = formal_lambda_block_linear_correction( &
      accumulator, 1, width_blocks, lo_widths, nlo_widths)
  if (abs(value - 9.9d0) > 1d-11) stop 8
  value = formal_lambda_block_linear_correction( &
      accumulator, 2, width_blocks, lo_widths, nlo_widths)
  if (abs(value - 9.9d0) > 1d-11) stop 9
  value = formal_lambda_block_linear_correction( &
      accumulator, 3, width_blocks, lo_widths, nlo_widths)
  if (abs(value - 39.6d0) > 1d-11) stop 10

  value = formal_lambda_block_mixed_coefficient( &
      accumulator, 1, 2, width_blocks, lo_widths, nlo_widths)
  if (abs(value - 0.99d0) > 1d-11) stop 11
  value = formal_lambda_block_mixed_coefficient( &
      accumulator, 1, 3, width_blocks, lo_widths, nlo_widths)
  if (abs(value - 3.96d0) > 1d-11) stop 12
  value = formal_lambda_block_mixed_coefficient( &
      accumulator, 2, 3, width_blocks, lo_widths, nlo_widths)
  if (abs(value - 3.96d0) > 1d-11) stop 13

  epsilon = 1d-5
  finite_difference = (formal_lambda_global_weight( &
      accumulator, epsilon, lo_widths, nlo_widths) - &
      formal_lambda_global_weight( &
      accumulator, -epsilon, lo_widths, nlo_widths))/(2d0*epsilon)
  if (abs(finite_difference - formal_lambda_linear_correction( &
      accumulator, lo_widths, nlo_widths)) > 1d-7) stop 14

  call require_formal_lambda_closure( &
      accumulator, 172d0, lo_widths, nlo_widths)
  call require_formal_lambda_linear_closure( &
      accumulator, width_blocks, lo_widths, nlo_widths)
end program test_lambda_validation
'''
        with tempfile.TemporaryDirectory() as directory:
            program_path = os.path.join(directory, 'test.f90')
            executable = os.path.join(directory, 'test')
            with open(program_path, 'w') as stream:
                stream.write(program)
            subprocess.check_call([
                shutil.which('gfortran'), '-std=f2008',
                '-ffree-line-length-none', '-J', directory, '-I', directory,
                source, program_path, '-o', executable])
            subprocess.check_call([executable])


if __name__ == '__main__':
    unittest.main()
