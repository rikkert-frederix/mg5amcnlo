from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class FactorizedMintPolicyTest(unittest.TestCase):

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_process_policy_selects_finite_block_local_proxies(self):
        repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        source = os.path.join(
            repository, 'Template', 'fNLO', 'SubProcesses',
            'factorized_mint_policy.f90')
        stubs = r'''
module policy_state
  implicit none
  logical :: bundle = .false.
  logical :: multiplicative = .false.
end module policy_state

module nlo_contribution_bundle
  use policy_state
  implicit none
contains
  logical function has_nlo_contribution_bundle()
    has_nlo_contribution_bundle = bundle
  end function has_nlo_contribution_bundle
  integer function factorized_radiation_block(total_dimension, dimension)
    integer, intent(in) :: total_dimension, dimension
    factorized_radiation_block = 0
    if (dimension > total_dimension - 6) &
         factorized_radiation_block = (dimension - (total_dimension - 6) - 1)/3 + 1
  end function factorized_radiation_block
  integer function bundle_nlo_component(block)
    integer, intent(in) :: block
    bundle_nlo_component = block + 1
  end function bundle_nlo_component
end module nlo_contribution_bundle

module decay_chain_parameters
  use policy_state
  implicit none
contains
  logical function uses_multiplicative_nlo_combination()
    uses_multiplicative_nlo_combination = multiplicative
  end function uses_multiplicative_nlo_combination
end module decay_chain_parameters
'''
        program = r'''
program test_policy
  use policy_state
  use factorized_mint_policy
  implicit none
  double precision :: integrals(8), weight

  integrals = (/2d0, -3d0, 5d0, 7d0, -11d0, 13d0, 17d0, 19d0/)
  call factorized_mint_grid_weight(10, 1, integrals, 2, 4, weight)
  if (weight /= 2d0) stop 1

  multiplicative = .true.
  call factorized_mint_grid_weight(10, 2, integrals, 2, 4, weight)
  if (weight /= 3d0) stop 2
  if (.not. factorized_mint_uses_uniform_channels()) stop 3
  if (factorized_mint_shows_multiplicative_validation()) stop 4

  bundle = .true.
  call factorized_mint_grid_weight(10, 3, integrals, 2, 4, weight)
  if (weight /= 3d0) stop 5
  call factorized_mint_grid_weight(10, 5, integrals, 2, 4, weight)
  if (weight /= -11d0) stop 6
  call factorized_mint_grid_weight(10, 8, integrals, 2, 4, weight)
  if (weight /= 13d0) stop 7
  if (.not. factorized_mint_shows_multiplicative_validation()) stop 8
end program test_policy
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
                shutil.which('gfortran'), '-std=f2008',
                '-ffree-line-length-none', '-J', directory, '-I', directory,
                stub_path, source, program_path, '-o', executable])
            subprocess.check_call([executable])


if __name__ == '__main__':
    unittest.main()
