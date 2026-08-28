from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class MultiplicativeContributionBundleTest(unittest.TestCase):

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_independent_channels_and_event_targets(self):
        repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        source = os.path.join(
            repository, 'Template', 'fNLO', 'SubProcesses',
            'nlo_contribution_bundle.f90')
        stubs = r'''
module process_dimensions
  implicit none
  integer, parameter :: nexternal = 7
  integer, parameter :: nincoming = 2
  integer, parameter :: fks_configs = 6
  integer, parameter :: nsplitorders = 1
  integer, parameter :: amp_split_size = 1
  integer, parameter :: amp_split_orders(1, 1) = reshape((/4/), (/1, 1/))
contains
  subroutine validate_process_dimensions()
  end subroutine validate_process_dimensions
end module process_dimensions

module fks_metadata
  implicit none
contains
  integer function fks_j_d(configuration)
    integer, intent(in) :: configuration
    integer, parameter :: values(6) = (/1, 3, 4, 2, 5, 6/)
    fks_j_d = values(configuration)
  end function fks_j_d
end module fks_metadata

module fnlo_process_common
  implicit none
  integer :: nfksprocess = 1
end module fnlo_process_common
'''
        program = r'''
program test_multiplicative_bundle
  use fnlo_process_common, only: nfksprocess
  use nlo_contribution_bundle
  implicit none
  integer :: unit_number

  open(newunit=unit_number, file='nlo_contribution_info.dat', &
       status='replace', action='write')
  write(unit_number, '(a)') 'FORMAT 3'
  write(unit_number, '(a)') 'COUNT 3'
  write(unit_number, '(a)') 'VIRTUAL_GRIDS 0'
  write(unit_number, '(a)') 'CONTRIBUTION 1 PRODUCTION 1 2 1 0 0 0 0'
  write(unit_number, '(a)') 'CONTRIBUTION 2 NLO_DECAY 3 4 3 0 6 1 1'
  write(unit_number, '(a)') 'CONTRIBUTION 3 NLO_DECAY 5 6 5 0 24 1 2'
  write(unit_number, '(a)') 'END'
  close(unit_number)

  if (.not. has_nlo_contribution_bundle()) stop 1
  if (nlo_contribution_count() /= 3) stop 2
  if (contribution_fks_channel_count(1, 0) /= 2) stop 3
  if (contribution_fks_channel_count(1, 1) /= 1) stop 4
  if (contribution_fks_channel_count(1, 2) /= 1) stop 5
  if (contribution_fks_channel_configuration(1, 2, 1) /= 1) stop 6
  if (contribution_fks_channel_configuration(1, 1, 1) /= 2) stop 7
  if (contribution_fks_channel_count(2, 1) /= 1) stop 8
  if (contribution_fks_channel_count(2, 2) /= 1) stop 9
  if (contribution_fks_channel_count(3, 1) /= 2) stop 10
  if (contribution_fks_channel_count(3, 2) /= 0) stop 11
  if (multiplicative_event_capacity() /= 9) stop 12
  if (multiplicative_emission_target(1) /= 7) stop 13
  if (multiplicative_emission_target(2) /= 8) stop 14
  if (multiplicative_emission_target(3) /= 9) stop 15
  if (multiplicative_mc_integer_dimension(1,0) /= 1) stop 17
  if (multiplicative_mc_integer_dimension(1,2) /= 3) stop 18
  if (multiplicative_mc_integer_dimension(3,0) /= 7) stop 19
  nfksprocess = 4
  if (active_nlo_contribution() /= 2) stop 16
end program test_multiplicative_bundle
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
            subprocess.check_call([executable], cwd=directory)


if __name__ == '__main__':
    unittest.main()
