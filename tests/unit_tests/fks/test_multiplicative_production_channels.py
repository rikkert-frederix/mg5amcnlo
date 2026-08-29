from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class MultiplicativeProductionChannelsTest(unittest.TestCase):

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_amplitude_partition_is_normalized_inside_each_fks_category(self):
        repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        source = os.path.join(
            repository, 'Template', 'fNLO', 'SubProcesses',
            'multiplicative_production_channels.f90')
        program = r'''
program test_production_channels
  use multiplicative_production_channels
  implicit none
  integer :: categories(4), configurations(4), configuration_map(0:3)
  integer :: channel
  double precision :: diagram_weights(2), symmetry_factors(4), total

  ! Configuration 3 is symmetry-equivalent to retained configuration 1.
  ! Both FKS categories integrate the same two retained production maps.
  categories = (/1, 2, 1, 2/)
  configurations = (/1, 1, 2, 2/)
  configuration_map = (/3, 1, 2, 1/)
  diagram_weights = (/2d0, 5d0/)
  symmetry_factors = (/2d0, 2d0, 1d0, 1d0/)

  do channel = 1, 2
    total = 0d0
    if (channel == 1) then
      total = multiplicative_production_channel_partition( &
           categories, configurations, 4, 1, diagram_weights, &
           configuration_map, symmetry_factors(1)) + &
           multiplicative_production_channel_partition( &
           categories, configurations, 4, 3, diagram_weights, &
           configuration_map, symmetry_factors(3))
    else
      total = multiplicative_production_channel_partition( &
           categories, configurations, 4, 2, diagram_weights, &
           configuration_map, symmetry_factors(2)) + &
           multiplicative_production_channel_partition( &
           categories, configurations, 4, 4, diagram_weights, &
           configuration_map, symmetry_factors(4))
    end if
    if (abs(total - 1d0) > 1d-12) stop 1
  end do

  diagram_weights = 0d0
  if (multiplicative_production_channel_partition( &
       categories, configurations, 4, 1, diagram_weights, &
       configuration_map, symmetry_factors(1)) /= 0d0) stop 2
end program test_production_channels
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
