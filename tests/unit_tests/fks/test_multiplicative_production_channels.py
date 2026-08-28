from __future__ import absolute_import

import os
import shutil
import subprocess
import tempfile
import unittest


class MultiplicativeProductionChannelsTest(unittest.TestCase):

    @unittest.skipUnless(shutil.which('gfortran'), 'gfortran is unavailable')
    def test_partition_is_normalized_inside_each_fks_category(self):
        repository = os.path.abspath(os.path.join(
            os.path.dirname(__file__), '..', '..', '..'))
        source = os.path.join(
            repository, 'Template', 'fNLO', 'SubProcesses',
            'multiplicative_production_channels.f90')
        program = r'''
program test_production_channels
  use multiplicative_production_channels
  implicit none
  integer :: categories(7), channel
  double precision :: total

  categories = (/1, 2, 1, 0, 2, 1, 0/)
  do channel = 1, 7
    select case (categories(channel))
    case (0)
      if (abs(multiplicative_production_channel_partition( &
           categories, 7, channel) - 0.5d0) > 1d-12) stop 1
    case (1)
      if (abs(multiplicative_production_channel_partition( &
           categories, 7, channel) - 1d0/3d0) > 1d-12) stop 2
    case (2)
      if (abs(multiplicative_production_channel_partition( &
           categories, 7, channel) - 0.5d0) > 1d-12) stop 3
    end select
  end do

  do channel = 0, 2
    total = 0d0
    if (channel == 0) then
      total = multiplicative_production_channel_partition( &
           categories, 7, 4) + &
           multiplicative_production_channel_partition(categories, 7, 7)
    else if (channel == 1) then
      total = multiplicative_production_channel_partition( &
           categories, 7, 1) + &
           multiplicative_production_channel_partition(categories, 7, 3) + &
           multiplicative_production_channel_partition(categories, 7, 6)
    else
      total = multiplicative_production_channel_partition( &
           categories, 7, 2) + &
           multiplicative_production_channel_partition(categories, 7, 5)
    end if
    if (abs(total - 1d0) > 1d-12) stop 4
  end do
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
