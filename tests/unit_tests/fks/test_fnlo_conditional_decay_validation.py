"""Check the complete production-real plus decay-real validation buffer."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR


class TestConditionalDecayValidation(unittest.TestCase):
    def test_missing_production_parton_is_restored_for_validation(self):
        compiler = shutil.which('gfortran')
        if compiler is None:
            self.skipTest('gfortran is required')
        root = Path(MG5DIR)

        def procedure(filename, name, kind='subroutine'):
            source = (root / 'Template/fNLO/SubProcesses' / filename).read_text()
            marker = f'  {kind} {name}('
            start = source.index(marker)
            marker = f'  end {kind.split()[-1]} {name}'
            end = source.index(marker, start) + len(marker)
            return source[start:end]

        helper = procedure('genps_fks.f90', 'check_nlo_decay_visible_event')
        selector = procedure('nlo_decay_kinematics.f90',
                             'real_production_decay_embedding', 'logical function')
        fetch = procedure('nlo_decay_kinematics.f90',
                          'fetch_real_production_decay_parents')
        # The actual check, cache selection and cache accessors run together.
        # Only unrelated driver state and fatal-error reporting are stubbed.
        harness = '''
module process_dimensions
  implicit none
  integer, parameter :: nexternal=5,nincoming=2,max_particles=5,max_branch=5
end module
module conditional_decay_fixture
  use process_dimensions
  use kin_functions_module, only: dot=>dot_impl
  use phase_space_kinematics, only: phspncheck_nocms
  use hard_real_kinematics
  implicit none
  logical :: pi_enabled=.true.,conditional_real=.false.
contains
''' + helper + '\n' + selector + '\n' + fetch + '''
  logical function multiplicative_nlo_enabled()
    multiplicative_nlo_enabled=pi_enabled
  end function
  logical function production_conditional_real_branch()
    production_conditional_real_branch=conditional_real
  end function
  subroutine fail_genps_fks(message)
    character(len=*),intent(in)::message
    print *,message
    stop 90
  end subroutine
end module
program check_conditional_decay
  use conditional_decay_fixture
  implicit none
  double precision :: full(0:3,6),reduced(0:3,5),masses(-5:5),extended(-5:6)
  double precision :: before(0:3,5),energy,pz
  logical :: pass
  integer :: i
  full(:,1)=(/50d0,0d0,0d0,50d0/)
  full(:,2)=(/50d0,0d0,0d0,-50d0/)
  full(:,3)=(/25d0,25d0,0d0,0d0/)
  full(:,4)=(/25d0,-25d0,0d0,0d0/)
  full(:,5)=(/25d0,0d0,25d0,0d0/)
  full(:,6)=(/25d0,0d0,-25d0,0d0/)
  ! Exercise the actual non-CMS production-real embedding.
  do i=1,6
    energy=full(0,i); pz=full(3,i)
    full(0,i)=1.25d0*(energy+0.6d0*pz)
    full(3,i)=1.25d0*(pz+0.6d0*energy)
  end do
  reduced=full(:,1:5); before=reduced; masses=0d0; extended=0d0
  call phspncheck_nocms(5,100d0,masses,reduced,pass)
  if (pass) stop 1
  ! The checker must accept a buffer larger than the subprocess nexternal.
  call phspncheck_nocms(6,100d0,extended,full,pass)
  if (.not.pass) stop 2
  call cache_hard_real_point(full,reduced,(/0.5d0,0.5d0/),(/0.4d0,0.4d0/), &
       0.2d0,0d0,0d0,1d0,1d0,0d0,6,21)
  ! 80 is the inverse-Born energy, not the invariant of these incoming beams.
  call check_nlo_decay_visible_event(reduced,masses,80d0,pass)
  if (.not.pass .or. any(reduced/=before)) stop 3
  call clear_hard_real_point()
  conditional_real=.true.
  call cache_production_real_parent_point(full,6,21)
  call check_nlo_decay_visible_event(reduced,masses,80d0,pass)
  if (.not.pass .or. any(reduced/=before)) stop 4
  ! The legacy/Born branch still rejects the incomplete physical event.
  conditional_real=.false.
  call check_nlo_decay_visible_event(reduced,masses,100d0,pass)
  if (pass) stop 5
  ! Restore conservation in the legacy five-leg event and check it passes.
  reduced(:,5)=full(:,5)+full(:,6)
  masses(5)=50d0
  call check_nlo_decay_visible_event(reduced,masses,100d0,pass)
  if (.not.pass) stop 6
  ! A corrupted omitted production parton must still fail validation.
  conditional_real=.true.
  full(1,6)=3d0
  call cache_production_real_parent_point(full,6,21)
  call check_nlo_decay_visible_event(before,0d0*masses,80d0,pass)
  if (pass) stop 7
  print *,'PASS: conditional decay complete-event conservation'
end program
'''
        with tempfile.TemporaryDirectory(prefix='fnlo_conditional_check_') as directory:
            directory = Path(directory)
            dimensions, fixture = harness.split('module conditional_decay_fixture', 1)
            (directory / 'dimensions.f90').write_text(dimensions)
            (directory / 'check.f90').write_text('module conditional_decay_fixture' + fixture)
            executable = directory / 'check'
            result = subprocess.run([
                compiler, '-O0', '-g', '-fcheck=all',
                '-ffpe-trap=invalid,zero,overflow', '-o', str(executable),
                str(directory / 'dimensions.f90'),
                str(root / 'Template/fNLO/Source/kin_functions.f90'),
                str(root / 'Template/fNLO/SubProcesses/phase_space_kinematics.f90'),
                str(root / 'Template/fNLO/SubProcesses/hard_real_kinematics.f90'),
                str(directory / 'check.f90'),
            ], cwd=directory, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            result = subprocess.run([str(executable)], cwd=directory,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn('PASS: conditional decay complete-event conservation', result.stdout)
