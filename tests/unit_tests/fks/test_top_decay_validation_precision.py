"""The independent virtual comparison must request matching precision."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR


class TestTopDecayValidationPrecision(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler=shutil.which('gfortran')
        if not compiler:
            raise unittest.SkipTest('gfortran is required')
        cls.temporary=tempfile.TemporaryDirectory(prefix='top_virtual_precision_')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.directory=Path(cls.temporary.name)
        stub=cls.directory/'bundle.f90'
        stub.write_text('module nlo_contribution_bundle\ncontains\n'
            'integer function nlo_contribution_count()\nnlo_contribution_count=3\n'
            'end function\nend module\n')
        driver=cls.directory/'check.f90'
        driver.write_text('''program check
  use top_decay_virtual_dispatch
  implicit none
  double precision :: p(0:3,4)
  complex(kind=8) :: analytic(3,2,2), reference(3,2,2)
  character(len=16) :: mode
  integer :: point
  if (abs(tdv_validation_precision(1d-3)-1d-10)>1d-25) stop 10
  if (abs(tdv_validation_precision(-1d0)-1d-10)>1d-25) stop 11
  if (abs(tdv_validation_precision(0d0)-1d-10)>1d-25) stop 12
  if (tdv_validation_precision(1d-12)/=1d-12) stop 13
  call get_command_argument(1,mode)
  analytic=cmplx(1d0,0d0,kind=8)
  reference=analytic
  if (mode=='bad_density') reference(1,1,1)=cmplx(1d0+2d-8,0d0,kind=8)
  do point=1,tdv_required_validation_points
    if (.not.tdv_madloop_required(2,.true.)) stop 14
    p=dble(point)
    call tdv_validate_against_madloop(2,p,analytic,reference,1d-12,217)
  end do
  if (tdv_madloop_required(2,.true.)) stop 15
  if (.not.tdv_madloop_required(2,.false.)) stop 16
  write(*,*) 'VALIDATION_PRECISION_PASS'
end program
''')
        source=Path(MG5DIR)/'Template/fNLO/SubProcesses'
        cls.executable=cls.directory/'check'
        result=subprocess.run([compiler,'-O0','-g','-fcheck=all','-ffree-line-length-none',
            str(stub),str(source/'top_decay_virtual_cdr.f90'),
            str(source/'top_decay_virtual_dispatch.f90'),str(Path(MG5DIR)/'HELAS/vxxxxx.F'),
            str(driver),'-o',str(cls.executable)],
            cwd=cls.directory,text=True,capture_output=True)
        if result.returncode:
            raise AssertionError(result.stdout+result.stderr)

    def test_requested_precision_and_three_point_switch(self):
        result=subprocess.run([str(self.executable)],text=True,capture_output=True)
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        self.assertIn('VALIDATION_PRECISION_PASS',result.stdout)

    def test_density_mismatch_still_fails_at_original_tolerance(self):
        result=subprocess.run([str(self.executable),'bad_density'],text=True,capture_output=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('analytic top-decay virtual validation failed',result.stdout)
        self.assertNotIn('VALIDATION_PRECISION_PASS',result.stdout)
