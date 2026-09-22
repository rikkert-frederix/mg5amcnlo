"""Small massive shells survive boosts in the loop reduction interface."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR
from madgraph.iolibs.file_writers import FortranWriter


class TestCTKinematicMatrix(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler = shutil.which('gfortran')
        if not compiler:
            raise unittest.SkipTest('gfortran is required')
        temporary = tempfile.TemporaryDirectory(prefix='ct_small_mass_')
        cls.addClassCleanup(temporary.cleanup)
        folder = Path(temporary.name)
        source = (Path(MG5DIR)/'madgraph/iolibs/template_files/loop_optimized/CT_interface.inc').read_text()
        start = source.index('SUBROUTINE %(proc_prefix)sBUILD_KINEMATIC_MATRIX')
        source = source[start:source.index('## if(samurai_available)', start)]
        source = source % dict(proc_prefix='TEST_', real_dp_format='REAL*8',
            complex_dp_format='COMPLEX*16', real_mp_format='REAL*16', complex_mp_format='COMPLEX*32')
        writer = FortranWriter(str(folder/'matrix.f'))
        writer.writelines(source.splitlines())
        writer.close()
        (folder/'MadLoopParams.inc').write_text(
            '      REAL*8 OSTHRES\n      COMMON /TEST_THRESHOLD/ OSTHRES\n')
        (folder/'check.f90').write_text('''program check
  implicit none
  real*8 :: p(3,0:3), osthres
  real*16 :: qp(3,0:3)
  complex*16 :: m(3), s(3,3)
  complex*32 :: qm(3), qs(3,3)
  common /test_threshold/ osthres
  character(8) :: mode
  integer :: i
  call get_command_argument(1,mode)
  osthres=1d-8
  ! The stopped mb=0.1 GeV decay, with loop routing (0,pt,pW).
  p(1,:)=0d0
  p(2,:)=[3366.09504132084339d0,-1415.31485370776841d0, &
          1400.95796287141047d0,2708.32793958673301d0]
  p(3,:)=[2428.78109674828511d0,-1056.18021336619813d0, &
          1050.28628430388153d0,1916.74149159287572d0]
  m=cmplx([172.5d0**2,0d0,0.1d0**2],0d0,kind=8)
  qp=p
  qm=cmplx([172.5_16**2,0.0_16,0.1_16**2],0.0_16,kind=16)
  ! MadLoop repairs the on-shell momenta before its QP evaluation.
  qp(2,0)=qp(3,0)+sqrt(sum((qp(2,1:3)-qp(3,1:3))**2)+real(qm(3),kind=16))
  if (mode=='dp') then
    call test_build_kinematic_matrix(3,p,m,s)
    if (abs(s(2,3))>1d-14.or.abs(s(3,2))>1d-14) stop 11
    if (maxval(abs(s-transpose(s)))>8d0*epsilon(1d0)*maxval(abs(s))) stop 12
    do i=1,3
      if (s(i,i)/=-2d0*m(i)) stop 13
    end do
  else
    call test_mp_build_kinematic_matrix(3,qp,qm,qs)
    if (abs(qs(2,3))>1e-28_16.or.abs(qs(3,2))>1e-28_16) stop 21
    if (maxval(abs(qs-transpose(qs)))>1e-25_16) stop 22
    do i=1,3
      if (qs(i,i)/=-2.0_16*qm(i)) stop 23
    end do
  end if
  ! An exactly massless difference still receives the massless shell.
  p=0d0
  p(2,:)=[1000d0,0d0,0d0,1000d0]
  qp=p
  call test_build_kinematic_matrix(3,p,m,s)
  call test_mp_build_kinematic_matrix(3,qp,qm,qs)
  if (abs(s(2,3)+m(2)+m(3))>1d-14) stop 31
  if (abs(qs(2,3)+qm(2)+qm(3))>1e-28_16) stop 32
  ! A resolved off-shell invariant must remain off shell.
  p(2,:)=[2d0,0d0,0d0,0d0]
  qp=p
  call test_build_kinematic_matrix(3,p,m,s)
  call test_mp_build_kinematic_matrix(3,qp,qm,qs)
  if (abs(s(2,3)+m(2)+m(3)-4d0)>1d-14) stop 41
  if (abs(qs(2,3)+qm(2)+qm(3)-4.0_16)>1e-28_16) stop 42
  write(*,*) 'MASS_SHELL_PASS'
end program
''')
        cls.executable = folder/'check'
        result = subprocess.run([compiler, '-O0', '-fcheck=all', '-ffixed-line-length-132',
            'matrix.f', 'check.f90', '-o', str(cls.executable)], cwd=folder,
            text=True, capture_output=True)
        if result.returncode:
            raise AssertionError(result.stdout+result.stderr)

    def test_boosted_small_mass_double_precision(self):
        self.check_mode('dp')

    def test_boosted_small_mass_quadruple_precision(self):
        self.check_mode('qp')

    def check_mode(self, mode):
        result = subprocess.run([str(self.executable), mode], text=True, capture_output=True)
        self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
        self.assertIn('MASS_SHELL_PASS', result.stdout)
