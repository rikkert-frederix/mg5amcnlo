"""The independent virtual comparison must request matching precision."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR
from madgraph.iolibs.export_spin_density import SpinDensityExporter


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
  if (tdv_reference_needs_rescue(analytic,reference)) stop 17
  reference(1,1,2)=cmplx(1d0,2d-8,kind=8)
  if (.not.tdv_reference_needs_rescue(analytic,reference)) stop 18
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
        layout=dict(mode='two_body',open_size=2,parent_pdg=6,parent_leg=1,
                    bottom_leg=2,vector_leg=3,bottom_mass_fortran='0D0')
        variant=dict(fortran_name='MOCK_MADLOOP',contribution_id=2)
        generated=SpinDensityExporter._analytic_top_decay_lines(variant,layout,'1D-3')
        program='''PROGRAM CHECK_RESCUE
USE TOP_DECAY_VIRTUAL_DISPATCH
IMPLICIT NONE
REAL*8 MDL_MT,MDL_MW,MU_R,G,PBMAG
REAL*8 SDM_INSERTION_P(0:3,3),TDV_VALIDATION_P(0:3,4)
REAL*8 TDV_PRECISION_ASKED,SDM_PRECISION
COMPLEX*16 GC_11,TDV_ANALYTIC_RHO(3,2,2),SDM_INSERTION_RHO(3,2,2)
INTEGER CTMODEINIT,CTMODERUN,TDV_SAVED_CT_MODE_INIT,TDV_SAVED_CT_MODE_RUN
INTEGER SDM_RET_CODE,CALLS
LOGICAL TDV_ANALYTIC_AVAILABLE,TDV_NEEDS_MADLOOP
CHARACTER(16) MODE
CALL GET_COMMAND_ARGUMENT(1,MODE)
MDL_MT=172.5D0
MDL_MW=80.385D0
MU_R=MDL_MT
G=1D0
GC_11=(0D0,0.46D0)
PBMAG=(MDL_MT**2-MDL_MW**2)/(2D0*MDL_MT)
SDM_INSERTION_P(:,1)=[MDL_MT,0D0,0D0,0D0]
SDM_INSERTION_P(:,2)=[PBMAG,0D0,0D0,PBMAG]
SDM_INSERTION_P(:,3)=[MDL_MT-PBMAG,0D0,0D0,-PBMAG]
CTMODEINIT=0
CTMODERUN=0
CALLS=0
'''.splitlines()+generated+'''
IF (CTMODEINIT.NE.1.OR.CTMODERUN.NE.-1) STOP 31
IF (MODE.EQ.'good'.AND.CALLS.NE.1) STOP 32
IF (MODE.NE.'good'.AND.CALLS.NE.2) STOP 33
WRITE(*,*) 'GENERATED_REFERENCE_RESCUE_PASS'
CONTAINS
SUBROUTINE MOCK_MADLOOP(P,RHO,REQUESTED,PRECISION,RET_CODE)
REAL*8 P(0:3,3),REQUESTED,PRECISION
COMPLEX*16 RHO(3,2,2)
INTEGER RET_CODE
CALLS=CALLS+1
IF (CALLS.EQ.1) THEN
C The first call initializes MadLoop's card before settings can be saved.
CTMODEINIT=1
CTMODERUN=-1
ELSE
IF (CTMODEINIT.NE.4.OR.CTMODERUN.NE.4) STOP 34
ENDIF
IF (ABS(REQUESTED-1D-10).GT.1D-25) STOP 35
RHO=TDV_ANALYTIC_RHO
IF ((MODE.NE.'good'.AND.CALLS.EQ.1).OR.MODE.EQ.'bad_density') THEN
RHO(1,1,2)=RHO(1,1,2)+(0D0,0.01D0)*MAXVAL(ABS(RHO))
ENDIF
PRECISION=0D0
RET_CODE=140
END SUBROUTINE
END PROGRAM
'''.splitlines()
        generated_driver=cls.directory/'generated.f'
        generated_driver.write_text('\n'.join(
            line if line.startswith('C ') else
            '     $ '+line.lstrip()[1:].lstrip() if line.lstrip().startswith('$') else
            '      '+line.strip() for line in program)+'\n')
        cls.generated_executable=cls.directory/'generated'
        result=subprocess.run([compiler,'-O0','-g','-fcheck=all','-ffixed-line-length-132',
            '-ffree-line-length-none',str(stub),str(source/'top_decay_virtual_cdr.f90'),
            str(source/'top_decay_virtual_dispatch.f90'),str(Path(MG5DIR)/'HELAS/vxxxxx.F'),
            str(generated_driver),'-o',str(cls.generated_executable)],
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

    def test_generated_rescue_is_conditional_and_restores_modes(self):
        for mode in ('good','rescue'):
            with self.subTest(mode=mode):
                result=subprocess.run([str(self.generated_executable),mode],text=True,capture_output=True)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                self.assertIn('GENERATED_REFERENCE_RESCUE_PASS',result.stdout)
                self.assertEqual('uniform quadruple-precision' in result.stdout,mode=='rescue')

    def test_failed_quadruple_reference_is_still_fatal(self):
        result=subprocess.run([str(self.generated_executable),'bad_density'],text=True,capture_output=True)
        self.assertNotEqual(result.returncode,0)
        self.assertIn('analytic top-decay virtual validation failed',result.stdout)
        self.assertNotIn('GENERATED_REFERENCE_RESCUE_PASS',result.stdout)
