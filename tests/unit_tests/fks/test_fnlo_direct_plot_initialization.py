"""Direct real leaves must book histograms without an n-body evaluation."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR


class TestDirectPlotInitialization(unittest.TestCase):
    def test_first_multiplicative_leaf_initializes_once(self):
        compiler = shutil.which('gfortran')
        if compiler is None:
            self.skipTest('gfortran is required')
        root = Path(MG5DIR)
        source = (root / 'Template/fNLO/SubProcesses/madfks_plot.f90').read_text()
        start = source.index('  subroutine outfun_multiplicative_impl(')
        end_marker = '  end subroutine outfun_multiplicative_impl'
        end = source.index(end_marker, start) + len(end_marker)
        init_start = source.index('  subroutine initplot_impl()')
        init_end_marker = '  end subroutine initplot_impl'
        init_end = source.index(init_end_marker, init_start) + len(init_end_marker)
        # Exercise the actual plotting entry and actual HwU implementation;
        # the fixture analysis books one bin and fills its central weight.
        harness = '''
module direct_plot_fixture
  use boostwdir2_module, only: boostwdir2
  use HwU_module
  implicit none
  integer, parameter :: nincoming=2
  logical :: plot_initialized=.false.
  integer :: initialization_count=0
  logical :: do_rwgt_scale=.false.,do_rwgt_decay_scale=.false.,do_rwgt_pdf=.false.
  integer :: dyn_scale(0:1)=(/1,0/),lhaPDFid(0:1)=(/0,0/)
  integer :: nmemPDF(1)=0
  integer, parameter :: maxPDFs=1
  logical :: lpdfvar(1)=.false.
  character(len=32) :: LHAPDFsetname(1)='unused'
  double precision :: xsecPDFr_acc(0:1,1)=0d0
  double precision, allocatable :: xsec_scale_points(:,:)
contains
''' + source[start:end] + '\n' + source[init_start:init_end] + '''
  subroutine validate_process_dimensions()
  end subroutine
  subroutine analysis_begin(nwgt,labels)
    integer, intent(in) :: nwgt
    character(len=*), intent(in) :: labels(*)
    call HwU_inithist(nwgt,labels)
    call HwU_book(1,'direct hard leaf',1,0d0,1d0)
    initialization_count=initialization_count+1
  end subroutine
  integer function fnlo_scale_max_point_count()
    fnlo_scale_max_point_count=1
  end function
  integer function fnlo_scale_point_count(scale)
    integer, intent(in) :: scale
    fnlo_scale_point_count=0
  end function
  subroutine fnlo_scale_point_label(scale,point,label)
    integer, intent(in) :: scale,point
    character(len=*), intent(out) :: label
    label='unused'
  end subroutine
  subroutine initpdfsetbynamem(set,name)
    integer, intent(in) :: set
    character(len=*), intent(in) :: name
    stop 12
  end subroutine
  subroutine numberPDFm(set,count)
    integer, intent(in) :: set
    integer, intent(out) :: count
    stop 13
  end subroutine
  subroutine InitPDFm(set,member)
    integer, intent(in) :: set,member
    stop 14
  end subroutine
  subroutine analysis_fill_multiplicative(p,count,status,pdg,weights,body)
    integer, intent(in) :: count,status(count),pdg(count),body
    double precision, intent(in) :: p(0:4,count),weights(*)
    if (count/=3 .or. body/=1) stop 10
    call HwU_fill(1,0.5d0,weights)
  end subroutine
  subroutine accumulate_analysis_weights(weights)
    double precision, intent(in) :: weights(*)
  end subroutine
end module
program direct_plot_check
  use direct_plot_fixture
  implicit none
  double precision :: p(0:3,3),weights(1),values(2)
  integer :: status(3),pdg(3)
  p=0d0
  p(:,1)=(/1d0,0d0,0d0,1d0/)
  p(:,2)=(/1d0,0d0,0d0,-1d0/)
  p(:,3)=(/2d0,0d0,0d0,0d0/)
  status=(/-1,-1,1/); pdg=(/1,-1,21/)
  weights=2d0
  call outfun_multiplicative_impl(p,3,0d0,weights,pdg,status,1)
  call HwU_add_points()
  ! A later first n-body prefactor still calls initplot unconditionally.
  ! It must preserve the already accumulated direct-hard contribution.
  call initplot_impl()
  weights=4d0
  call outfun_multiplicative_impl(p,3,0d0,weights,pdg,status,1)
  call HwU_add_points()
  if (initialization_count/=1) stop 11
  values=1d0
  call HwU_accum_iter(.true.,2,values)
  open(17,file='histogram.dat',status='replace')
  call HwU_output(17,1d0)
  close(17)
end program
'''
        with tempfile.TemporaryDirectory(prefix='fnlo_direct_plot_') as directory:
            directory = Path(directory)
            fixture = directory / 'check.f90'
            fixture.write_text(harness)
            executable = directory / 'check'
            command = [compiler, '-O0', '-g', '-fcheck=all',
                       '-ffpe-trap=invalid,zero,overflow', '-o', str(executable),
                       str(root / 'Template/fNLO/SubProcesses/boostwdir2.f90'),
                       str(root / 'Template/fNLO/FixedOrderAnalysis/HwU.f90'),
                       str(fixture)]
            result = subprocess.run(command, cwd=directory,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            result = subprocess.run([str(executable)], cwd=directory,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            bins = [line.split() for line in (directory / 'histogram.dat').read_text().splitlines()
                    if line.lstrip().startswith('0.')]
            self.assertEqual(len(bins), 1)
            self.assertAlmostEqual(float(bins[0][2]), 3.0)
