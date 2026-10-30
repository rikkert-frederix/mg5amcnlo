"""IR-safe hard separation, including independent phase-space measures."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR
import tests.unit_tests.fks.test_fnlo_decay_card as card_fixture


class TestProductionHardSplit(unittest.TestCase):
    def test_raw_hard_matrix_and_density(self):
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is required')
        root = Path(MG5DIR)
        template = root/'Template/fNLO/SubProcesses'
        original = (template/'fks_singular.f90').read_text()
        routines = []
        for name in ('evaluate_hard_real_matrix', 'sreal'):
            start = original.index('  subroutine '+name+'(')
            end = original.index('  end subroutine '+name, start)
            end = original.index('\n', end)
            routines.append(original[start:end])
        dimensions = '''module process_dimensions
  implicit none
  integer,parameter :: nexternal=3, nincoming=2, amp_split_size=1, nsplitorders=1, qcd_pos=1
  integer,parameter :: amp_split_orders(1,1)=reshape([2],[1,1])
end module
'''
        harness = '''module hard_matrix_fixture
  use process_dimensions
  use spin_density_fks_matrices
  use kin_functions_module, only: dot => dot_impl
  implicit none
  logical :: softtest=.false.,colltest=.false.,need_color_links=.true.
  integer :: i_fks=3,j_fks=1,nfksprocess=1,limit_calls=0
  double precision :: test_xi=0.2d0,test_y=-0.4d0,amp_split(1)
contains
  function fixture_density_shape() result(shape)
    complex(kind=8) :: shape(2,2)
    shape=reshape([(2d0,0d0),(0d0,-1d0),(0d0,1d0),(3d0,0d0)],[2,2])
  end function
  subroutine evaluate_real_matrix(event_slot,weight)
    integer,intent(in) :: event_slot
    double precision,intent(out) :: weight
    complex(kind=8) :: density(2,2,2)
    weight=3.1d0/(test_xi**2*(1d0-test_y))
    amp_split(1)=weight
    density=(0d0,0d0)
    density(1,:,:)=weight*fixture_density_shape()
    call set_spin_density_real_matrix(1,density)
  end subroutine
  logical function has_nlo_decay()
    has_nlo_decay=.false.
  end function
  integer function nlo_decay_fks_i(configuration)
    integer,intent(in) :: configuration
    nlo_decay_fks_i=3
  end function
  integer function nlo_decay_fks_j(configuration)
    integer,intent(in) :: configuration
    nlo_decay_fks_j=1
  end function
  subroutine select_kernel_event(event_slot,momenta,masses)
    integer,intent(in) :: event_slot
    double precision,intent(out) :: momenta(0:3,3),masses(3)
    momenta=0d0; masses=0d0
    momenta(:,1)=[500d0,0d0,0d0,500d0]
    momenta(:,2)=[500d0,0d0,0d0,-500d0]
  end subroutine
  double precision function fks_subtraction_shat(event_slot)
    integer,intent(in) :: event_slot
    fks_subtraction_shat=1d6
  end function
  subroutine sbornsoft(event_slot,momenta,xi,y,q,weight)
    integer,intent(in) :: event_slot
    double precision,intent(in) :: momenta(0:3,3),xi,y,q
    double precision,intent(out) :: weight
    limit_calls=limit_calls+1; weight=777d0; amp_split=777d0
  end subroutine
  subroutine sborncol_isr(momenta,xi,y,shat,weight)
    double precision,intent(in) :: momenta(0:3,3),xi,y,shat
    double precision,intent(out) :: weight
    limit_calls=limit_calls+1; weight=777d0; amp_split=777d0
  end subroutine
  subroutine sborncol_fsr(momenta,y,shat,weight)
    double precision,intent(in) :: momenta(0:3,3),y,shat
    double precision,intent(out) :: weight
    limit_calls=limit_calls+1; weight=777d0; amp_split=777d0
  end subroutine
'''+ '\n'.join(routines)+'\nend module hard_matrix_fixture\n'
        with tempfile.TemporaryDirectory(prefix='hard_raw_matrix_') as directory:
            dimensions_file = Path(directory)/'dimensions.f90'
            dimensions_file.write_text(dimensions)
            harness_file = Path(directory)/'harness.f90'
            harness_file.write_text(harness)
            executable = Path(directory)/'check'
            result = subprocess.run([
                compiler, '-O0', '-g', '-fcheck=all', '-ffree-line-length-none',
                '-ffpe-trap=invalid,zero,overflow', '-o', str(executable),
                str(dimensions_file), str(root/'Template/fNLO/Source/kin_functions.f90'),
                str(template/'spin_density_fks_matrices.f90'),
                str(template/'production_hard_split.f90'), str(harness_file),
                str(root/'tests/input_files/fks_decay/hard_raw_matrix_checks.f90')],
                cwd=directory, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            result = subprocess.run([str(executable)], cwd=directory,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertIn('native real scalar/density agree away from IR', result.stdout)
            self.assertIn('raw hard IR weight agrees with quadruple precision', result.stdout)

    def test_massive_fsr_inverse_against_native_forward_map(self):
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is required')
        root = Path(MG5DIR)
        template = root/'Template/fNLO/SubProcesses'
        fixture = root/'tests/input_files/fks_decay'
        source = (template/'genps_fks.f90').read_text()
        routines = []
        for name in ('generate_momenta_massive_final', 'getangles'):
            start = source.index('  subroutine '+name+'(')
            end = source.index('  end subroutine '+name, start)
            end = source.index('\n', end)
            routines.append(source[start:end])
        # Compile the actual forward kernel, rather than another transcription
        # of its algebra, against a minimal version of its module globals.
        native_source = '''module native_massive_fks
  use process_dimensions, only: nexternal, nincoming
  use kin_functions_module, only: rho => rho_impl
  use phase_space_kinematics, only: rotate_invar
  use boostwdir2_module, only: boostwdir2_in_place
  implicit none
  integer, parameter :: real_event=3, soft_counterevent=0
  integer, parameter :: collinear_counterevent=1, soft_collinear_counterevent=2
  logical :: softtest=.false., colltest=.false.
  double precision :: xi_i_fks_fix=-2d0, y_ij_fks_fix=-2d0
  double precision :: massive_xjac_cache=1d0
contains
'''+ '\n'.join(routines)+'\nend module native_massive_fks\n'
        with tempfile.TemporaryDirectory(prefix='native_fsr_inverse_') as directory:
            native_file = Path(directory)/'native_massive_fks.f90'
            native_file.write_text(native_source)
            executable = Path(directory)/'check'
            result = subprocess.run([
                compiler, '-O0', '-g', '-fcheck=all', '-ffree-line-length-none',
                '-ffpe-trap=invalid,zero,overflow', '-o', str(executable),
                str(fixture/'phase_space_test_dimensions.f90'),
                str(root/'Template/fNLO/Source/kin_functions.f90'),
                str(template/'phase_space_kinematics.f90'),
                str(template/'boostwdir2.f90'),
                str(template/'hard_real_kinematics.f90'), str(native_file),
                str(fixture/'hard_fsr_native_roundtrip_checks.f90')],
                cwd=directory, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            result = subprocess.run([str(executable)], cwd=directory,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertIn('PASS: 60 native massive FSR inverse round trips', result.stdout)

    def test_card_writer_and_runtime(self):
        card_fixture.TestFNLODecayCard.setUpClass()
        fixture = card_fixture.TestFNLODecayCard()
        try:
            base = fixture.card()
            self.assertEqual(base, fixture.card(production_hard_split='NONE'))
            self.assertIn('NONE = production_hard_split ! NONE or REAL_DIAGRAMS.', base)
            self.assertIn('1.0000000000000001e-01 = production_hard_scale', base)
            default = fixture.run_card(base, 'hard_proposal')['HARD_PROPOSAL']
            self.assertEqual(default[0], 'F')
            self.assertEqual(float(default[1]), 0.1)
            for scale in (0.03, 0.1, 0.3):
                card = fixture.card(production_hard_split='REAL_DIAGRAMS',
                                    production_hard_scale=scale)
                row = fixture.run_card(card, 'hard_proposal')['HARD_PROPOSAL']
                self.assertEqual(row[0], 'T')
                self.assertEqual(float(row[1]), scale)
                self.assertEqual(fixture.run_card(base), fixture.run_card(card))
            for options in ({'production_hard_split': 'BAD'},
                            {'production_hard_scale': 0.},
                            {'production_hard_scale': -1.},
                            {'production_hard_scale': float('inf')},
                            {'production_hard_scale': float('nan')},
                            {'production_hard_split': 'REAL_DIAGRAMS',
                             'production_order': 'LO'}):
                with self.assertRaises(ValueError):
                    fixture.card(**options)
            card = fixture.card(production_hard_split='REAL_DIAGRAMS')
            for text, message in [
                    (card+'\nNONE = production_hard_split', 'duplicate PRODUCTION_HARD_SPLIT'),
                    (card+'\n0.03 = production_hard_scale', 'duplicate PRODUCTION_HARD_SCALE'),
                    (card.replace('REAL_DIAGRAMS =', 'BAD ='), 'NONE or REAL_DIAGRAMS'),
                    (card.replace('1.0000000000000001e-01 = production_hard_scale',
                                  '0. = production_hard_scale'), 'positive and finite'),
                    (card.replace('NLO = production_order', 'LO = production_order'),
                     'requires NLO production')]:
                fixture.run_card(text, 'hard_proposal', error=message)
        finally:
            card_fixture.TestFNLODecayCard.doClassCleanups()

    def test_partition_and_subtracted_density_integration(self):
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is required')
        root = Path(MG5DIR)
        with tempfile.TemporaryDirectory(prefix='production_hard_split_') as directory:
            executable = Path(directory)/'check'
            result = subprocess.run([
                compiler, '-O0', '-g', '-fcheck=all',
                '-ffpe-trap=invalid,zero,overflow', '-o', str(executable),
                str(root/'Template/fNLO/SubProcesses/production_hard_split.f90'),
                str(root/'tests/input_files/fks_decay/production_hard_split_checks.f90')],
                cwd=directory, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            result = subprocess.run([str(executable)], cwd=directory,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertIn('PASS: IR limits, massive sister and Lorentz invariance', result.stdout)
            self.assertIn('PASS: subtracted scalar/density closure with independent proposals', result.stdout)
            self.assertIn('PASS: disabled split', result.stdout)
            for mode, message in [('bad_scale', 'finite and positive'),
                                  ('nan_scale', 'finite and positive'),
                                  ('bad_channel', 'no real diagram channel')]:
                result = subprocess.run([str(executable), mode], cwd=directory,
                                        capture_output=True, text=True)
                self.assertNotEqual(result.returncode, 0, result.stdout+result.stderr)
                self.assertIn(message, result.stdout)


if __name__ == '__main__':
    unittest.main()
