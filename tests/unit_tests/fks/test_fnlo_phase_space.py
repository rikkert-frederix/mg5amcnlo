"""Compile the actual phase-space kernel with floating-point traps enabled."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR


class TestFNLOPhaseSpace(unittest.TestCase):
    def test_bw_support_and_measures(self):
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is required')
        root = Path(MG5DIR)
        fixture = root / 'tests/input_files/fks_decay'
        with tempfile.TemporaryDirectory(prefix='fnlo_bw_support_') as build:
            executable = Path(build) / 'check'
            result = subprocess.run([
                compiler, '-O0', '-g', '-fcheck=all',
                '-ffpe-trap=invalid,zero,overflow', '-o', str(executable),
                str(fixture / 'phase_space_test_dimensions.f90'),
                str(root / 'Template/fNLO/Source/kin_functions.f90'),
                str(root / 'Template/fNLO/SubProcesses/phase_space_kinematics.f90'),
                str(root / 'Template/fNLO/SubProcesses/factorized_block_kinematics.f90'),
                str(fixture / 'bw_support_checks.f90')],
                cwd=build, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            result = subprocess.run([str(executable)], cwd=build,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn('CORE_SUPPORT_POINTS 15', result.stdout)
            self.assertIn('DECAY_SUPPORT_POINTS 20', result.stdout)
            self.assertIn('PASS: full-support', result.stdout)

    def test_lambda_boundaries(self):
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is required')
        root = Path(MG5DIR)
        fixture = root / 'tests/input_files/fks_decay'
        with tempfile.TemporaryDirectory(prefix='fnlo_lambda_') as build:
            executable = Path(build) / 'check'
            result = subprocess.run([
                compiler, '-O0', '-g', '-fcheck=all',
                '-ffpe-trap=invalid,zero,overflow', '-o', str(executable),
                str(fixture / 'phase_space_test_dimensions.f90'),
                str(root / 'Template/fNLO/Source/kin_functions.f90'),
                str(root / 'Template/fNLO/SubProcesses/phase_space_kinematics.f90'),
                str(fixture / 'phase_space_lambda_checks.f90')],
                cwd=build, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            result = subprocess.run([str(executable)], cwd=build,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertNotIn('NaN', result.stdout)
            self.assertIn('PASS: zero invariant', result.stdout)
