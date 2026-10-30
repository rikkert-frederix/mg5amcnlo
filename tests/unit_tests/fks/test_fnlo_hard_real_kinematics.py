"""Numerical measures and inverse-map checks for production hard channels."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR


class TestHardRealKinematics(unittest.TestCase):
    def test_measures_and_inverse_maps(self):
        compiler = shutil.which('gfortran')
        if compiler is None:
            self.skipTest('gfortran is required')
        root = Path(MG5DIR)
        with tempfile.TemporaryDirectory(prefix='fnlo_hard_real_') as directory:
            executable = Path(directory) / 'check'
            result = subprocess.run([
                compiler, '-O0', '-g', '-fcheck=all',
                '-ffpe-trap=invalid,zero,overflow', '-o', str(executable),
                str(root / 'Template/fNLO/SubProcesses/hard_real_kinematics.f90'),
                str(root / 'tests/input_files/fks_decay/hard_real_kinematics_checks.f90'),
            ], cwd=directory, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            result = subprocess.run([str(executable)], cwd=directory,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn('PASS: real diagram measures', result.stdout)
