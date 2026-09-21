"""Test the actual histogram accumulator, including sparse iterations."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR


class TestFNLOHwU(unittest.TestCase):
    def test_empty_iterations_and_zero_variance_are_not_missing_data(self):
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is required')
        root = Path(MG5DIR)
        with tempfile.TemporaryDirectory(prefix='fnlo_hwu_') as build:
            executable = Path(build) / 'check'
            result = subprocess.run([
                compiler, '-O0', '-g', '-fcheck=all',
                '-ffpe-trap=invalid,zero,overflow', '-o', str(executable),
                str(root / 'Template/fNLO/FixedOrderAnalysis/HwU.f90'),
                str(root / 'tests/input_files/fks_decay/hwu_iteration_checks.f90')],
                cwd=build, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            result = subprocess.run([str(executable)], cwd=build,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        rows = []
        control_rows = []
        current = None
        for line in result.stdout.splitlines():
            if line.startswith('<histogram>'):
                current = rows if 'iteration averaging' in line else control_rows
            elif line.startswith('<'):
                current = None
            elif current is not None and line.strip():
                current.append([float(value) for value in line.split()])
        self.assertEqual(len(rows), 16, result.stdout)
        self.assertEqual(len(control_rows), 24, result.stdout)
        self.assertTrue(all(row[2:] == [0., 0., 0.] for row in control_rows))
        first = [(2/3, 20/3, 4/3), (1/3, 10/3, 2/3),
                 (4/3, 0., 4.), (0., 0., 0.)]
        expected = first + [tuple(value*9/11 for value in row) for row in first]
        for row, reference in zip(rows, expected*2):
            for actual, target in zip(row[2:], reference):
                # HwU output intentionally stores seven significant figures.
                self.assertAlmostEqual(actual, target, delta=max(1.e-8, abs(target)*5.1e-7))
