"""The actual leaf refresh preserves normalization and fixes spin bases."""

from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR


class TestRefreshMultiplicativeLOFrames(unittest.TestCase):
    def test_polarized_vector_decay_in_changed_leaf_frame(self):
        compiler = shutil.which('gfortran')
        if compiler is None:
            self.skipTest('gfortran is required')
        root = Path(MG5DIR)
        driver = (root/'Template/fNLO/SubProcesses/driver_mintFO.f90').read_text()
        start = driver.index('  subroutine refresh_multiplicative_lo_frames(')
        end = driver.index('  end subroutine refresh_multiplicative_lo_frames', start)
        end = driver.index('\n', end)
        routine = driver[start:end]
        fixture = (root/'tests/input_files/fks_decay/refresh_lo_frame_checks.f90').read_text()
        self.assertEqual(fixture.count('! INSERT_ACTUAL_REFRESH'), 1)
        with tempfile.TemporaryDirectory(prefix='refresh_lo_frames_') as directory:
            directory = Path(directory)
            source = directory/'check.f90'
            source.write_text(fixture.replace('! INSERT_ACTUAL_REFRESH', routine))
            executable = directory/'check'
            result = subprocess.run([compiler, '-O0', '-g', '-fcheck=all',
                '-ffree-line-length-none', '-ffpe-trap=invalid,zero,overflow',
                str(source), '-o', str(executable)], cwd=directory,
                capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            result = subprocess.run([str(executable)], cwd=directory,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout+result.stderr)
            self.assertIn('PASS: physical spin-one LO leaf frame refresh', result.stdout)
