"""Actual recursive materializer exposes the frames of assembled decay blocks."""
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR


class TestLeafBlockFrames(unittest.TestCase):
    def test_recursive_recoil_and_optional_api(self):
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is required')
        root = Path(MG5DIR)
        fixture = root / 'tests/input_files/fks_decay'
        # Compile the actual boost and materialization code; only the generated
        # process layout and workspace storage are replaced by a three-node tree.
        source = (fixture / 'leaf_block_frames_checks.f90').read_text()
        declarations, program = source.split('program check_leaf_frames', 1)
        with tempfile.TemporaryDirectory(prefix='fnlo_leaf_frames_') as build:
            build = Path(build)
            (build / 'layout.f90').write_text(declarations)
            (build / 'check.f90').write_text('program check_leaf_frames' + program)
            executable = build / 'check'
            result = subprocess.run([
                compiler, '-O0', '-g', '-fcheck=all',
                '-ffpe-trap=invalid,zero,overflow', '-o', str(executable),
                str(fixture / 'phase_space_test_dimensions.f90'),
                str(root / 'Template/fNLO/Source/kin_functions.f90'),
                str(root / 'Template/fNLO/SubProcesses/phase_space_kinematics.f90'),
                str(root / 'Template/fNLO/SubProcesses/factorized_block_kinematics.f90'),
                str(build / 'layout.f90'),
                str(root / 'Template/fNLO/SubProcesses/multiplicative_event_materialization.f90'),
                str(build / 'check.f90')], cwd=build, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            result = subprocess.run([str(executable)], cwd=build, capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn('PASS: recursive leaf block helicity frames', result.stdout)
