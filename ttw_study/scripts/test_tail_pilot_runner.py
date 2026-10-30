"""Regression checks for bounded process execution and replay isolation."""
import argparse
import asyncio
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import run_tail_pilot as pilot


class TestTailPilotRunner(unittest.TestCase):
    def fixture(self, root, binary):
        process = root/'process'
        parent = process/'SubProcesses'/pilot.SUBPROCESS
        parent.mkdir(parents=True)
        (process/'Cards').mkdir()
        (process/'Cards/decay_card.dat').write_text(
            'MULTIPLICATIVE = nlo_decay_combination\nNONE = decay_scale_variation_mode\n')
        for name in pilot.RUNTIME_FILES:
            (parent/name).write_text('runtime data\n')
        proposal = process/'pilot/proposal'
        proposal.mkdir(parents=True)
        for name in ('mint_grids', 'grid.MC_integer'):
            (proposal/name).write_text('unchanged initial proposal\n')
        executable = parent/'madevent_mintFO'
        executable.write_text('#!'+sys.executable+'\n'+binary)
        executable.chmod(0o755)
        return process, proposal

    def test_fixed_count_replay_and_variant_cards_are_isolated(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            process, proposal = self.fixture(root, 'from pathlib import Path\n'
                'Path("res.dat").write_text("finished")\n')
            jobs = [pilot.make_worker(process, proposal, 'test', variant, 0, seed, 32)
                for variant, seed in (('S', 910001), ('Pi', 910002))]
            replay = pilot.make_worker(process, proposal, 'replay', 'S', 0, 910001, 32)
            for name in ('randinit', 'input_app.txt', 'mint_grids', 'grid.MC_integer', 'decay_card.dat'):
                self.assertEqual((Path(jobs[0]['directory'])/name).read_bytes(),
                                 (Path(replay['directory'])/name).read_bytes())
                self.assertFalse((Path(replay['directory'])/name).is_symlink())
            self.assertIn('ADDITIVE', (Path(jobs[0]['directory'])/'decay_card.dat').read_text())
            self.assertIn('MULTIPLICATIVE', (process/'Cards/decay_card.dat').read_text())
            self.assertIn('NITERATIONS = 1', (Path(jobs[0]['directory'])/'input_app.txt').read_text())
            runtime = pilot.runtime_provenance(process)
            pilot.write_worker_manifest(jobs[0], runtime, {'description':'saved proposal'})
            manifest = json.loads(Path(jobs[0]['manifest']).read_text())
            self.assertEqual(manifest['runtime']['binary_sha256'],
                pilot.digest(process/'SubProcesses'/pilot.SUBPROCESS/'madevent_mintFO'))
            self.assertEqual(manifest['inputs']['seed'], 910001)
            self.assertEqual(manifest['input_sha256']['mint_grids'],
                pilot.digest(proposal/'mint_grids'))
            execution = dict(jobs=jobs)
            data = dict(process=str(process), executions=[execution])
            with patch.object(pilot, 'lhapdf_environment', return_value={}):
                asyncio.run(pilot.run_workers(root/'record.json', data, execution, 2, 10))
            self.assertTrue(all(job['status'] == 'finished' for job in jobs))
            self.assertEqual((proposal/'mint_grids').read_text(), 'unchanged initial proposal\n')

    def test_timeout_is_saved_and_process_group_terminated(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            process, proposal = self.fixture(root, 'import time\ntime.sleep(30)\n')
            job = pilot.make_worker(process, proposal, 'timeout', 'S', 0, 910001, 32)
            execution = dict(jobs=[job])
            data = dict(process=str(process), executions=[execution])
            with patch.object(pilot, 'lhapdf_environment', return_value={}):
                asyncio.run(pilot.run_workers(root/'record.json', data, execution, 1, .05))
            saved = json.loads((root/'record.json').read_text())['executions'][0]['jobs'][0]
            self.assertEqual(saved['status'], 'timed out')
            self.assertLess(saved['returncode'], 0)
            self.assertLess(saved['elapsed_s'], 5)

    def test_virtual_floor_experiment_does_not_modify_shared_grid_or_parameters(self):
        with tempfile.TemporaryDirectory() as directory:
            process, proposal = self.fixture(Path(directory), '')
            shared = process/'SubProcesses'/pilot.SUBPROCESS/'FKS_params.dat'
            original_parameters = '#VirtualFraction\n1.0d0\n#MinVirtualFraction\n0.005d0\n'
            shared.write_text(original_parameters)
            original_grid = ('AVE 0.0 0.1 0.2 0.3\nAVE 0.005 -0.2\n')*4
            (proposal/'mint_grids').write_text(original_grid)
            job = pilot.make_worker(process, proposal, 'floor', 'S', 0, 910001, 20000,
                iterations=4, restart=0, virtual_fraction=.05)
            worker = Path(job['directory'])
            self.assertEqual(shared.read_text(), original_parameters)
            self.assertEqual((proposal/'mint_grids').read_text(), original_grid)
            self.assertFalse((worker/'FKS_params.dat').is_symlink())
            lines = (worker/'mint_grids').read_text().splitlines()
            fractions = [float(line.split()[1]) for line in lines if len(line.split()) == 3]
            self.assertEqual(fractions, [.05]*4)
            self.assertEqual(sum(line == 'AVE 0.0 0.1 0.2 0.3' for line in lines), 4)
            self.assertEqual(job['previous_grid_virtual_fractions'], [.005]*4)
            self.assertIn('NITERATIONS = 4', (worker/'input_app.txt').read_text())
            self.assertIn('RESTART = 0', (worker/'input_app.txt').read_text())


if __name__ == '__main__':
    unittest.main()
