"""Production hard channels are added once to new fixed-order NLO jobs."""
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest

from madgraph.interface.amcatnlo_run_interface import aMCatNLOCmd, aMCatNLOError


class TestProductionHardChannels(unittest.TestCase):
    def test_native_jobs_include_each_real_diagram_once(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'Cards').mkdir()
            process = root / 'SubProcesses' / 'P0_test'
            process.mkdir(parents=True)
            card = root / 'Cards' / 'decay_card.dat'
            card.write_text('NONE = production_hard_split\n')
            (process / 'channels.txt').write_text('1.1 1.2 2.1 2.2\n')
            (process / 'hard_symfact.dat').write_text(
                ''.join('%d.3 1\n' % i for i in range(1, 37)))
            context = SimpleNamespace(me_dir=str(root), run_card={
                'npoints_FO_grid': 1000, 'niters_FO_grid': 3})
            context.production_hard_channels = lambda p, mode: \
                aMCatNLOCmd.production_hard_channels(context, p, mode)

            def jobs(mode, accuracy=-1):
                result, collect, stage = aMCatNLOCmd.create_jobs_to_run(
                    context, {'only_generation': False}, ['P0_test'], accuracy,
                    'all' if mode == 'NLO' else 'born', -1, mode)
                self.assertEqual(result, collect)
                self.assertEqual(stage, -1)
                self.assertTrue(all(job['nchans'] <= 20 for job in result))
                self.assertEqual([job['channel'] for job in result],
                                 [str(i + 1) for i in range(len(result))])
                for job in result:
                    hard = [token.endswith('.3') for token in job['configs'].split()]
                    self.assertTrue(all(hard) or not any(hard), job['configs'])
                return result

            original_jobs = jobs('NLO')
            self.assertEqual([job['nchans'] for job in original_jobs], [4])
            original = original_jobs[0]['configs'].split()
            self.assertEqual(original, ['1.1', '1.2', '2.1', '2.2'])
            card.write_text("'real_diagrams' = production_hard_split ! sampling\n")
            for accuracy in (-1, 0.005):
                enabled_jobs = jobs('NLO', accuracy)
                self.assertEqual([job['nchans'] for job in enabled_jobs], [4, 18, 18])
                self.assertEqual(enabled_jobs[0]['configs'].split(), original)
                enabled = [token for job in enabled_jobs for token in job['configs'].split()]
                self.assertEqual(enabled, original + ['%d.3' % i for i in range(1, 37)])
                self.assertEqual(len(enabled), len(set(enabled)))
            self.assertEqual(jobs('LO')[0]['configs'].split(), original)

            # Neither family needs to fill the other family's last batch.
            (process / 'channels.txt').write_text(
                ' '.join('%d.1' % i for i in range(1, 26)) + '\n')
            large = jobs('NLO')
            self.assertEqual([job['nchans'] for job in large], [13, 12, 18, 18])
            (process / 'channels.txt').write_text('\n')
            hard_only = jobs('NLO')
            self.assertEqual([job['nchans'] for job in hard_only], [18, 18])

    def test_old_exports_fail_before_a_partial_split_is_launched(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / 'Cards').mkdir()
            card = root / 'Cards' / 'decay_card.dat'
            context = SimpleNamespace(me_dir=str(root))
            card.write_text('REAL_DIAGRAMS = production_hard_split\n')
            with self.assertRaisesRegex(aMCatNLOError, 'Re-export'):
                aMCatNLOCmd.production_hard_channels(context, 'P0_test', 'NLO')
            card.write_text('INVALID = production_hard_split\n')
            with self.assertRaisesRegex(aMCatNLOError, 'NONE or REAL_DIAGRAMS'):
                aMCatNLOCmd.production_hard_channels(context, 'P0_test', 'NLO')


if __name__ == '__main__':
    unittest.main()
