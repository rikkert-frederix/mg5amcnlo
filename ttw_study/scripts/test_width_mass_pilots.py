import unittest

from run_sampler_pilots import run_settings
from run_width_mass_pilots import cases, export_key, process_path, run_name


class TestWidthMassPilots(unittest.TestCase):
    def test_successor_seed_range_does_not_reuse_failed_campaign(self):
        self.assertEqual({r['seed'] for r in cases(60017)},set(range(60017,60033)))
        self.assertFalse({r['seed'] for r in cases()} & {r['seed'] for r in cases(60017)})
        for invalid in (0,-1,900000000,1.5):
            with self.assertRaises(ValueError):
                cases(invalid)

    def test_fresh_complete_companions_with_unique_seeds(self):
        rows=cases()
        self.assertEqual(len(rows),16)
        self.assertEqual({row['seed'] for row in rows},set(range(60001,60017)))
        self.assertEqual(len({export_key(row) for row in rows}),8)
        self.assertEqual(len({run_name(row) for row in rows}),16)
        for row in rows:
            self.assertEqual(row['production_sampling'],'w-current' if row['w_treatment']=='all-bw' else 'flat')
            settings=run_settings(process_path(row,'test'),row,.03)
            self.assertIsNone(settings.grid_reference)
            self.assertFalse(settings.main_run)
            self.assertEqual(settings.accuracy,.03)
            self.assertEqual(settings.production_scale,'core-w-ht-half')
            self.assertEqual(run_name(row).endswith('_mb4p8'),row['decay_bottom_mass']==4.8)
        for charge in ('plus','minus'):
            for mode in ('onshell','all-bw'):
                self.assertEqual({row['variant'] for row in rows if row['charge']==charge
                                  and row['w_treatment']==mode and row['decay_bottom_mass']==4.8},{'S','Pi'})

    def test_massless_references_exist_before_massive_exports(self):
        available={('plus','onshell'),('plus','all-bw')}
        for row in cases():
            key=(row['charge'],row['w_treatment'])
            if row['decay_bottom_mass']:
                self.assertIn(key,available)
            else:
                available.add(key)


if __name__=='__main__':
    unittest.main()
