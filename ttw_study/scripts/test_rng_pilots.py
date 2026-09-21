import copy
import unittest

from run_rng_pilots import cases, reference_gate
from run_sampler_pilots import run_name, run_settings


class TestRNGPilots(unittest.TestCase):
    def test_successor_uses_distinct_explicit_seed_range(self):
        self.assertEqual({r['seed'] for r in cases(57005)},set(range(57005,57009)))
        self.assertFalse({r['seed'] for r in cases()} & {r['seed'] for r in cases(57005)})
        for invalid in (0,-1,900000000,1.5):
            with self.assertRaises(ValueError):
                cases(invalid)

    def test_four_fresh_independent_cases_without_grid_transfer(self):
        rows=cases()
        self.assertEqual(len(rows),4)
        self.assertEqual(len({run_name(r) for r in rows}),4)
        self.assertEqual({r['seed'] for r in rows},set(range(57001,57005)))
        for mode in ('all-bw','onshell'):
            self.assertEqual({r['variant'] for r in rows if r['w_treatment']==mode},{'S','Pi'})
        for row in rows:
            self.assertIsNone(run_settings('fresh',row,.03).grid_reference)
            self.assertEqual(row['production_sampling'],'w-current' if row['w_treatment']=='all-bw' else 'flat')

    def test_unresolved_reference_stops_further_launches(self):
        comparison=dict(status='independent LO current normalization pilot; inspect convergence',
                        associated_branching_factor_applied=False,decayed_pb=[.001]*183,
                        expected_pb=[.001]*183,nominal_combined_mc_error_pb=.000005,nominal_pull=1.)
        original=copy.deepcopy(comparison)
        self.assertEqual(reference_gate(comparison)['nominal_pull'],1.)
        self.assertEqual(comparison,original)
        for changes in (dict(nominal_pull=3.1),dict(nominal_pull=float('nan')),
                        dict(nominal_combined_mc_error_pb=.0001),
                        dict(associated_branching_factor_applied=True),dict(decayed_pb=[.001]),
                        dict(expected_pb=[0.]*183)):
            with self.assertRaises(ValueError):
                reference_gate(dict(comparison,**changes))


if __name__=='__main__':
    unittest.main()
