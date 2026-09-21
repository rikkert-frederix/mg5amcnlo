import collections
import itertools
import unittest

from run_main_campaign import cases


class TestMainCoverage(unittest.TestCase):
    def test_every_replica_has_every_physical_assignment_and_prescription_once(self):
        rows = cases()
        self.assertEqual(len(rows), 480)
        self.assertEqual(len({r['seed'] for r in rows}), 480)
        populations = collections.defaultdict(list)
        for row in rows:
            populations[row['replica'], row['charge'], row['variant']].append(tuple(row['flavours']))
        expected = set(itertools.product(('e', 'mu'), repeat=3))
        self.assertEqual(len(populations), 5*2*6)
        for flavours in populations.values():
            self.assertEqual(len(flavours), 8)
            self.assertEqual(set(flavours), expected)

    def test_insufficient_retraining_and_invalid_seed_ranges_are_rejected(self):
        for replicas, seed in [(4, 100001), (5, 0), (5, 900000000)]:
            with self.assertRaises(ValueError):
                cases(replicas, seed)


if __name__ == '__main__':
    unittest.main()
