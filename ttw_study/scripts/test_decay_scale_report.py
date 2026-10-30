"""Scientific reduction checks: coordinate order, envelopes and native bins."""
import itertools
from pathlib import Path
import tempfile
import unittest

import numpy as np

import collect_decay_scale_report as collect
from render_decay_scale_report import point_matched_ratio


class DecayScaleReportTests(unittest.TestCase):
    def setUp(self):
        self.points = sorted(itertools.product((.5, 1., 2.), repeat=4))

    def test_scale_groups_keep_the_declared_correlations(self):
        for name, count in [('common3', 3), ('correlated7', 7),
                            ('combined21', 21), ('combined63', 63),
                            ('all81_diagnostic', 81)]:
            self.assertEqual(len(collect.GROUPS[name]), count)
        response = np.array([t-r for r, f, t, a in self.points])
        self.assertEqual(collect.band(response, self.points, 'correlated7'), (0., 0.))
        self.assertEqual(collect.band(response, self.points, 'combined63'), (-1.5, 1.5))
        self.assertIn((.5, 1., .5, .5), collect.GROUPS['correlated7'])
        self.assertNotIn((.5, 1., 1., 1.), collect.GROUPS['correlated7'])

    def test_signed_runtime_labels_preserve_top_antitop_coordinates(self):
        label = 'dyn=3 muR= 0.500 muF= 1.000 d-6= 2.000 d6= 1.000'
        self.assertEqual(collect.scale_tools.scale_point(label), (.5, 1., 1., 2.))

    def test_point_matched_ratio_cancels_common_responses(self):
        denominator = np.array([10+r+f+2*t+3*a for r, f, t, a in self.points])
        numerator = 2*denominator
        ratios = point_matched_ratio(numerator, denominator)
        self.assertEqual(collect.band(ratios, self.points, 'combined63'), (2., 2.))
        self.assertLess(numerator.min()/denominator.max(), 2.)
        self.assertGreater(numerator.max()/denominator.min(), 2.)

    def test_nonpositive_denominator_remains_undefined(self):
        actual = point_matched_ratio(np.array([2., 3., 4.]), np.array([1., 0., -2.]))
        self.assertEqual(actual[0], 2.)
        self.assertTrue(np.isnan(actual[1:]).all())

    def test_replicas_are_averaged_without_inverse_variance_bias(self):
        central = np.array([[1., 4.], [3., 8.]])
        error = np.array([[.1, .2], [1., 2.]])
        result = collect.average_replicas(central, error, central[..., None])
        np.testing.assert_allclose(result['central'], [2., 6.])
        np.testing.assert_allclose(result['scales'][:,0], [2., 6.])
        np.testing.assert_allclose(result['error'], np.sqrt([1.01,4.04])/2)
        np.testing.assert_allclose(result['empirical_error'], [1., 2.])

    def test_paired_jackknife_preserves_shared_replica_fluctuations(self):
        denominator = np.array([[1., 3.], [2., 4.], [3., 5.]])
        np.testing.assert_allclose(collect.paired_ratio_error(
            2*denominator, denominator), 0.)
        self.assertAlmostEqual(collect.paired_ratio_error(
            np.array([1., 3.]), np.array([2., 2.])), .5)

    def test_extreme_production_points_are_only_in_the_full_grid(self):
        values = np.ones(81)
        values[self.points.index((.5, 2., .5, 2.))] = 1e6
        self.assertEqual(collect.band(values, self.points, 'combined63'), (1., 1.))
        self.assertEqual(collect.band(values, self.points, 'all81_diagnostic'), (1., 1e6))

    def test_hwu_reader_keeps_bin_integrals_and_all_coordinates(self):
        columns = ['xmin', 'xmax', 'central value', 'dy']+[
            f'dyn=3 muR={r} muF={f} d-6={a} d6={t}' for r, f, t, a in self.points]
        text = '##& '+' & '.join(columns)+'\n'
        text += '<histogram> 2 "test observable |T@NLO"\n'
        text += ' '.join(map(str, [0., 5., 7., .2]+[7.]*81))+'\n'
        text += ' '.join(map(str, [5., 20., 11., .3]+[11.]*81))+'\n'
        text += '<\\histogram>\n'
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary)/'test.HwU'
            path.write_text(text)
            actual = collect.read_hwu(path)
        self.assertEqual(actual['titles'], ['test observable'])
        np.testing.assert_array_equal(actual['edges'], [[0., 5.], [5., 20.]])
        np.testing.assert_array_equal(actual['central'], [7., 11.])
        self.assertEqual(actual['scales'].shape, (2, 81))


if __name__ == '__main__':
    unittest.main()
