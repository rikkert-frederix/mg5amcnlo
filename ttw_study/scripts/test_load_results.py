import itertools
import copy
from pathlib import Path
import tempfile
import unittest

import numpy as np

from load_results import FACTORS, POINTS, canonical_columns, parameter_signature, strict_identity


class TestCanonicalWeights(unittest.TestCase):
    def columns(self,variant):
        columns = ['central value','dy']
        if variant in ('LO','P'):
            columns += ['dyn=3 muR=%g muF=%g' % point
                        for point in itertools.product(FACTORS,repeat=2)]
        else:
            columns += ['dyn=3 muR=%g muF=%g d-6=%g d6=%g' % (r,f,a,t)
                        for r,f,t,a in POINTS]
        return columns+['PDF=%d NNPDF40_nlo_as_01180' % i for i in range(331700,331801)]

    def test_broadcast_inactive_decay_axes_only(self):
        for variant in ('LO','P','D','S','PiD','Pi'):
            indices = canonical_columns(self.columns(variant),variant)
            self.assertEqual(len(indices),183)
            self.assertEqual(len(set(indices[1:82])),9 if variant in ('LO','P') else 81)
            if variant in ('LO','P'):
                self.assertEqual(indices[1:10],[indices[1]]*9)
            with self.assertRaises(ValueError):
                canonical_columns(self.columns(variant)[0:-1],variant)
        with self.assertRaises(ValueError):
            canonical_columns(self.columns('LO'),'S')
        with self.assertRaises(ValueError):
            canonical_columns(self.columns('S'),'LO')

    def test_incomplete_or_duplicate_grid_is_not_silently_used(self):
        columns = self.columns('S')
        with self.assertRaises(ValueError):
            canonical_columns(columns+[columns[2]],'S')
        with self.assertRaises(ValueError):
            canonical_columns(columns+[columns[-1]],'S')
        with self.assertRaises(ValueError):
            canonical_columns(columns[:2]+columns[3:],'S')


class TestMatchedIdentity(unittest.TestCase):
    def test_parameter_comments_and_order_are_not_physics(self):
        with tempfile.TemporaryDirectory() as temporary:
            a,b,c = [Path(temporary)/name for name in ('a.dat','b.dat','c.dat')]
            a.write_text('BLOCK MASS # initial\n 6 172.5 # top\n 24 80.385\nDECAY 6 1.35\n')
            b.write_text('BLOCK MASS #   changed comment\n 24 80.385\n 6 1.725e2\nDECAY 6 1.35\n')
            c.write_text('BLOCK MASS\n 6 173.5\n 24 80.385\nDECAY 6 1.35\n')
            self.assertEqual(parameter_signature(a),parameter_signature(b))
            self.assertNotEqual(parameter_signature(a),parameter_signature(c))

    def samples(self):
        samples = {}
        for seed, (variant, value) in enumerate((('LO',1.),('P',2.),('D',3.),('S',4.)),1):
            manifest = dict(settings=dict(iseed=seed,lhaid=331700,ebeam1=6500.,
                                           req_acc_fo=.01*seed,fo_job_target_time=60.),
                            top_width_lo=1.48,top_width_nlo=1.35,w_width=2.09767,
                            top_width_reference_scale=172.5,top_width_w_treatment='onshell',
                            top_width_bottom_mass=0.,production_scale='core-w-ht-half',
                            decay_scale_grouping='separate',
                            hashes={'analysis_source':'same analysis','param_card.dat':'same parameters'})
            args = dict(charge='plus',flavours=['e','e','mu'],w_treatment='onshell',
                        decay_bottom_mass=0.,corrected='both')
            samples[variant] = dict(values=np.full((2,183),value),errors=np.full(2,.1*seed),
                                    edges=np.array([[0.,1.],[1.,2.]]),offsets=np.array([0,2]),
                                    titles=['test histogram'],
                                    report=dict(variant=variant,manifest=manifest,
                                                execution={'export_generation':{'args':args}}))
        return samples

    def test_identity_and_independent_error_propagation(self):
        result = strict_identity(self.samples())
        np.testing.assert_array_equal(result['difference'],np.zeros((2,183)))
        np.testing.assert_allclose(result['central_mc_error'],np.sqrt(.3))
        np.testing.assert_array_equal(result['central_pull'],[0.,0.])

    def test_different_physics_is_rejected(self):
        for category, key, value in (
                ('settings','lhaid',14400),('settings','ebeam1',6800.),
                ('hashes','param_card.dat','different weak parameters'),
                ('hashes','analysis_source','different cuts')):
            samples = self.samples()
            samples['S']['report']['manifest'][category][key] = value
            with self.assertRaises(ValueError):
                strict_identity(samples)
        for key, value in (('top_width_nlo',1.4),('top_width_bottom_mass',4.8),
                           ('production_scale','fixed'),
                           ('physical_scale_coordinates',{'numerator_factors':[1,1,2,.5]}),
                           ('custom_user_hook',{'sha256':'different scale hook'})):
            samples = self.samples()
            samples['S']['report']['manifest'][key] = value
            with self.assertRaises(ValueError):
                strict_identity(samples)
        samples = self.samples()
        samples['S']['report']['execution']['export_generation']['args']['flavours'] = ['mu','e','mu']
        with self.assertRaises(ValueError):
            strict_identity(samples)

    def test_corrupted_layout_labels_and_reused_seeds_are_rejected(self):
        original = self.samples()
        for change in ('bin','label','seed'):
            samples = copy.deepcopy(original)
            if change == 'bin':
                samples['S']['edges'][0,1] = 1.5
            elif change == 'label':
                samples['S']['report']['variant'] = 'Pi'
            else:
                samples['S']['report']['manifest']['settings']['iseed'] = 1
            with self.assertRaises(ValueError):
                strict_identity(samples)


if __name__ == '__main__':
    unittest.main()
