import unittest
import copy
import numpy as np
from mass_effects import contrasts, estimate, validate_physics
from replica_statistics import ratio


class TestMassEffects(unittest.TestCase):
    def test_actual_four_sample_physics_gate_and_negative_controls(self):
        from campaign import STUDY
        from load_results import load
        names=('TTWplus_onshell_eemu_mb0p0_both_rng_v1_S_onshell_core-w-ht-half_separate_57003',
               'TTWplus_onshell_eemu_mb0p0_both_rng_v1_Pi_onshell_core-w-ht-half_separate_57004',
               'TTWplus_onshell_eemu_mb4p8_both_width_mass_v3_S_onshell_core-w-ht-half_separate_60001_mb4p8',
               'TTWplus_onshell_eemu_mb4p8_both_alignment_v1_Pi_onshell_core-w-ht-half_separate_61001_mb4p8')
        paths=[STUDY/'results/audits'/(name+'.json') for name in names]
        if not all(p.exists() for p in paths):
            self.skipTest('The four archived study pilots are not available')
        identities=[STUDY/'inputs'/name for name in
                    ('width_mass_v3_plus_onshell_mb4p8_production_identity.json',
                     'alignment_v1_production_identity.json')]
        results=[load(p) for p in paths]
        self.assertEqual(validate_physics(results,identities),4.8)
        # A separate width table cannot silently replace archived provenance.
        custom=STUDY/'inputs/small_mass_inputs_v1.json'
        if custom.exists():
            with self.assertRaisesRegex(ValueError,'not recorded'):
                validate_physics(results,identities,width_inputs=custom)
            from campaign import digest
            declared=[dict(r,report=copy.deepcopy(r['report'])) for r in results]
            for row in declared:
                row['report']['manifest']['small_mass_inputs_sha256']=digest(custom)
            self.assertEqual(validate_physics(declared,identities,width_inputs=custom),4.8)
        with self.assertRaises(ValueError):
            validate_physics(results,[])
        for field,value in [('top_width_nlo',1.),('w_width',1.),('decay_bottom_mass',0.)]:
            changed=[dict(r,report=copy.deepcopy(r['report'])) for r in results]
            for row in changed[2:]:
                row['report']['manifest'][field]=value
            with self.assertRaises(ValueError):
                validate_physics(changed,identities)
        changed=[dict(r,report=copy.deepcopy(r['report'])) for r in results]
        for row in changed[2:]:
            row['report']['manifest']['settings']['mur_over_ref']=2.
        with self.assertRaises(ValueError):
            validate_physics(changed,identities)

    def test_definitions_before_weight_reduction(self):
        np.testing.assert_allclose(contrasts(10.,12.,11.,15.),
            [1.1,1.25,2.,4.,2.,1.2,15./11.,15./11.-1.2,(15./11.)/1.2])
        self.assertTrue(np.isnan(contrasts(0.,1.,2.,3.)[0]))

    def test_four_independent_predictions_not_arbitrarily_paired(self):
        base=np.arange(1.,7.)[:,None]/6
        vectors=[scale*base for scale in (1.,2.,3.,5.)]
        strata=[np.zeros(6,dtype=int)]*4
        result=estimate(vectors,strata)
        shuffled=estimate(vectors[:3]+[vectors[3][::-1]],strata)
        np.testing.assert_allclose(result['value'],shuffled['value'])
        np.testing.assert_allclose(result['mc_error'],shuffled['mc_error'])
        expected=np.sqrt(sum(np.var((v*6)[:,0],ddof=1)/6 for v in vectors))
        self.assertAlmostEqual(result['mc_error'][4,0],expected)

    def test_bin_parent_rate_correlations_are_retained(self):
        base=np.arange(1.,7.)/6
        vectors=[np.stack([scale*base,denominator*base],axis=1)
                 for scale,denominator in ((1.,2.),(2.,3.),(3.,4.),(4.,5.))]
        result=estimate(vectors,[np.zeros(6,dtype=int)]*4,
            lambda *rows:contrasts(*[ratio(row[0],row[1]) for row in rows]))
        np.testing.assert_allclose(result['mc_error'],0.,atol=1.e-14)


if __name__=='__main__':
    unittest.main()
