import unittest

import numpy as np

from run_literature_reference import G0,GN,MU,decay_card,unexpanded


class TestReferenceConventions(unittest.TestCase):
    def test_reference_card_with_actual_compiled_decay_runtime(self):
        import tests.unit_tests.fks.test_fnlo_decay_card as fixture
        case = fixture.TestFNLODecayCard()
        fixture.TestFNLODecayCard.setUpClass()
        try:
            for variant in ('LO','P','S'):
                # The fixture's detailed scale mode assumes active NLO
                # decay axes; LO/P deliberately have none. Use its order
                # and grid modes for those prescriptions.
                orders = case.run_card(decay_card(variant),'orders')
                counterterm = -2*(GN-G0)/G0 if variant == 'S' else 0.
                self.assertAlmostEqual(float(orders['COUNTERTERM'][0]),counterterm)
                grid = case.run_card(decay_card(variant),'grid')
                self.assertEqual(int(grid['GRID_COUNT'][0]),9)
            data = case.run_card(decay_card('S'))
            self.assertEqual(list(map(float,data['SCALES'])),[MU,MU])
            np.testing.assert_allclose(list(map(float,data['WIDTHS'])),[GN,GN],rtol=1.e-13)
        finally:
            fixture.TestFNLODecayCard.doClassCleanups()

    def test_explicit_fixed_widths_and_common_numerator_scales(self):
        card = decay_card('S')
        self.assertIn('EXPLICIT = decay_width_scale_mode(6)',card)
        self.assertIn('CORRELATED = decay_scale_variation_mode',card)
        self.assertIn('NLO = production_order',card)
        self.assertIn('NLO = decay_order',card)
        self.assertNotIn('AUTO = decay_width_scale_mode(6)',card)
        self.assertAlmostEqual(MU,212.6925)
        self.assertIn('LO = decay_order',decay_card('P'))
        self.assertIn('LO = production_order',decay_card('LO'))
        with self.assertRaises(ValueError):
            decay_card('Pi')

    def test_unexpanded_reconstruction_undoes_only_width_counterterm(self):
        born = np.array([2.,3.])
        numerator_nlo = np.array([5.,7.])
        g = (GN-G0)/G0
        strict = numerator_nlo-2*g*born
        np.testing.assert_allclose(unexpanded(strict,born),numerator_nlo/(1+g)**2)


if __name__ == '__main__':
    unittest.main()
