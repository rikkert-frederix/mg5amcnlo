import unittest

import numpy as np

from run_ttbar_reference import check_histogram, decay_card, expanded_normalized, normalized, run_settings


class TestTTbarReference(unittest.TestCase):
    def physics(self):
        return dict(mt=172.5,MW=80.385,W_width=2.0928,top_width_LO=1.48063,top_width_NLO=1.3535)

    def test_normalized_expansion_and_bin_integrals(self):
        born = np.ones((22,10)); born[0]=10.; born[1]=10.
        strict = born*1.3
        np.testing.assert_allclose(normalized(strict),np.ones((10,10)))
        np.testing.assert_allclose(expanded_normalized(born,strict),np.ones((10,10)))
        strict[2]*=2
        strict[0]+=1.3
        expected = (strict[2:12]/born[0]-born[2:12]*(strict[0]-born[0])/born[0]**2)/.1
        np.testing.assert_allclose(expanded_normalized(born,strict),expected)
        check_histogram(dict(offsets=np.array([0,2,12,22]),values=born))
        born[2,1]*=2
        with self.assertRaises(ValueError):
            check_histogram(dict(offsets=np.array([0,2,12,22]),values=born))

    def test_fast_channels_receive_enough_covariance_batches(self):
        from madgraph.various.banner import RunCardNLO
        from madgraph.interface.amcatnlo_run_interface import aMCatNLOCmd
        interface=object.__new__(aMCatNLOCmd)
        interface.stop_for_runweb=True
        interface.options=dict(run_mode=2,nb_core=64)
        interface.run_card=RunCardNLO()
        job=dict(p_dir='P0_uxu',channel='1',dirname='/tmp/P0_uxu/all_G1',resultABS=1.,
                 time_spend=.04,niters=2,npoints=26880,niters_done=4,npoints_done=200,
                 combined=1,accuracy=.01)
        counts=[]
        for target in (1.,.05):
            settings=run_settings(dict(settings={}),self.physics(),49501,.005,target)
            for name,value in settings.items():interface.run_card.set(name,value,user=True,raiseerror=True)
            jobs,collected=interface.split_jobs_fixed_order([job],[job])
            counts.append(len(jobs))
            self.assertEqual(len(jobs),len(collected))
            self.assertAlmostEqual(sum(x['wgt_mult'] for x in jobs),1.)
            self.assertTrue(all(x['npoints']>=1000 and x['niters']==1 for x in jobs))
        self.assertLess(counts[0],5)
        self.assertGreaterEqual(counts[1],5)
        self.assertLessEqual(counts[1],64)
        for target in (0.,-1.,float('nan'),float('inf')):
            with self.assertRaises(ValueError):run_settings(dict(settings={}),self.physics(),49501,.005,target)

    def test_runtime_card_fixed_width_and_scale_convention(self):
        import tests.unit_tests.fks.test_fnlo_decay_card as fixture
        case = fixture.TestFNLODecayCard()
        fixture.TestFNLODecayCard.setUpClass()
        try:
            physics = self.physics()
            for variant in ('LO','S','Pi'):
                card = decay_card(variant,physics)
                grid = case.run_card(card,'grid')
                self.assertEqual(int(grid['GRID_COUNT'][0]),9)
                orders = case.run_card(card,'orders')
                # This accessor returns the expansion coefficient even in
                # product mode; the multiplicative integrator does not add
                # that coefficient. Test the mode separately from the value.
                counterterm = -2*(physics['top_width_NLO']/physics['top_width_LO']-1.) if variant != 'LO' else 0.
                self.assertAlmostEqual(float(orders['COUNTERTERM'][0]),counterterm)
                self.assertEqual(orders['ORDERS'][-1],'T' if variant == 'Pi' else 'F')
            data = case.run_card(decay_card('S',physics))
            np.testing.assert_allclose(list(map(float,data['SCALES'])),[172.5,172.5])
            np.testing.assert_allclose(list(map(float,data['WIDTHS'])),[physics['top_width_NLO']]*2)
            product = case.run_card(decay_card('Pi',physics))
            self.assertAlmostEqual(float(product['PRODUCT_WIDTH'][0]),
                                   (physics['top_width_LO']/physics['top_width_NLO'])**2)
        finally:
            fixture.TestFNLODecayCard.doClassCleanups()
