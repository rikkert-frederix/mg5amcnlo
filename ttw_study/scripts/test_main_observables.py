import unittest

import numpy as np

from full_flavour_statistics import FLAVOURS,VARIANTS,estimate
from main_observables import (ABSOLUTE_LABELS,LABELS,charge_observables,
    normalize_spectrum,prescription_observables)


class TestMainObservables(unittest.TestCase):
    def totals(self):
        # S=P+D-LO in each bin, but not for the separately normalized shapes.
        values=dict(LO=[10.,20.],P=[12.,25.],D=[9.,19.],S=[11.,24.],Pi=[11.5,24.8],PiD=[9.2,19.2])
        return {(charge,variant):np.asarray(row)[:,None]*scale
                for charge,scale in (('plus',1.),('minus',.6)) for variant,row in values.items()}

    def test_linear_identity_is_not_imposed_after_normalization(self):
        totals=self.totals()
        absolute=prescription_observables(totals,'plus')
        normalized=prescription_observables(totals,'plus',normalize_spectrum)
        self.assertEqual(len(absolute),len(ABSOLUTE_LABELS))
        self.assertEqual(len(normalized),len(LABELS))
        np.testing.assert_allclose(absolute[-1],0.)
        shape={v:normalize_spectrum(totals['plus',v]) for v in VARIANTS}
        self.assertGreater(abs((shape['S']-shape['P']-shape['D']+shape['LO'])[0,0]),.001)
        np.testing.assert_allclose(normalized[LABELS.index('Pi_over_S')],shape['Pi']/shape['S'])

    def test_charge_ratio_and_asymmetry_have_the_correct_denominators(self):
        total=self.totals()
        result=charge_observables(total)
        np.testing.assert_allclose(result[3],1/.6)
        np.testing.assert_allclose(result[4],.4/1.6)
        normalized=charge_observables(total,normalize_spectrum)
        # The combined-charge shape is normalized AFTER adding cross sections.
        np.testing.assert_allclose(normalized[2],normalized[0])
        np.testing.assert_allclose(normalized[3],1.)
        np.testing.assert_allclose(normalized[4],0.,atol=1.e-14)
        zero={key:np.zeros_like(row) for key,row in total.items()}
        self.assertTrue(np.isnan(charge_observables(zero)[3:]).all())

    def test_complete_flavour_covariance_propagates_to_strict_identity(self):
        total=self.totals()
        spread=np.arange(1.,6.)/3
        samples={(charge,variant,flavour):spread[:,None,None]*total[charge,variant][None,...]/8
                 for charge in ('plus','minus') for variant in VARIANTS for flavour in FLAVOURS}
        result=estimate(samples,transform=lambda sums:prescription_observables(sums,'plus'))
        np.testing.assert_allclose(result['value'][-1],0.,atol=1.e-13)
        self.assertTrue((result['mc_error'][-1]>0.).all())
        expected=sum(np.var(spread,ddof=1)/5*total['plus',v]**2/8 for v in ('S','P','D','LO'))
        np.testing.assert_allclose(result['mc_error'][-1]**2,expected)


if __name__=='__main__':
    unittest.main()
