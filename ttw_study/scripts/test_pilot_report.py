import unittest

import numpy as np

from load_results import POINTS
from pilot_report import clean, describe, independent_linear
from replica_statistics import ratio


class TestPilotDescriptions(unittest.TestCase):
    def test_four_prescription_contrast_and_independence(self):
        results = {name:dict(values=np.full((1,183),value),errors=np.array([.1*seed]),
                             report={'manifest':{'settings':{'iseed':seed}}})
                   for seed,(name,value) in enumerate((('Pi',5.),('S',3.),('PiD',2.),('D',1.)),1)}
        contrast = independent_linear(results,0,{'Pi':1.,'S':-1.,'PiD':-1.,'D':1.})
        self.assertEqual(contrast['scale']['central'],1.)
        self.assertAlmostEqual(contrast['central_mc_error_pb'],np.sqrt(.3))
        results['D']['report']['manifest']['settings']['iseed'] = 1
        with self.assertRaises(ValueError):
            independent_linear(results,0,{'Pi':1.,'D':-1.})

    def test_expected_scale_slices_and_pdf_separation(self):
        weights = np.ones(183)
        for i,(r,f,t,a) in enumerate(POINTS,1):
            weights[i] = r+2*f+3*t+5*a
        weights[82:] = 1000.
        result = describe(weights)
        for group, count in (('production7',7),('top3',3),('antitop3',3),
                             ('decay3',3),('decay9',9),('combined21',21),
                             ('combined63',63),('common3',3),('all81_diagnostic',81)):
            self.assertEqual(result['scale'][group]['npoints'],count)
            self.assertEqual(result['scale'][group]['nvalid'],count)
        self.assertLess(result['scale']['combined63']['envelope'][1],1000.)
        self.assertEqual(float(result['pdf']['nominal']),1000.)

    def test_ratios_precede_envelopes_and_no_missing_point_is_dropped(self):
        denominator = np.ones(183)
        denominator[1:82] = np.arange(1.,82.)
        numerator = 3*denominator
        result = describe(ratio(numerator,denominator))
        self.assertEqual(result['scale']['combined63']['envelope'],[3.,3.])
        corner = 1+POINTS.index((.5,2.,1.,1.))
        denominator[corner] = 0.
        result = describe(ratio(numerator,denominator))
        self.assertIsNone(result['scale']['all81_diagnostic']['envelope'])
        self.assertEqual(result['scale']['all81_diagnostic']['nvalid'],80)
        self.assertEqual(result['scale']['combined63']['envelope'],[3.,3.])

    def test_json_cleaner_handles_numpy_and_undefined_ratios(self):
        value = clean(dict(x=np.array([1.,np.nan,np.inf]),b=np.bool_(True)))
        self.assertEqual(value,dict(x=[1.,None,None],b=True))


if __name__ == '__main__':
    unittest.main()
