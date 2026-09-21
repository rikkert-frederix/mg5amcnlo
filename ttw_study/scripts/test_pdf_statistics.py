import math
import unittest

import numpy as np

from pdf_statistics import uncertainty
from replica_statistics import ratio


class TestPDFStatistics(unittest.TestCase):
    def test_replica_standard_deviation_not_error_of_mean(self):
        members = np.r_[1000.,np.arange(100.)]
        result = uncertainty(members)
        self.assertEqual(float(result['nominal']),1000.)
        self.assertAlmostEqual(float(result['library_central']),49.5)
        self.assertAlmostEqual(float(result['error_symmetric']),np.std(members[1:],ddof=1),places=5)
        self.assertAlmostEqual(float(result['library_central_minus_nominal']),-950.5)

    def test_member_wise_ratio_preserves_cancellation(self):
        numerator = np.linspace(1.,2.,101)
        result = uncertainty(ratio(numerator,2*numerator))
        self.assertEqual(float(result['nominal']),.5)
        self.assertEqual(float(result['error_symmetric']),0.)

    def test_nonlinear_member_mean_is_not_substituted_for_nominal(self):
        denominator = np.linspace(1.,2.,101)
        result = uncertainty(ratio(np.ones(101),denominator))
        self.assertEqual(float(result['nominal']),1.)
        self.assertAlmostEqual(float(result['library_central']),np.mean(1/denominator[1:]))
        self.assertNotAlmostEqual(float(result['library_central']),1/np.mean(denominator[1:]))

    def test_undefined_ratio_does_not_drop_pdf_members(self):
        members = np.ones((2,101))
        members[1,20] = np.nan
        result = uncertainty(members)
        np.testing.assert_array_equal(result['valid'],[True,False])
        self.assertTrue(np.isnan(result['error_symmetric'][1]))
        with self.assertRaises(ValueError):
            uncertainty(np.ones(100))

    def test_ct18_90_percent_errors_are_rescaled(self):
        import lhapdf
        members = np.r_[10.,np.tile([11.,9.],29)]
        result = uncertainty(members,'CT18NLO')
        self.assertEqual(result['error_type'],'hessian')
        self.assertAlmostEqual(result['confidence_level_percent'],68.2689492137,places=9)
        original = lhapdf.getPDFSet('CT18NLO').uncertainty(members.tolist(),90.,False)
        self.assertGreater(original.errsymm,float(result['error_symmetric']))
        self.assertAlmostEqual(original.errsymm/float(result['error_symmetric']),1.64485362695,places=5)


if __name__ == '__main__':
    unittest.main()
