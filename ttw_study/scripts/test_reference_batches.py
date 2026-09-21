import unittest

from reference_batches import POINTS, columns


class TestReferenceBatches(unittest.TestCase):
    def labels(self, decay=True):
        return ['central value','dy']+[
            'dyn=3 muR= %6.3f muF= %6.3f' % (r,f)+(' d6= %6.3f' % r if decay else '')
            for r,f in POINTS]

    def test_common_scale_labels(self):
        for decay in (False,True):
            self.assertEqual(columns(self.labels(decay)),[0]+list(range(2,11)))

    def test_reject_mismatched_decay_or_duplicate(self):
        labels = self.labels()
        labels[2] = labels[2].replace('d6=  1.000','d6=  2.000')
        with self.assertRaises(ValueError):
            columns(labels)

    def test_native_combined_output_auxiliary_bands_are_not_samples(self):
        labels = self.labels()
        aux = ['delta_mu_'+kind+' -2 @aux' for kind in ('cen', 'min', 'max')]
        self.assertEqual(columns(labels[:2]+aux+labels[2:]), [0]+list(range(5, 14)))
        with self.assertRaises(ValueError):
            columns(labels+['unknown @aux'])
        labels = self.labels()
        labels[3] = labels[2]
        with self.assertRaises(ValueError):
            columns(labels)
