import unittest

import numpy as np

from full_flavour_statistics import FLAVOURS,estimate
from replica_statistics import ratio


class TestFullFlavourStatistics(unittest.TestCase):
    def samples(self):
        base=np.arange(1.,6.)[:,None,None]
        return {(charge,variant,flavour):(i+1)*scale*base for charge in ('plus','minus')
                for variant,scale in (('S',1.),('Pi',1.2)) for i,flavour in enumerate(FLAVOURS)}

    def test_average_retrainings_before_disjoint_sum_and_independent_variance(self):
        samples=self.samples()
        result=estimate(samples,variants=('S','Pi'))
        scale=sum(range(1,9))
        np.testing.assert_allclose(result['value'][:,0,0],[3*scale,3*scale*1.2]*2)
        expected=sum(i*i for i in range(1,9))*np.var(np.arange(1.,6.),ddof=1)/5
        self.assertAlmostEqual(result['mc_error'][0,0,0]**2,expected)
        # Equal-valued simulated inputs are explicitly independent ensembles.
        # Their charge errors add; they cannot cancel through replica-index pairing.
        difference=estimate(samples,variants=('S','Pi'),transform=lambda totals:
            totals['plus','S']-totals['minus','S'])
        self.assertAlmostEqual(difference['mc_error'][0,0]**2,2*expected)

    def test_independent_order_does_not_manufacture_covariance(self):
        samples=self.samples()
        callback=lambda totals:np.asarray([ratio(totals['plus','Pi'],totals['plus','S']),
            ratio(totals['plus','S'],totals['minus','S'])])
        a=estimate(samples,variants=('S','Pi'),transform=callback)
        shuffled={key:(value[::-1] if i%2 else value) for i,(key,value) in enumerate(samples.items())}
        b=estimate(shuffled,variants=('S','Pi'),transform=callback)
        np.testing.assert_allclose(a['value'],b['value'])
        np.testing.assert_allclose(a['mc_error'],b['mc_error'])
        np.testing.assert_allclose(a['value'][:,0,0],[1.2,1.])
        self.assertTrue((a['mc_error']>0).all())

    def test_normalization_uses_complete_flavour_means_with_correlated_parent(self):
        samples={key:np.concatenate([value*.3,value],axis=1) for key,value in self.samples().items()}
        result=estimate(samples,variants=('S','Pi'),transform=lambda totals:
            np.asarray([ratio(row[0],row[1]) for row in totals.values()]))
        np.testing.assert_allclose(result['value'],.3)
        np.testing.assert_allclose(result['mc_error'],0.,atol=1.e-14)

    def test_incomplete_flavours_retraining_and_mixed_shapes_are_rejected(self):
        for kind in ('flavour','replica','shape','nonfinite'):
            samples=self.samples(); key=next(iter(samples))
            if kind=='flavour':del samples[key]
            elif kind=='replica':samples[key]=samples[key][:4]
            elif kind=='shape':samples[key]=np.ones((5,2,1))
            else:samples[key][0]=np.nan
            with self.assertRaises(ValueError):estimate(samples,variants=('S','Pi'))


if __name__=='__main__':
    unittest.main()
