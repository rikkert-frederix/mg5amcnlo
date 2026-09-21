import itertools
from types import SimpleNamespace
import unittest

from audit_native_current import native_vectors


class TestNativeCurrent(unittest.TestCase):
    def histograms(self):
        labels=['central value','dy']
        labels+=['muR= %g muF= %g' % pair for pair in itertools.product((.5,1.,2.),repeat=2)]
        labels+=['PDF= %d' % p for p in range(331700,331801)]
        return [SimpleNamespace(title=title,bins=[
            SimpleNamespace(boundaries=(i+.5,i+1.5),wgts=dict.fromkeys(labels,.1 if i==0 else 0.))
            for i in range(5)]) for title in ('total rate','total rate Born')]

    def test_preserve_all_template_bins_and_born_duplicate(self):
        result=native_vectors(self.histograms())
        self.assertEqual(result['values'].shape,(10,112))
        self.assertEqual(len(result['indices']),183)
        self.assertEqual(result['offsets'],[0,5,10])

    def test_unexpected_support_or_born_difference_rejected(self):
        for h,b,key,value in [(0,1,'central value',.2),(1,0,'central value',.2),
                              (0,0,'dy',-.1),(0,0,'dy',float('nan'))]:
            histograms=self.histograms()
            histograms[h].bins[b].wgts[key]=value
            with self.assertRaises(ValueError):
                native_vectors(histograms)
        with self.assertRaises(ValueError):
            native_vectors(self.histograms()[:1])


if __name__=='__main__':
    unittest.main()
