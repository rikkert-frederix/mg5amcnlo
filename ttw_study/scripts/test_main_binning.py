import copy
import unittest

import numpy as np

from main_binning import apply_plan, layout_signature


class TestMainBinning(unittest.TestCase):
    def fixture(self):
        titles,edges,offsets=[],[],[0]
        for bank in range(10):
            for local in range(1,22):
                size=9 if local==1 else 4
                titles.append('bank_%d'%bank+(' rates: stages' if local==1 else ' observable_%d'%local))
                edges.extend([[i,i+1] for i in range(size)])
                offsets.append(offsets[-1]+size)
        layout=dict(titles=titles,edges=np.array(edges),offsets=np.array(offsets),weights=['nominal','varied'])
        values=np.full((5,len(edges),2),100.)
        for bank in range(10):
            first=offsets[bank*21+2]
            noise=np.arange(-2.,3.)
            values[:,first:first+4,0]=np.stack([4+noise,10-noise,20+noise,30-noise],axis=1)
            values[:,first:first+4,1]=2*values[:,first:first+4,0]
        vectors={'sample':values}
        conditional={'sample':np.full((5,len(edges)),3.)}
        plan=dict(schema='main_common_binning_v1',layout_sha256=layout_signature(layout),histograms={'3':[0,2,4]})
        return layout,vectors,conditional,plan

    def test_signed_cancellation_is_recomputed_before_mc_errors(self):
        layout,vectors,conditional,plan=self.fixture()
        old=copy.deepcopy(vectors)
        result,values,diagonal,record=apply_plan(layout,vectors,conditional,plan)
        self.assertEqual(record['original_bins']-record['retained_bins'],20)
        for bank in range(10):
            a,b=result['offsets'][bank*21+2:bank*21+4]
            # Select only the two merged bins of continuous observable 3.
            np.testing.assert_allclose(values['sample'][:,a:a+2,0],[[14,50]]*5)
            np.testing.assert_allclose(values['sample'][:,a:a+2,1],[[28,100]]*5)
            np.testing.assert_allclose(np.var(values['sample'][:,a:a+2],axis=0,ddof=1),0.)
            self.assertTrue(np.isnan(diagonal['sample'][:,a:a+2]).all())
            start=result['offsets'][bank*21]
            np.testing.assert_allclose(diagonal['sample'][:,start:start+9],3.)
        np.testing.assert_array_equal(vectors['sample'],old['sample'])
        self.assertGreater(np.var(old['sample'][:,13,0],ddof=1),0.)
        np.testing.assert_allclose(values['sample'].sum(axis=1),old['sample'].sum(axis=1))

    def test_rate_merging_partial_ranges_and_split_bins_rejected(self):
        layout,vectors,conditional,plan=self.fixture()
        for definitions in ({'1':[0,9]},{'2':[0,2,4]},{'14':[0,2,4]},
                            {'3':[0,2,3]},{'3':[0,.5,4]}, {'3':[0,1,2,3,4]},{}):
            changed=dict(plan,histograms=definitions)
            with self.subTest(definitions=definitions),self.assertRaises(ValueError):
                apply_plan(layout,vectors,conditional,changed)
        with self.assertRaisesRegex(ValueError,'original main histogram layout'):
            apply_plan(layout,vectors,conditional,dict(plan,layout_sha256='wrong'))

    def test_mismatched_cut_binnings_and_diagnostics_rejected(self):
        layout,vectors,conditional,plan=self.fixture()
        changed=copy.deepcopy(layout)
        a,b=changed['offsets'][23:25]
        changed['edges'][a:a+4]*=2
        updated=dict(plan,layout_sha256=layout_signature(changed))
        with self.assertRaisesRegex(ValueError,'different original bins'):
            apply_plan(changed,vectors,conditional,updated)
        with self.assertRaisesRegex(ValueError,'groups differ'):
            apply_plan(layout,vectors,{},plan)


if __name__=='__main__':
    unittest.main()
