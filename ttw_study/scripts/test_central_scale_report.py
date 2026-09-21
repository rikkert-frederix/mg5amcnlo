import copy
import json
from pathlib import Path
import tempfile
import unittest
from unittest import mock

import numpy as np

from campaign import STUDY, digest
from central_scale_report import LABELS, estimate, observables, run, validate_physics
from load_results import load
from replica_statistics import ratio
from run_scale_validation import DONE, cases, physical_scale_coordinates


class TestCentralScaleReport(unittest.TestCase):
    def test_shared_controls_and_independent_deletions(self):
        base=np.arange(1.,7.)[:,None]/6.
        vectors=[scale*base for scale in (1.,2.,3.,5.,4.,7.)]
        strata=[np.zeros(6,dtype=int)]*6
        result=estimate(vectors,strata)
        shuffled=estimate(vectors[:-1]+[vectors[-1][::-1]],strata)
        np.testing.assert_allclose(result['value'],observables([vector.sum(axis=0) for vector in vectors]))
        np.testing.assert_allclose(result['mc_error'],shuffled['mc_error'])
        a=LABELS.index('native__change_in_Pi_minus_S')
        b=LABELS.index('fixed__change_in_Pi_minus_S')
        factor=result['covariance_factor'][:,:,0]
        expected=sum(np.var((vector*6)[:,0],ddof=1)/6 for vector in vectors[:2])
        self.assertAlmostEqual(factor[:,a]@factor[:,b],expected)
        self.assertGreater(expected,0.)

    def test_normalization_uses_correlated_parent_rate(self):
        noise=np.arange(1.,7.)/6.
        vectors=[np.stack([numerator*noise,denominator*noise],axis=1)[:,:,None]
                 for numerator,denominator in ((2,10),(3,10),(4,20),(5,20),(6,30),(7,30))]
        result=estimate(vectors,[np.zeros(6,dtype=int)]*6,lambda row:ratio(row[:-1],row[-1]))
        np.testing.assert_allclose(result['value'][:6,0,0],[.2,.3,.2,.25,.2,7/30])
        np.testing.assert_allclose(result['mc_error'],0.,atol=1.e-14)
        with self.assertRaises(ValueError):
            estimate(vectors[:4],[np.zeros(6,dtype=int)]*4)

    def test_only_the_declared_production_scale_changes(self):
        queue_path=STUDY/'inputs/rng_pilot_queue_alignment_controls_v1.json'
        if not queue_path.exists():
            self.skipTest('Requires the preserved physical all-BW card fixture')
        queue=json.loads(queue_path.read_text())
        base=[load(job['audit']) for job in queue['jobs'] if job['case']['w_treatment']=='all-bw']
        self.assertEqual(len(base),2)
        samples=[]
        for choice in (None,cases()[-2],cases()[-1]):
            for result in base:
                row=dict(result,report=copy.deepcopy(result['report']))
                manifest=row['report']['manifest']
                if choice:
                    manifest.update(production_scale=choice['production_scale'],production_scale_grouping='NONE',
                        technical_test=choice,physical_scale_coordinates=physical_scale_coordinates(choice))
                    manifest['settings']['qes_over_ref']=1.
                    for key in ('fixed_ren_scale','fixed_fac_scale','fixed_qes_scale'):
                        manifest['settings'][key]=choice['production_scale']=='fixed'
                samples.append(row)
        validate_physics(samples)
        for key,value in (('lhaid',14400),('qes_over_ref',.5),('fixed_ren_scale',True)):
            changed=copy.deepcopy(samples)
            changed[2]['report']['manifest']['settings'][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):
                validate_physics(changed)
        changed=copy.deepcopy(samples)
        changed[4]['report']['manifest']['top_width_nlo']=1.
        with self.assertRaises(ValueError):
            validate_physics(changed)

    def test_pending_queue_produces_no_results(self):
        with tempfile.TemporaryDirectory() as directory:
            queue=Path(directory)/'queue.json'
            queue.write_text(json.dumps(dict(status='running technical pilot')))
            output=Path(directory)/'report.json'
            with self.assertRaisesRegex(ValueError,'complete twelve-case'):
                run(queue,output)
            self.assertFalse(output.exists())
            self.assertFalse(output.with_suffix('.npz').exists())

    def test_complete_synthetic_report_keeps_shared_comparisons_and_all_weights(self):
        from narrow_w_report import SPECTRA
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            def job(name,variant,case=None):
                audit=root/(name+'_audit.json')
                manifest={'w_treatment':'all-bw'} if case is None else {'technical_test':case}
                audit.write_text(json.dumps(dict(manifest=manifest)))
                worker=root/(name+'_worker.json')
                worker.write_text(json.dumps(dict(audit=str(audit))))
                batches=root/(name+'_batches.npz')
                batches.write_bytes(b'synthetic; array reader is mocked')
                batches.with_suffix('.json').write_text(json.dumps(dict(
                    input=str(worker),input_sha256=digest(worker))))
                virtual=root/(name+'_virtual.json')
                virtual.write_text(json.dumps(dict(input=str(worker),input_sha256=digest(worker),
                    status='all archived top/antitop analytic virtual checks passed',
                    maximum_relative_difference=1.e-12,tolerance=1.e-8,checked_points=30)))
                return dict(variant=variant,audit=str(audit),batches=str(batches),
                            analytic_virtual_checks=str(virtual))
            controls=root/'controls.json'
            controls.write_text(json.dumps(dict(jobs=[job('grouped_'+v,v) for v in ('S','Pi')])))
            queue=root/'queue.json'
            queue.write_text(json.dumps(dict(status=DONE,cases=cases(),source_hashes={},
                validation_queue=str(controls),reference_queues_sha256={str(controls):digest(controls)},
                jobs=[job(case['name']+'_'+v,v,case) for case in cases() for v in ('S','Pi')])))
            titles=['histogram %d'%i for i in range(1,22)]
            sources=[dict(titles=titles) for _ in range(6)]
            selected=[]
            noise=np.array([.7,.9,1.,1.1,1.3])/5.
            for scale in (1.,1.1,1.2,1.3,1.4,1.5):
                rates=noise[:,None,None]*np.array([100,80,70,50,25,12,8,4,1])[None,:,None]
                rates=rates*scale*np.ones((1,1,183))
                spectra={}
                for local in SPECTRA:
                    parent=rates[:,4 if local in (11,12) else 3,:]
                    spectra[local]=np.stack([.25*parent,.5*parent,parent],axis=1)
                selected.append(dict(strata=np.zeros(5,dtype=int),spectra=spectra,
                    rates={'cut_%d'%i:rates for i in range(5)},
                    edges={local:np.array([[0.,1.],[1.,2.]]) for local in SPECTRA}))
            output=root/'report.json'
            with mock.patch('central_scale_report.load',side_effect=sources), \
                 mock.patch('central_scale_report.validate_physics'), \
                 mock.patch('central_scale_report.verified_pairs',side_effect=[{(i,1)} for i in range(6)]), \
                 mock.patch('central_scale_report.selected_vectors',side_effect=selected):
                run(queue,output)
            record=json.loads(output.read_text())
            self.assertEqual(len(record['rates']),5)
            self.assertEqual(len(record['spectra']),6)
            with np.load(output.with_suffix('.npz')) as arrays:
                self.assertEqual(arrays['rates_0'].shape,(24,9,183))
                np.testing.assert_allclose(arrays['rates_0'][:6,3,0],[50,55,60,65,70,75])
                np.testing.assert_allclose(arrays['normalized_3'][:6,:,0],[[.25,.5]]*6)
                index=LABELS.index('native__double_ratio')
                self.assertAlmostEqual(arrays['rates_0'][index,3,0],(1.3/1.2)/1.1)
                self.assertEqual(arrays['rates_0_covariance_factor_nominal'].shape,(30,24,9))
            self.assertEqual(digest(output.with_suffix('.npz')),record['arrays_sha256'])
            with self.assertRaisesRegex(ValueError,'overwrite'):
                run(queue,output)


if __name__=='__main__':
    unittest.main()
