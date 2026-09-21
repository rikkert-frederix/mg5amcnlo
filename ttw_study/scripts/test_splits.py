import unittest
from pathlib import Path

import numpy as np

from read_splits import split_ensembles, validate_jobs, validate_refinement_history, validate_exclusion_evidence
from load_results import canonical_columns, FACTORS
from replica_statistics import independent_jackknife_many, ratio


class SplitChecks(unittest.TestCase):
    def test_repeated_native_refinements_cannot_certify_independence(self):
        for step in (-1,0,1):
            validate_refinement_history('INFO: Refining results, step %d\n' % step)
        for log in ('INFO: no refinement',
                    'Refining results, step 1\nRefining results, step 2',
                    'Refining results, step 3'):
            with self.assertRaises(ValueError):
                validate_refinement_history(log)

    def test_raw_fortran_pdf_header_padding(self):
        columns = ['central value','dy']
        columns += ['muR= %g muF= %g' % (r,f) for r in FACTORS for f in FACTORS]
        columns += ['PDF= %6d' % member for member in range(331700,331801)]
        self.assertEqual(len(canonical_columns(columns,'LO')),183)
        columns[-1] = 'PDF= invalid'
        with self.assertRaises(ValueError):
            canonical_columns(columns,'LO')

    def workers(self):
        return [dict(job=dict(p_dir=group,channel='1',mint_mode=-1,niters=1,
                             niters_done=1,npoints=1000,npoints_done=1000,
                             split=i,configs='1.1 1.2',nchans=2,accuracy=.1,
                             run_mode='all',wgt_frac=1.,wgt_mult=.2))
                for group in ('first','second') for i in range(1,6)]

    def test_equal_counts_and_weights(self):
        self.assertEqual([len(group) for group in validate_jobs(self.workers()).values()],[5,5])
        for key,value in [('split',2),('wgt_mult',1.),('npoints_done',500),('niters_done',2)]:
            workers = self.workers()
            workers[0]['job'][key] = value
            with self.assertRaises(ValueError):
                validate_jobs(workers)

    def test_strata_sum_and_conditional_covariance(self):
        first = np.arange(1.,6.)
        second = np.arange(20.,25.)
        rows = np.r_[first,second]
        values = np.stack([rows,2*rows],axis=1)/5.
        ensembles = split_ensembles(values,np.repeat([0,1],5))
        result = independent_jackknife_many(ensembles,lambda *means: sum(means))
        np.testing.assert_allclose(result['value'],values.sum(axis=0))
        expected = np.sqrt(np.var(first,ddof=1)/5+np.var(second,ddof=1)/5)
        np.testing.assert_allclose(result['mc_error'],[expected,2*expected])
        fraction = independent_jackknife_many(ensembles,lambda *means: ratio(sum(means)[0],sum(means)[1]))
        self.assertAlmostEqual(float(fraction['value']),.5)
        self.assertLess(float(fraction['mc_error']),1.e-14)

    def test_excluded_split_requires_matching_metadata_and_normalization(self):
        workers=self.workers()[:5]
        for w in workers:
            w['job'].update(wgt_mult=1./6.,split_result_scale=6./5.,
                split_outlier=dict(original_count=6,excluded_splits=[6],audit='/audit.json'))
        self.assertEqual(len(next(iter(validate_jobs(workers).values()))),5)
        for key,value in [('split_result_scale',1.),('wgt_mult',.2),('split',6)]:
            old=workers[0]['job'][key]
            workers[0]['job'][key]=value
            with self.assertRaises(ValueError):validate_jobs(workers)
            workers[0]['job'][key]=old

    def test_exclusion_requires_archived_evidence_for_the_same_raw_result(self):
        job=self.workers()[0]['job']
        job.update(split_result_scale=6./5.,result=1.2,resultABS=1.2,error=.12,errorABS=.12,
            split_outlier=dict(original_count=6,excluded_splits=[6],audit='/process/Events/run/split_outliers/step_1/audit.json'))
        group=dict(p_dir=job['p_dir'],channel=job['channel'],run_mode=job['run_mode'],
                   original_count=6,excluded_split=6,retained_scale=1.2,
                   raw_jobs=[dict(split=1,result=1.,resultABS=1.,error=.1,errorABS=.1)])
        evidence={'Events/run/split_outliers/step_1/audit.json':dict(groups=[group])}
        validate_exclusion_evidence([dict(job=job)],evidence,Path('/process'))
        with self.assertRaisesRegex(ValueError,'Missing archived'):
            validate_exclusion_evidence([dict(job=job)],{},Path('/process'))
        group['raw_jobs'][0]['result']=2.
        with self.assertRaisesRegex(ValueError,'Retained result differs'):
            validate_exclusion_evidence([dict(job=job)],evidence,Path('/process'))


if __name__ == '__main__':
    unittest.main()
