import copy
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from pilot_report import RATES
from run_scale_validation import cases, overlapping_points
from technical_inspection import frozen_source_drift, residual_summary, run, summarize_identity


class TestTechnicalInspection(unittest.TestCase):
    def fixture(self,case):
        points=[row[2] for row in overlapping_points(case['factors'])]
        values=np.ones((4,9,len(points)))
        values[0]=2.
        errors=np.ones_like(values)
        values[:,7:]=0.
        values[3,7:]=np.nan
        errors[:,7:]=0.
        errors[3,7:]=np.nan
        record=dict(case=case,reference_points=points,rates={},spectra={})
        arrays=dict(rate_values=values,rate_mc_errors=errors)
        for i,label in enumerate(RATES):
            record['rates'][label]=dict(direct_pb=values[0,i].tolist(),reference_pb=values[1,i].tolist(),
                difference_pb=values[2,i].tolist(),difference_mc_error_pb=errors[2,i].tolist())
        for local in (3,4,6,11,12,20):
            value=np.full((4,2,len(points)),.3)
            value[2]=.02
            error=np.full_like(value,.01)
            arrays['normalized_%d'%local]=value
            arrays['normalized_mc_errors_%d'%local]=error
            record['spectra']['histogram_%d'%local]=dict(finite_nominal_pulls=2,
                maximum_absolute_nominal_pull=2.,nominal_bins_above_three=0,normalized_to='fiducial_1b')
        return record,arrays

    def test_coordinates_nominal_duplicate_and_structural_zeros(self):
        for case,entries,distinct in ((cases()[0],37,36),(cases()[1],82,81)):
            record,arrays=self.fixture(case)
            result=summarize_identity(record,arrays)
            self.assertEqual(result['weight_comparison_entries'],entries)
            self.assertEqual(result['distinct_physical_scale_coordinates'],distinct)
            self.assertEqual(result['all_rate_weights']['finite_comparisons'],7*entries)
            self.assertEqual(result['all_rate_weights']['maximum_absolute_residual'],1.)
            self.assertIsNone(result['nominal_rates']['2b_2extra']['residual_over_mc_error'])
            self.assertEqual(result['normalized_spectra']['histogram_3']['nominal']['maximum_absolute_residual'],2.)
        self.assertEqual(residual_summary([0.,.4,-.9],[0.,.1,.1])['comparisons_above_three'],2)

    def test_mismatched_numerical_evidence_is_rejected(self):
        record,arrays=self.fixture(cases()[0])
        changed=copy.deepcopy(record)
        changed['reference_points'][0]=[1.]*4
        with self.assertRaisesRegex(ValueError,'coordinates'):
            summarize_identity(changed,arrays)
        changed=copy.deepcopy(record)
        changed['rates']['fiducial_1b']['difference_pb'][0]=.5
        with self.assertRaisesRegex(ValueError,'Rate summary'):
            summarize_identity(changed,arrays)
        changed=copy.deepcopy(record)
        changed['spectra']['histogram_3']['maximum_absolute_nominal_pull']=1.
        with self.assertRaisesRegex(ValueError,'Spectrum summary'):
            summarize_identity(changed,arrays)

    def test_incomplete_queue_rejected_before_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            queue=Path(directory)/'pending.json'
            queue.write_text(json.dumps(dict(cases=cases(),status='running technical pilot')))
            output=Path(directory)/'inspection.json'
            with self.assertRaisesRegex(ValueError,'complete technical queue'):
                run(queue,output)
            self.assertFalse(output.exists())
            self.assertFalse(output.with_suffix('.queue.json').exists())

    def test_source_drift_requires_explicit_historical_mode(self):
        source='ttw_study/scripts/technical_inspection.py'
        with self.assertRaisesRegex(ValueError, 'Frozen technical source'):
            frozen_source_drift({source: 'not-the-current-hash'})
        drift=frozen_source_drift({source: 'not-the-current-hash'}, allow_source_drift=True)
        self.assertEqual(drift[source]['frozen_sha256'], 'not-the-current-hash')
        self.assertNotEqual(drift[source]['current_sha256'], 'not-the-current-hash')


if __name__=='__main__':
    unittest.main()
