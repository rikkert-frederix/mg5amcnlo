"""Three-treatment covariance identities and actual pilot-input checks."""
import copy
import json
import unittest

import numpy as np

from campaign import STUDY
from load_results import load
from w_treatment_pilot_report import LABELS, MODES, estimate, validate_physics


def pilot_jobs(charge):
    inputs=STUDY/'inputs'
    jobs=[]
    for name in ('rng_pilot_queue_alignment_controls_v1.json','width_mass_pilot_queue_width_mass_v4.json'):
        jobs.extend(json.loads((inputs/name).read_text())['jobs'])
    selected={}
    for job in jobs:
        report=json.loads((STUDY/'results/audits'/ (job['process'].split('/')[-1]+'_'+job['run']+'.json')).read_text())
        generated=report['execution']['export_generation']['args']
        if generated['charge']==charge and generated['decay_bottom_mass']==0.:
            key=generated['w_treatment'],job['variant']
            if key in selected: raise ValueError('Duplicate physical pilot context')
            selected[key]=job
    return [selected[mode,variant] for mode in MODES for variant in ('S','Pi')]


class TestWTreatmentPilotReport(unittest.TestCase):
    def test_shared_intermediate_mode_covariance_and_additive_closure(self):
        # Two independent five-batch strata, with their integration weights.
        spread=np.asarray([.9,.95,1.,1.05,1.1])
        base=np.r_[spread,2.*spread][:,None,None]/5.
        vectors=[factor*base*np.asarray([1.,2.])[None,:,None]*np.linspace(.8,1.2,183)[None,None,:]
                 for factor in (10.,12.,11.,14.,13.,16.)]
        strata=[np.repeat([0,1],5)]*6
        result=estimate(vectors,strata)
        change=LABELS.index('product_shift_change')
        ratios=LABELS.index('product_ratio_change')
        for index in (change,ratios):
            np.testing.assert_allclose(result['value'][2,index],result['value'][0,index]+result['value'][1,index])
            np.testing.assert_allclose(result['covariance_factor'][:,2,index],
                result['covariance_factor'][:,0,index]+result['covariance_factor'][:,1,index],atol=2.e-13)
        # The shared top-bw S/Pi contribution has opposite signs in these two changes.
        factor=result['covariance_factor'][:,:,change]
        covariance=np.sum(factor[:,0]*factor[:,1],axis=0)
        expected=sum(sum(np.var((row[label==k]*5),axis=0,ddof=1)/5. for k in (0,1))
                     for row,label in zip(vectors[2:4],strata[2:4]))
        np.testing.assert_allclose(covariance,-expected,rtol=1.e-12)
        permuted=list(vectors); permuted[2]=permuted[2][[4,0,2,1,3,9,5,7,6,8]]
        changed=estimate(permuted,strata)
        np.testing.assert_allclose(changed['mc_error'],result['mc_error'],atol=1.e-12)
        normalized=estimate(vectors,strata,lambda row:row[:-1]/row[-1])
        np.testing.assert_allclose(normalized['value'][:,:4],.5)
        self.assertLess(float(np.nanmax(normalized['mc_error'])),1.e-12)

    def test_actual_both_charge_inputs_and_wrong_width_or_scale_controls(self):
        for charge in ('plus','minus'):
            results=[load(job['audit']) for job in pilot_jobs(charge)]
            self.assertEqual(validate_physics(results),dict(charge=charge,flavours=['e','e','mu']))
            for key,value in (('top_width_nlo',1.),('w_width',1.),('production_scale','core-ht-half'),
                              ('top_width_w_treatment','onshell')):
                changed=[dict(row,report=copy.deepcopy(row['report'])) for row in results]
                changed[2]['report']['manifest'][key]=value
                with self.assertRaisesRegex(ValueError,'unmatched benchmark'): validate_physics(changed)
            changed=[dict(row,report=copy.deepcopy(row['report'])) for row in results]
            changed[4]['report']['manifest']['settings']['mur_over_ref']=2.
            with self.assertRaisesRegex(ValueError,'common dynamic'): validate_physics(changed)
            changed=[dict(row,report=copy.deepcopy(row['report'])) for row in results]
            changed[3]['report']['execution']['export_generation']['args']['flavours']=['e','mu','e']
            with self.assertRaisesRegex(ValueError,'Physical inputs'): validate_physics(changed)


if __name__=='__main__':
    unittest.main()
