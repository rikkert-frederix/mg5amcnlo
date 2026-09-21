import copy
import unittest

from run_sampler_pilots import (cases, cancelled_control, drained_control, run_name, run_settings,
                               sampler_source_guard)


class TestSamplerPilots(unittest.TestCase):
    def test_successor_requires_exact_sampler_source_guard_and_completed_tail(self):
        case=dict(w_treatment='all-bw',variant='Pi',production_sampling='w-current',seed=53004)
        job=dict(process='old',run='Pi',case=case,workers='workers',audit='audit')
        previous=dict(status='stopped',jobs=[job],current_job=dict(process='old',run='Pi',case=case),
                      error=repr(ValueError('Frozen sampler source/input changed; use a fresh inspected queue')))
        self.assertEqual(sampler_source_guard(previous),job)
        self.assertIsNone(sampler_source_guard(dict(status='running')))
        for changes in (dict(error='integration failed'),dict(jobs=[]),dict(jobs=[job,job]),
                        dict(current_job=dict(process='old',run='unfinished',case=case))):
            with self.assertRaises(ValueError):
                sampler_source_guard(dict(previous,**changes))

    def test_distinct_seeds_and_same_source_linear_control(self):
        rows=cases()
        self.assertEqual(len(rows),6)
        self.assertEqual(len({r['seed'] for r in rows}),6)
        self.assertEqual(len({run_name(r) for r in rows}),6)
        self.assertEqual([r['production_sampling'] for r in rows if r['variant']=='LO'],
                         ['w-current','flat'])
        self.assertTrue(all(r['w_treatment']=='all-bw' for r in rows
                            if r['production_sampling']=='w-current'))
        for mode in ('onshell','all-bw'):
            self.assertEqual({r['variant'] for r in rows if r['w_treatment']==mode
                              and r['variant']!='LO'},{'S','Pi'})
        for case in rows:
            settings=run_settings('test',case,.03)
            self.assertIsNone(settings.grid_reference)
            if case['variant']=='LO' and case['production_sampling']=='flat':
                self.assertEqual((settings.accuracy,settings.points,settings.iterations),(-1.,262144,1))
            else:
                self.assertEqual(settings.accuracy,.03)

    def test_only_accept_completed_archived_control_and_intentional_stop(self):
        self.assertIsNone(drained_control(dict(status='running')))
        job=dict(process='old',run='S',variant='S',audit='audit',workers='workers',
                 case=dict(charge='plus',w_treatment='all-bw',decay_bottom_mass=0.))
        previous=dict(status='stopped',error=repr(RuntimeError(
            'Source changed after queue preparation; inspect and start a fresh queue')),
                      current_job=dict(process='old',run='S'),jobs=[job])
        self.assertEqual(drained_control(previous),job)
        for changes in (dict(error='failed MC'),dict(jobs=[]),dict(jobs=[job,job]),
                        dict(current_job=dict(process='old',run='Pi'))):
            with self.assertRaises(ValueError):
                drained_control(dict(previous,**changes))
        invalid=copy.deepcopy(previous)
        invalid['jobs'][0]['case']['w_treatment']='onshell'
        with self.assertRaises(ValueError):
            drained_control(invalid)

    def test_failed_control_never_becomes_a_completed_physics_sample(self):
        case=dict(charge='plus',w_treatment='all-bw',decay_bottom_mass=0.)
        previous=dict(status='stopped',jobs=[],current_job=dict(process='old',run='S',case=case),
                      error="RuntimeError('Integration failed or lacks HwU output; see log')")
        diagnostic=dict(status='frozen nonconverged adaptive diagnostic; excluded from physics',
                        process='old',run='S',archive='snapshot.tar.gz')
        self.assertTrue(cancelled_control(previous,diagnostic)['excluded_from_physics'])
        self.assertIsNone(cancelled_control(dict(status='running'),diagnostic))
        for changes in (dict(process='unrelated'),dict(run='Pi'),dict(status='finished')):
            with self.assertRaises(ValueError):
                cancelled_control(previous,dict(diagnostic,**changes))
        with self.assertRaises(ValueError):
            cancelled_control(dict(previous,jobs=[dict(process='old',run='S')]),diagnostic)


if __name__=='__main__':
    unittest.main()
