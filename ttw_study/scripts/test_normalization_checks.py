import itertools
from pathlib import Path
from types import SimpleNamespace
import unittest

from campaign import ROOT, lhapdf_environment
from native_inclusive import native_vectors
from normalization_report import physical_settings
from run_normalization_checks import adoption_plan, cases, native_commands, native_settings, recovery_plan


class TestNormalizationChecks(unittest.TestCase):
    def histograms(self):
        labels=['central value','dy']
        labels+=['muR= %g muF= %g'%p for p in itertools.product((.5,1.,2.),repeat=2)]
        labels+=['PDF= %d'%p for p in range(331700,331801)]
        return [SimpleNamespace(title=title,bins=[SimpleNamespace(
            boundaries=(i+.5,i+1.5),wgts=dict.fromkeys(labels,.1 if i==0 else 0.))
            for i in range(5)]) for title in ('total rate','total rate Born')]

    def test_nlo_total_is_not_replaced_by_born_and_signed_weights_survive(self):
        histograms=self.histograms()
        histograms[0].bins[0].wgts['central value']=-.3
        row=native_vectors(histograms,'P')
        self.assertEqual(row['values'][0,0],-.3)
        self.assertEqual(row['values'][5,0],.1)
        self.assertEqual(len(row['indices']),183)
        with self.assertRaisesRegex(ValueError,'LO total and Born'):
            native_vectors(histograms,'LO')
        histograms[0].bins[2].wgts['central value']=1.
        with self.assertRaisesRegex(ValueError,'unused native bins'):
            native_vectors(histograms,'P')

    def test_complete_distinct_cases_and_charge_correct_native_currents(self):
        rows=cases()
        self.assertEqual(len(rows),14)
        self.assertEqual(len({r['seed'] for r in rows}),14)
        expected=set(itertools.product(('stable','current'),('plus','minus'),('LO','P')))
        self.assertEqual({(r['kind'],r['charge'],r['variant']) for r in rows[:8]},expected)
        for charge,current in [('plus','mu+ vm'),('minus','mu- vm~')]:
            command=native_commands(Path('/unused'), 'current',charge)
            self.assertIn('t t~ '+current+' QCD=2 QED=2 [QCD]',command)
            self.assertNotIn('decay',command)

    def test_fixed_current_and_dynamic_stable_scales_and_qes_matching(self):
        reference=dict(w_treatment='onshell',production_scale='core-w-ht-half',settings=dict(lhaid=331700))
        for case in cases()[:8]:
            settings=native_settings(reference,case,.005)
            self.assertEqual(settings['fixed_ren_scale'],case['kind']=='current')
            self.assertEqual(settings['fixed_fac_scale'],case['kind']=='current')
            self.assertEqual(settings['fixed_qes_scale'],case['kind']=='current')
            self.assertEqual(settings['mur_ref_fixed'],212.6925)
            self.assertEqual(settings['fo_job_min_splits'],5)
        a=dict(lhaid=331700,iseed=1,req_acc_fo=.02,fo_job_min_splits=5)
        b=dict(lhaid=331700,iseed=2,req_acc_fo=.001,qes_over_ref=1.)
        self.assertEqual(physical_settings(a),physical_settings(b))
        b['qes_over_ref']=.5
        self.assertNotEqual(physical_settings(a),physical_settings(b))

    def test_stopped_queue_recovery_replaces_only_the_failed_case(self):
        selected=cases()
        previous=dict(status='stopped',cases=selected,current_case=selected[2],
                      jobs=[dict(case=selected[0]),dict(case=selected[1])])
        plan=recovery_plan(previous,82015)
        self.assertEqual(plan['carried_jobs'],previous['jobs'])
        self.assertEqual(plan['cases'][:2],selected[:2])
        self.assertEqual(plan['cases'][2],dict(selected[2],seed=82015))
        self.assertEqual(plan['cases'][3:],selected[3:])
        for altered in (dict(previous,status='running'),
                        dict(previous,jobs=[dict(case=selected[1])]),
                        dict(previous,current_case=dict(selected[2],seed=99999))):
            with self.assertRaises(ValueError):
                recovery_plan(altered,82015)
        with self.assertRaises(ValueError):
            recovery_plan(previous,selected[0]['seed'])

    def test_reader_recovery_can_adopt_only_the_stopped_case(self):
        selected=cases()
        previous=dict(status='stopped',cases=selected,current_case=selected[2],
                      jobs=[dict(case=selected[0]),dict(case=selected[1])])
        adopted=dict(case=selected[2],batches='/immutable/batches.npz')
        plan=adoption_plan(previous,adopted)
        self.assertEqual(plan['cases'],selected)
        self.assertEqual(plan['carried_jobs'],previous['jobs']+[adopted])
        self.assertEqual(plan['adopted_case'],selected[2])
        with self.assertRaises(ValueError):
            adoption_plan(previous,dict(case=selected[3]))

    def test_local_lhapdf_data_path_precedes_inherited_paths(self):
        local=str(ROOT/'LHAPDF/share/LHAPDF')
        environment=lhapdf_environment({'LHAPDF_DATA_PATH':'/old/pdfsets:'+local})
        self.assertEqual(environment['LHAPDF_DATA_PATH'],local+':/old/pdfsets')
        self.assertTrue((Path(local)/'lhapdf.conf').is_file())
        self.assertTrue((Path(local)/'NNPDF40_nlo_as_01180/NNPDF40_nlo_as_01180.info').is_file())


if __name__=='__main__':
    unittest.main()
