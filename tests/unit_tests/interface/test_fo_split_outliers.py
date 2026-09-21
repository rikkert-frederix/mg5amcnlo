"""Opt-in exclusion of an isolated extreme fixed-order split replica."""
import copy
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest import mock

from madgraph.interface import amcatnlo_run_interface as run
from madgraph.various import banner
from madgraph.madevent import sum_html


class TestSplitOutliers(unittest.TestCase):
    def interface(self, root):
        interface = object.__new__(run.aMCatNLOCmd)
        interface.stop_for_runweb = True
        interface.me_dir = str(root)
        interface.run_name = 'test'
        interface.run_card = banner.RunCardNLO()
        interface.run_card['fo_split_outlier_threshold'] = 10.
        interface.options = {'lhapdf': 'lhapdf-config'}
        interface.get_randinit_seed = lambda: 1
        return interface

    def jobs(self, root, count=8, bad=1, sign=1.):
        jobs = []
        for split in range(1, count+1):
            directory = Path(root)/'SubProcesses'/'P0_test'/('all_G1_%d' % split)
            value, error = (sign*100., 90.) if split == bad else (1., .1)
            jobs.append(dict(dirname=str(directory),p_dir='P0_test',p_label='P0',
                channel='1',split=split,mint_mode=-1,configs='1.1',nchans=1,
                run_mode='all',accuracy=.1,wgt_frac=1.,wgt_mult=1./count,
                niters=1,niters_done=1,npoints=1000,npoints_done=1000,time_spend=2.,
                result=value,resultABS=abs(value),error=error,errorABS=error,
                err_perc=100.*error/value,err_percABS=100.*error/abs(value),
                contribution_results=dict(components=[dict(id=1,label='BORN',value=value,error=error)],
                    sum=value,total=value,error=error,closure=0.)))
        return jobs

    def test_disabled_by_default_and_invalid_configuration(self):
        with tempfile.TemporaryDirectory() as root:
            interface = self.interface(root)
            interface.run_card = banner.RunCardNLO()
            jobs = self.jobs(root)
            before = copy.deepcopy(jobs)
            self.assertEqual(interface.exclude_fixed_order_split_outliers(jobs,jobs[:],1),(jobs,jobs))
            self.assertEqual(jobs,before)
            for key,value in [('fo_split_outlier_threshold',float('nan')),
                              ('fo_split_outlier_threshold',-1.),
                              ('fo_split_outlier_variance_fraction',.5),
                              ('fo_split_outlier_variance_fraction',1.),
                              ('fo_split_outlier_min_splits',5)]:
                interface.run_card = {key:value}
                with self.assertRaises(run.aMCatNLOError):
                    interface.fixed_order_split_outlier_settings()

    def test_detects_both_signs_but_not_broad_or_incomparable_samples(self):
        find = run.aMCatNLOCmd.find_fixed_order_split_outliers
        for sign in (-1.,1.):
            jobs = self.jobs('/unused',sign=sign)
            self.assertEqual(find(jobs,10.)[0]['excluded']['split'],1)
        for change in ('normal','two_extremes','incomplete','unequal_points','unequal_accuracy',
                       'unsplit','different_channels','zero_peer_information','already_filtered'):
            jobs = self.jobs('/unused')
            if change=='normal':
                jobs[0].update(result=1.,error=.1)
            elif change=='two_extremes':
                jobs[1].update(result=-100.,error=90.)
            elif change=='incomplete': jobs.pop()
            elif change=='unequal_points': jobs[0]['npoints_done']=999
            elif change=='unequal_accuracy': jobs[0]['accuracy']=.2
            elif change=='unsplit':
                for j in jobs:j['split']=0
            elif change=='different_channels': jobs[0]['channel']='2'
            elif change=='zero_peer_information':
                for j in jobs[1:]:j['error']=0.
            else:jobs[0]['split_outlier']={}
            with self.subTest(change=change):
                self.assertEqual(find(jobs,10.),[])
        self.assertEqual(find(self.jobs('/unused',count=5),10.),[])

    def test_exclusion_normalizes_rates_components_and_archives_original(self):
        with tempfile.TemporaryDirectory() as root:
            interface = self.interface(root)
            jobs = self.jobs(root)
            other = self.jobs(root,count=1,bad=0)[0]
            other.update(p_dir='P1_other',p_label='P1',split=0,dirname=str(Path(root)/'other'))
            raw = copy.deepcopy(jobs)
            original = Path(jobs[0]['dirname'])/'MADatNLO.HwU'
            original.parent.mkdir(parents=True)
            original.write_text('preserve this raw output\n')
            selected,collected = interface.exclude_fixed_order_split_outliers(jobs,jobs[:]+[other],2)
            self.assertEqual([j['split'] for j in selected],list(range(2,9)))
            self.assertIn(other,collected)
            self.assertEqual(original.read_text(),'preserve this raw output\n')
            self.assertEqual(jobs[0],raw[0])
            self.assertAlmostEqual(sum(j['result'] for j in selected),8.)
            for job in selected:
                self.assertAlmostEqual(job['error'],.1*8./7.)
                self.assertAlmostEqual(job['contribution_results']['total'],8./7.)
                self.assertEqual(job['wgt_mult'],1./8.)
            audit = Path(selected[0]['split_outlier']['audit'])
            report = json.loads(audit.read_text())
            self.assertEqual(report['unfiltered_total']['value'],108.)
            self.assertAlmostEqual(report['filtered_total']['value'],9.)
            self.assertEqual((audit.parent/'P0_test'/'all_G1_1'/'MADatNLO.HwU').read_text(),original.read_text())
            # No repeated trimming/rescaling of a filtered group.
            before = copy.deepcopy(selected)
            interface.exclude_fixed_order_split_outliers(selected,collected,3)
            self.assertEqual(selected,before)
            # Dropping split 1 must not make the entire integration channel disappear.
            with mock.patch.object(interface,'combine_split_order_grids') as grids:
                combined = interface.combine_split_order_run(selected)
            grids.assert_called_once()
            self.assertEqual(len(combined),1)
            self.assertAlmostEqual(combined[0]['result'],8.)
            self.assertNotIn('split_outlier',combined[0])
            self.assertNotIn('split_result_scale',combined[0])
            result = interface.write_res_txt_file(collected,2)
            self.assertAlmostEqual(result['xsect'],9.)
            self.assertAlmostEqual(result['contributions']['total'],9.)
            self.assertAlmostEqual(result['errt'],math.sqrt(7*(.1*8./7.)**2+.1**2))

    def test_saved_filtered_result_is_scaled_once_on_each_reload(self):
        with tempfile.TemporaryDirectory() as root:
            interface = self.interface(root)
            interface.fnlo_multiplicative_enabled = lambda: False
            Path(root,'res_2.dat').write_text('2.0 0.2 1.0 0.1 1 1000 2.0\n')
            job = dict(dirname=root,split_result_scale=8./7.)
            for _ in range(2):
                interface.append_the_results([job],2)
                self.assertAlmostEqual(job['result'],8./7.)
                self.assertAlmostEqual(job['error'],.1*8./7.)

    def test_collection_filters_before_deciding_to_refine(self):
        with tempfile.TemporaryDirectory() as root:
            interface = self.interface(root)
            jobs = self.jobs(root,count=64,bad=34)
            for job in jobs:
                if job['split']!=34:
                    job['error']=job['errorABS']=.01
                    job['contribution_results']['components'][0]['error']=.01
            Path(root,'SubProcesses').mkdir()
            interface.results=mock.Mock()
            with mock.patch.object(interface,'append_the_results'), \
                    mock.patch.object(interface,'make_make_all_html_results',return_value=(64.,.1)), \
                    mock.patch.object(interface,'combine_split_order_grids'), \
                    mock.patch.object(interface,'collect_scale_pdf_info',return_value=[]), \
                    mock.patch.object(interface,'print_summary'):
                new,collected=interface.collect_the_results({},.01,jobs,jobs[:],2,'NLO','all')
            self.assertEqual(new,[])
            self.assertEqual(len(collected),63)
            self.assertAlmostEqual(interface.cross_sect_dict['xsect'],64.)
            self.assertTrue(Path(root,'SubProcesses','job_status.pkl').is_file())

    def test_histograms_and_scale_pdf_use_the_same_retained_normalization(self):
        with tempfile.TemporaryDirectory() as root:
            interface = self.interface(root)
            jobs = self.jobs(root,bad=0)
            jobs[0]['split_result_scale']=8./7.
            fake = mock.Mock()
            fake.poll.return_value=0
            with mock.patch.object(run.misc,'Popen',return_value=fake) as popen:
                interface.combine_plots_HwU(jobs,'output',normalisation=[2.]*8)
            command = popen.call_args[0][0]
            factors = next(x for x in command if x.startswith('--multiply=')).split('=',1)[1]
            self.assertEqual(list(map(float,factors.split(','))),[16./7.]+[2.]*7)
            interface.fnlo_decay_scale_variation_mode = lambda: 'NONE'
            with mock.patch.object(interface,'pdf_scale_from_reweighting',return_value=[]) as collect:
                interface.collect_scale_pdf_info({},jobs)
            self.assertEqual(collect.call_args[0][1],[8./7.]+[1.]*7)

    def test_html_uses_the_same_rate_and_error_normalization(self):
        with tempfile.TemporaryDirectory() as root:
            pdir = Path(root)/'P0_test'
            directory = pdir/'all_G1_2'
            directory.mkdir(parents=True)
            (directory/'results.dat').write_text('2.0 0.2 0.0 0 0 0 0 0 0 1.0\n')
            interface = mock.Mock()
            interface.results.current={'run_name':'test'}
            interface.get_Pdir.return_value=[str(pdir)]
            results = sum_html.collect_result(interface,jobs=[dict(p_dir='P0_test',dirname=str(directory),split_result_scale=8./7.)])
            self.assertAlmostEqual(results.xsec,8./7.)
            self.assertAlmostEqual(results.xerru,.2*8./7.)


if __name__ == '__main__':
    unittest.main()
