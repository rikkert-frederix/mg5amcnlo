"""Regression for fixed-count refinement from a saved fNLO grid."""
import unittest
from pathlib import Path
import tempfile
from unittest.mock import Mock

from madgraph.interface.amcatnlo_run_interface import aMCatNLOCmd
from madgraph.various.banner import RunCardNLO


class TestFixedCountGridResume(unittest.TestCase):
    def interface(self):
        interface = object.__new__(aMCatNLOCmd)
        interface.stop_for_runweb = True
        interface.options = dict(run_mode=2,nb_core=64)
        interface.run_card = RunCardNLO()
        interface.run_card['req_acc_fo'] = -1.
        interface.run_card['fo_job_target_time'] = 1.
        interface.run_card['npoints_fo'] = 65536
        interface.run_card['niters_fo'] = 1
        return interface

    def job(self):
        return dict(p_dir='P0_test',channel='1',dirname='/tmp/P0_test/all_G1',
                    resultABS=1.,errorABS=.01,time_spend=300.,niters=1,npoints=65536,
                    niters_done=1,npoints_done=200,combined=1,accuracy=0.,split=0)

    def test_controller_routes_opt_in_fixed_counts_through_splitter(self):
        interface = self.interface()
        job = self.job()
        interface.append_the_results = Mock()
        interface.write_res_txt_file = Mock(return_value=dict(errt=.1,xsect=1.,erra=.1,xseca=1.))
        interface.make_make_all_html_results = Mock(return_value=(1.,.1))
        interface.results = Mock()
        interface.combine_split_order_run = lambda jobs: jobs
        interface.print_summary = Mock()
        interface.prepare_directories = Mock()
        with tempfile.TemporaryDirectory(prefix='fnlo_fixed_count_controller_') as directory:
            interface.me_dir = directory
            (Path(directory)/'SubProcesses').mkdir()
            jobs,collected = interface.collect_the_results({},-1.,[job],[job],-1,'NLO','all')
        self.assertEqual(len(jobs),64)
        self.assertEqual(len(collected),64)
        self.assertTrue(all(row['accuracy'] == 0. and row['niters'] == 1 for row in jobs))
        self.assertEqual(sum(row['npoints'] for row in jobs),65536)

    def test_fixed_count_split_never_discards_remainder_points(self):
        interface = self.interface()
        for points,iterations,cores in ((65537,1,64),(1001,6,4),(12001,3,64),(12000,3,64)):
            with self.subTest(points=points,iterations=iterations,cores=cores):
                job = self.job()
                job.update(npoints=points,niters=iterations)
                interface.options['nb_core'] = cores
                jobs,collected = interface.split_jobs_fixed_order([job],[job])
                self.assertLessEqual(len(jobs),cores)
                self.assertEqual(sum(row['npoints']*row['niters'] for row in jobs),points*iterations)
                self.assertTrue(all(row['accuracy'] == 0. for row in jobs))

    def test_resume_requires_a_fresh_fixed_count_refinement(self):
        interface = object.__new__(aMCatNLOCmd)
        interface.stop_for_runweb = True
        interface.run_card = RunCardNLO()
        interface.run_card['npoints_fo'] = 12000
        interface.run_card['niters_fo'] = 3
        interface.cross_sect_dict = dict(errt=.1,xsect=10.,erra=.1,xseca=11.)
        job = dict(resultABS=11.,errorABS=.1,niters_done=2,npoints_done=1000,accuracy=.05)
        jobs = interface.update_jobs_to_run(-1.,-1,[job])
        self.assertEqual(len(jobs),1)
        self.assertEqual(jobs[0]['npoints'],12000)
        self.assertEqual(jobs[0]['niters'],3)
        self.assertEqual(jobs[0]['mint_mode'],-1)
        self.assertEqual(jobs[0]['accuracy'],0.)


if __name__ == '__main__':
    unittest.main()
