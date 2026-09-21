import hashlib
import unittest

from campaign import ROOT
from rng_history import RNG_SOURCE, MINT_SOURCE, validate_history


class TestRNGHistory(unittest.TestCase):
    log='Setting up grids\nRefining results, step 1\nRefining results, step 2\n'

    def fixture(self):
        raw={name:(ROOT/'Template/fNLO'/name).read_bytes() for name in (RNG_SOURCE,MINT_SOURCE)}
        final={}
        for stage in range(3):
            directory='/P0_udx/all_G1'+('_1' if stage else '')
            log=(' with seed 53003\nRanmar stream namespace refinement-v1 stage %d split %d\n'
                 'Ranmar initialization seeds %d 27\nTime spent in Total : 1\n' %
                 (stage,1 if stage else 0,200+stage))
            raw['Events/run/alllogs_%d.html' % stage]=(
                '<HTML><a name=%s></a>\n<font color="red">header</font>\n<PRE>\n%s\n</PRE></HTML>' %
                (directory,log)).encode()
            if stage:
                final['SubProcesses'+directory]=hashlib.sha256(log.encode()).hexdigest()
        return raw,final

    def test_complete_stages_and_final_worker_membership(self):
        raw,final=self.fixture()
        result=validate_history(self.log,53003,raw,final)
        self.assertEqual(result['stage_worker_counts'],{'0':1,'1':1,'2':1})
        self.assertEqual(result['distinct_initialization_pairs'],3)
        self.assertEqual(result['final_workers_checked'],1)
        with self.assertRaisesRegex(ValueError,'last complete'):
            validate_history(self.log,53003,raw,{next(iter(final)):'wrong'})

    def test_missing_source_stage_namespace_and_collision_fail_closed(self):
        for source in (RNG_SOURCE,MINT_SOURCE,'Events/run/alllogs_0.html','Events/run/alllogs_1.html'):
            raw,final=self.fixture()
            del raw[source]
            with self.assertRaises(ValueError):
                validate_history(self.log,53003,raw,final)
        for before,after in [(b'202 27',b'201 27'),(b'stage 2',b'stage 1'),
                             (b'split 1',b'split 2'),(b'53003',b'53004'),
                             (b'refinement-v1',b'legacy'),(b'Time spent in Total',b'incomplete'),
                             (b'202 27',b'40000 27')]:
            raw,final=self.fixture()
            raw['Events/run/alllogs_2.html']=raw['Events/run/alllogs_2.html'].replace(before,after)
            with self.assertRaises(ValueError):
                validate_history(self.log,53003,raw,final)
        raw,final=self.fixture()
        raw[RNG_SOURCE]+=b'changed'
        with self.assertRaisesRegex(ValueError,'source evidence'):
            validate_history(self.log,53003,raw,final)

    def test_imported_grid_and_unrecorded_refinement_not_assumed_safe(self):
        raw,final=self.fixture()
        for log in (self.log.replace('Setting up grids','Imported grids'),
                    self.log+'Refining results, step 3',
                    self.log.replace('step 1','step -1')):
            with self.assertRaises(ValueError):
                validate_history(log,53003,raw,final)


if __name__=='__main__':
    unittest.main()
