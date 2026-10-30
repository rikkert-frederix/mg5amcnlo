"""Compile the actual RANMAR module and exercise the actual worker launcher."""
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

from madgraph import MG5DIR


class TestRefinementSeeds(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        compiler=shutil.which('gfortran')
        if not compiler:
            raise unittest.SkipTest('gfortran is required')
        cls.temporary=tempfile.TemporaryDirectory(prefix='fnlo_ranmar_')
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.directory=Path(cls.temporary.name)
        cls.executable=cls.directory/'check'
        root=Path(MG5DIR)
        sources=[root/'Template/fNLO/Source/fnlo_runtime_common.f90',
                 root/'Template/fNLO/Source/ranmar.f90',
                 root/'tests/input_files/fks_decay/ranmar_stage_driver.f90']
        result=subprocess.run([compiler,'-O0','-g','-fcheck=all','-ftrapv',
                               '-o',str(cls.executable),*map(str,sources)],
                              cwd=cls.directory,text=True,capture_output=True)
        if result.returncode:
            raise AssertionError(result.stdout+result.stderr)

    def run_rng(self,stage=None,base=53003,process=0,split=1,configuration=1,moffset=None,group=None):
        with tempfile.TemporaryDirectory(dir=self.directory) as temporary:
            directory=Path(temporary)
            (directory/'randinit').write_text('r=%d\n' % base)
            (directory/'iproc.dat').write_text('%d\n' % process)
            if moffset is not None:
                (directory/'moffset.dat').write_text('%d\n' % moffset)
            env=dict(os.environ)
            env.pop('MG5AMC_REFINEMENT_STAGE',None)
            env.pop('MG5AMC_INTEGRATION_GROUP',None)
            if stage is not None:
                env['MG5AMC_REFINEMENT_STAGE']=str(stage)
            if group is not None:
                env['MG5AMC_INTEGRATION_GROUP']=str(group)
            return subprocess.run([str(self.executable)],cwd=directory,env=env,
                                  input='%d %d\n' % (configuration,split),
                                  text=True,capture_output=True)

    def pair(self,result):
        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
        found=re.findall(r'Ranmar initialization seeds\s+(\d+)\s+(\d+)',result.stdout)
        self.assertEqual(len(found),1,result.stdout)
        return tuple(map(int,found[0]))

    def expected_pair(self,base,process,split,configuration,stage,group=0):
        scaled=base*31300
        ij=1802+configuration+scaled%30081
        kl=9373+scaled//30081+process+3157*split
        ij=(ij-1)%31328+1 if ij>31328 else ij
        kl=(kl-1)%30081+1 if kl>30081 else kl
        rank=(ij*30082+kl+stage*32452843+group*49979687)%(31329*30082)
        return divmod(rank,30082)

    def test_legacy_default_and_stage_zero_are_bit_identical(self):
        default=self.run_rng()
        explicit=self.run_rng(0)
        self.assertEqual(default.stdout,explicit.stdout)
        self.assertEqual(default.stdout,self.run_rng(0,group=0).stdout)
        self.assertEqual(self.pair(default),(28553,7518))
        self.assertEqual(self.pair(self.run_rng(process=1,configuration=0)),(28552,7519))
        self.assertIn('BASE_SEED 53003',default.stdout)
        self.assertEqual(self.run_rng(2).stdout,self.run_rng(2).stdout)
        self.assertEqual(self.run_rng(2).stdout,self.run_rng(2,group=0).stdout)
        self.assertNotEqual(re.findall(r'^DRAW .*',default.stdout,re.M),
                            re.findall(r'^DRAW .*',self.run_rng(2).stdout,re.M))

    def test_training_and_three_refinements_have_distinct_study_worker_pairs(self):
        pairs=set()
        for stage in range(4):
            for process in (0,1):
                for split in range(65):
                    result=self.run_rng(stage,process=process,split=split)
                    pair=self.pair(result)
                    self.assertEqual(pair,self.expected_pair(53003,process,split,1,stage))
                    self.assertNotIn(pair,pairs)
                    pairs.add(pair)
                    self.assertIn('BASE_SEED 53003',result.stdout)
        self.assertEqual(len(pairs),520)
        self.assertEqual(math.gcd(32452843,31329*30082),1)

    def test_group_stage_split_namespaces_are_reproducible_and_distinct(self):
        pairs=set()
        draws=set()
        for group in range(3):
            for stage in (0,1):
                for split in (0,1):
                    for process in (0,1):
                        result=self.run_rng(stage,group=group,split=split,process=process)
                        self.assertEqual(result.stdout,self.run_rng(
                            stage,group=group,split=split,process=process).stdout)
                        pair=self.pair(result)
                        self.assertEqual(pair,self.expected_pair(
                            53003,process,split,1,stage,group))
                        self.assertNotIn(pair,pairs)
                        pairs.add(pair)
                        stream=tuple(re.findall(r'^DRAW .*',result.stdout,re.M))
                        self.assertNotIn(stream,draws)
                        draws.add(stream)
                        self.assertIn('integration-group-v1 group %d' % group,result.stdout)
        self.assertEqual(len(pairs),24)
        self.assertEqual(math.gcd(49979687,31329*30082),1)

    def test_large_integer_offsets_moffset_override_and_stage_bound(self):
        for stage in (0,1,31329*30082-1):
            pair=self.pair(self.run_rng(stage,base=2147483647,process=2147483647,
                                       split=2147483647,configuration=2147483647))
            self.assertEqual(pair,self.expected_pair(2147483647,2147483647,
                                                     2147483647,2147483647,stage))
        self.assertEqual(self.pair(self.run_rng(3,split=999,moffset=12)),
                         self.pair(self.run_rng(3,split=12)))
        for stage in ('','-1','1 2','1,2','1.0','+1','x','9'*80,31329*30082):
            result=self.run_rng(stage)
            self.assertNotEqual(result.returncode,0,stage)
            self.assertTrue('Invalid MG5AMC_REFINEMENT_STAGE' in result.stderr or
                            'stage exceeds its namespace' in result.stderr,result.stderr)
        for seed in (-1,9223372036854775807):
            self.assertNotEqual(self.run_rng(1,base=seed).returncode,0)
        pair=self.pair(self.run_rng(31329*30082-1,group=31329*30082-1,
                                   base=2147483647,process=2147483647,
                                   split=2147483647,configuration=2147483647))
        self.assertEqual(pair,self.expected_pair(2147483647,2147483647,
            2147483647,2147483647,31329*30082-1,31329*30082-1))
        for group in ('','-1','1 2','1,2','1.0','+1','x','9'*80,31329*30082):
            result=self.run_rng(1,group=group)
            self.assertNotEqual(result.returncode,0,group)
            self.assertTrue('Invalid MG5AMC_INTEGRATION_GROUP' in result.stderr or
                            'integration group exceeds its namespace' in result.stderr,
                            result.stderr)

    def test_actual_job_launch_sets_group_stage_split_and_preserves_inputs(self):
        template=(Path(MG5DIR)/'Template/fNLO/SubProcesses/ajob_template').read_text()
        with tempfile.TemporaryDirectory(dir=self.directory) as temporary:
            process=Path(temporary)
            (process/'randinit').write_text('r=53003\n')
            (process/'iproc.dat').write_text('0\n')
            shutil.copy2(self.executable,process/'madevent_mintFO')
            script=process/'ajob1'
            # write_ajob.f90 replaces the marker with this argument-driven loop.
            script.write_text(template.replace('TAGTAGTAGTAGTAGTAGTAG for i in 1 ; do',
                                               'for i in $channel ; do'))
            streams=[]
            for group in (1,2,3):
                for split in (1,2):
                    job=process/('all_G%d_%d' % (group,split))
                    job.mkdir()
                    source='1 %d\n' % split
                    (job/'input_app.txt').write_text(source)
                    for stage in (0,1):
                        result=subprocess.run(['bash',str(script),str(group),'all',
                                               str(split),str(stage)],
                            cwd=process,text=True,capture_output=True,
                            env=dict(os.environ,MG5AMC_REFINEMENT_STAGE='999',
                                     MG5AMC_INTEGRATION_GROUP='999'))
                        self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                        log=(job/('log_MINT%d.txt' % stage)).read_text()
                        self.assertIn('refinement-v1 stage %d split %d' % (stage,split),log)
                        self.assertIn('integration-group-v1 group %d' % (group-1),log)
                        pairs=re.findall(r'Ranmar initialization seeds\s+(\d+)\s+(\d+)',log)
                        self.assertEqual(len(pairs),1,log)
                        pair=tuple(map(int,pairs[0]))
                        self.assertEqual(pair,self.expected_pair(53003,0,split,1,stage,group-1))
                        streams.append(pair)
                        self.assertEqual((job/('input_app_MINT%d.txt' % stage)).read_text(),source)
                        self.assertEqual((process/'randinit').read_text(),'r=53003\n')
            self.assertEqual(len(set(streams)),12)
            invalid=subprocess.run(['bash',str(script),'1','all','1','bad'],
                                   cwd=process,text=True,capture_output=True)
            self.assertNotEqual(invalid.returncode,0)
            self.assertIn('Invalid fixed-order refinement stage',invalid.stderr)
            for group in ('','0','bad','1.0','942438979','9'*80):
                invalid=subprocess.run(['bash',str(script),group,'all','1','0'],
                                       cwd=process,text=True,capture_output=True)
                self.assertNotEqual(invalid.returncode,0)
                self.assertIn('Invalid fixed-order integration group',invalid.stderr)
            job=process/'GF1_1'
            job.mkdir()
            (job/'input_app.txt').write_text('1 1\n')
            result=subprocess.run(['bash',str(script),'1','F','1','0'],
                cwd=process,text=True,capture_output=True,
                env=dict(os.environ,MG5AMC_INTEGRATION_GROUP='999'))
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            self.assertIn('integration-group-v1 group 0',(job/'log_MINT0.txt').read_text())


if __name__=='__main__':
    unittest.main()
