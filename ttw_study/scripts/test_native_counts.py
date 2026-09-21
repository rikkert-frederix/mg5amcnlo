from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

from campaign import ROOT
from native_current_batches import count_evidence, effective_count


class TestNativeCounts(unittest.TestCase):
    def test_rule_against_actual_compiled_mint_routine(self):
        compiler=shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran is required')
        source=(ROOT/'Template/fNLO/SubProcesses/mint_module.f90').read_text()
        routine=re.search(r'  subroutine initialize_even_random_numbers\b.*?'
                          r'end subroutine initialize_even_random_numbers',source,re.S)[0]
        with tempfile.TemporaryDirectory(prefix='ttw_native_count_') as temporary:
            directory=Path(temporary)
            fixture=directory/'count.f90'
            fixture.write_text('program count\nimplicit none\n'
                'integer :: ncalls0,ndim,ng,k,npg,ncalls\nlogical :: firsttime\n'
                'read(*,*) ncalls0,ndim\ncall initialize_even_random_numbers\n'
                'write(*,*) ncalls\ncontains\n'+routine+'\nend program count\n')
            executable=directory/'count'
            result=subprocess.run([compiler,'-O0','-fcheck=all','-ftrapv',str(fixture),'-o',str(executable)],
                                  cwd=directory,text=True,capture_output=True)
            self.assertEqual(result.returncode,0,result.stdout+result.stderr)
            for requested in (200,1024,16384,38494,40963,65536,262144,1048576):
                for dimensions in (2,7,13,20,63):
                    result=subprocess.run([str(executable)],input='%d %d\n' % (requested,dimensions),
                                          text=True,capture_output=True)
                    self.assertEqual(result.returncode,0,result.stderr)
                    self.assertEqual(int(result.stdout),effective_count(requested,dimensions))
        self.assertEqual(effective_count(38494,13),32768)
        self.assertEqual(effective_count(40963,13),40960)

    def test_deterministic_change_must_precede_samples_and_match_complete_iteration(self):
        log=('about to integrate 13 38494 1\n------- iteration 1\n'
             'Update # PS points: 38494 --> 32768\nRanmar initialization seeds 1 2\n')
        job=dict(npoints=38494,npoints_done=32768)
        self.assertEqual(count_evidence(log,job)['effective'],32768)
        for altered in (log.replace('32768','32000'),log.replace('13 38494 1','13 38494 2'),
                        log+'------- iteration 2\n',log.replace('Update # PS','Missing # PS'),
                        'Ranmar initialization seeds 1 2\n'+log):
            with self.assertRaises(ValueError):
                count_evidence(altered,job)
        with self.assertRaises(ValueError):
            count_evidence(log,dict(job,npoints_done=32767))

    def test_explicit_noop_update_is_preserved_but_not_misread_as_a_draw_change(self):
        log=('about to integrate 13 15494 1\n------- iteration 1\n'
             'Update # PS points: 15494 --> 15494\nRanmar initialization seeds 1 2\n')
        evidence=count_evidence(log,dict(npoints=15494,npoints_done=15494))
        self.assertEqual(evidence['effective'],15494)
        self.assertTrue(evidence['logged_noop_update'])
        for altered in (log.replace('15494 --> 15494','15494 --> 15493'),
                        'about to integrate 13 15494 1\n------- iteration 1\n'
                        'Ranmar initialization seeds 1 2\nUpdate # PS points: 15494 --> 15494\n'):
            with self.assertRaises(ValueError):
                count_evidence(altered,dict(npoints=15494,npoints_done=15494))


if __name__=='__main__':
    unittest.main()
