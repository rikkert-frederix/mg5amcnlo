"""Run the actual fNLO job template's MadLoop-resource preamble in bash."""
import io
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest

import madgraph


class TestFNLOJobResources(unittest.TestCase):
    def run_preamble(self,directory,extra_shell=''):
        template = Path(madgraph.MG5DIR)/'Template/fNLO/SubProcesses/ajob_template'
        preamble = template.read_text().split('channel=$1')[0]
        # No sleeping in failure tests; preserve the template's extraction logic.
        script = 'function sleep { :; }\n'+extra_shell+preamble+'\nprintf "READY\\n"\n'
        return subprocess.run(['bash','-c',script],cwd=directory,text=True,capture_output=True)

    def archive(self,directory,names):
        with tarfile.open(directory/'MadLoop5_resources.tar.gz','w:gz') as archive:
            for name in names:
                member = tarfile.TarInfo('MadLoop5_resources/'+name)
                data = ('resource '+name).encode()
                member.size = len(data)
                archive.addfile(member,io.BytesIO(data))

    def test_namespaced_and_legacy_resources(self):
        for names in (['HelConfigs.dat'],['FNLOC1_HelConfigs.dat','FNLOC2_HelConfigs.dat',
                                         'FNLOC2_BornColorFlowCoefs.dat']):
            with self.subTest(names=names),tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                self.archive(directory,names)
                result = self.run_preamble(directory)
                self.assertEqual(result.returncode,0,result.stdout+result.stderr)
                self.assertIn('READY',result.stdout)
                for name in names:
                    self.assertTrue((directory/'MadLoop5_resources'/name).is_file())
                # Already extracted resources must be accepted unchanged.
                marker = directory/'MadLoop5_resources'/names[0]
                marker.write_text('keep this initialized resource')
                result = self.run_preamble(directory)
                self.assertIn('READY',result.stdout)
                self.assertEqual(marker.read_text(),'keep this initialized resource')

    def test_missing_component_does_not_pass_on_one_helicity_file(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            names = ['FNLOC1_HelConfigs.dat','FNLOC2_HelConfigs.dat']
            self.archive(directory,names)
            (directory/'MadLoop5_resources').mkdir()
            (directory/'MadLoop5_resources'/names[0]).write_text('one component only')
            result = self.run_preamble(directory)
            self.assertIn('READY',result.stdout)
            self.assertTrue((directory/'MadLoop5_resources'/names[1]).is_file())

    def test_corrupted_archive_returns_failure(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            (directory/'MadLoop5_resources.tar.gz').write_text('not a tar archive')
            result = self.run_preamble(directory)
            self.assertNotEqual(result.returncode,0)
            self.assertNotIn('READY',result.stdout)

    def test_no_archive_is_needed_for_normal_local_runs(self):
        with tempfile.TemporaryDirectory() as temporary:
            result = self.run_preamble(temporary)
            self.assertEqual(result.returncode,0)
            self.assertIn('READY',result.stdout)

    def test_retry_count_is_numeric_and_reaches_third_attempt(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            self.archive(directory,['FNLOC1_HelConfigs.dat'])
            shell = '''
extract_attempts=0
function tar {
    if [[ "$1" == -xzf ]] ; then
        extract_attempts=$((extract_attempts+1))
        if [[ $extract_attempts -lt 3 ]] ; then return 1; fi
    fi
    command tar "$@"
}
'''
            result = self.run_preamble(directory,shell)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertIn('READY',result.stdout)


if __name__ == '__main__':
    unittest.main()
